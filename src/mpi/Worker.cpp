//
// Created by simon on 26/10/22.
//
#define SCIP_DEBUG

#include "../branch_unrealistic.h"
#include "../Utils.h"
#include <iostream>
#include <scip/scipdefplugins.h>
#include <algorithm>
#include <scip/cons.h>
#include "Worker.h"

Worker* Worker::instance = nullptr;
Worker::Worker(unsigned int rank):
    rank(rank),
    directorRank(rank),
    nWorkers(0),
    workers(nullptr){
}

void Worker::setInstance(Worker *node) {
    instance = node;
}

Worker *Worker::getInstance() {
    return instance;
}

bool Worker::isMaster() {
    return !rank;
}

void Worker::work() {
    createScipInstance();
    retrieveInstance();
    while (!isFinished()) {
        getWorkersRange();
        SCIP *scip_copy = retrieveNode();
        SCIPsolve(scip_copy);

        int score = getScore(scip_copy);
        returnScore(score);
        SCIPfree(&scip_copy);
    }
    SCIPfree(&scipmain);
}

void Worker::returnScore(int score) {
    MPI_Send(&score, 1, MPI_INT, directorRank, SCORE_FLAG, MPI_COMM_WORLD);
}

void Worker::retrieveInstance() {
    int len;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    char* name =  new char[len+1];
    MPI_Bcast(name, len, MPI_CHAR, 0, MPI_COMM_WORLD);
    name[len] = '\0';

    SCIPreadProb(scipmain, name, NULL);
    delete[] name;
    MPI_Barrier(MPI_COMM_WORLD);
}

void Worker::createScipInstance() {
    Utils::create_scip_instance(&scipmain, true);
    //Utils::configure_slave_scip_instance(scipmain);
    //SCIPcreateProbBasic(scipmain, ("workerProb" + std::to_string(rank)).c_str());
}

SCIP *Worker::retrieveNode() {
    double leafTimeLimit;
    MPI_Recv(&leafTimeLimit, 1, MPI_DOUBLE, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);
    int depth, maxdepth;
    MPI_Recv(&depth, 1, MPI_INT, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&maxdepth, 1, MPI_INT, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);

    int nodeLimit;
    MPI_Recv(&nodeLimit, 1, MPI_INT, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);

    int n = SCIPgetNVars(scipmain);
    SCIP_Real *lb=new SCIP_Real[n];
    SCIP_Real *ub=new SCIP_Real[n];;
    SCIP_Real *bestSolVals=new SCIP_Real[n];;
    /*SCIPallocBufferArray(scipmain, &lb, n);
    SCIPallocBufferArray(scipmain, &ub, n);
    SCIPallocBufferArray(scipmain, &bestSolVals, n);*/
    // get lower bounds
    MPI_Recv(lb, n, MPI_DOUBLE, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);

    // get upper bounds
    MPI_Recv(ub, n, MPI_DOUBLE, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);

    int firstBrchId;
    MPI_Recv(&firstBrchId, 1, MPI_INT, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);

    double objlimit, left, right;
    MPI_Recv(&objlimit, 1, MPI_DOUBLE, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);
    MPI_Recv(bestSolVals, n, MPI_DOUBLE, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&left, 1, MPI_DOUBLE, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&right, 1, MPI_DOUBLE, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);

    int branchingMaxDepth;
    MPI_Recv(&branchingMaxDepth, 1, MPI_INT, directorRank, NODE_INFO_FLAG, MPI_COMM_WORLD, &status);

    SCIP* res = createScipInstance(leafTimeLimit, depth, maxdepth, nodeLimit, n, lb, ub, firstBrchId, objlimit,
                                   bestSolVals, left, right, branchingMaxDepth);


    /*SCIPfreeBufferArray(scipmain, &bestSolVals);
    SCIPfreeBufferArray(scipmain, &ub);
    SCIPfreeBufferArray(scipmain, &lb);*/
    delete[] bestSolVals;
    delete[] ub;
    delete[] lb;
    return res;
}

SCIP * Worker::createScipInstance(double leafTimeLimit, int depth, int maxdepth, int nodeLimit, int n, double *lb,
                                  double *ub,
                                  int firstBrchId, double objlimit, double *bestSolvals, double left, double right,
                                  int branchingMaxDepth) {
    SCIP *scip_copy=nullptr;
    SCIP_Bool valid;
    SCIPcreate(&scip_copy);
    SCIPenableDebugSol(scip_copy);

    SCIPcopy(scipmain, scip_copy, nullptr, nullptr, "sub_", TRUE, FALSE, TRUE, TRUE, &valid);
    assert(valid == TRUE);

    //SCIPcopyParamSettings(scipmain, scip_copy);

    auto *objbranchrule = new Branch_unrealistic(scip_copy, depth, maxdepth, leafTimeLimit);

    SCIPincludeObjBranchrule(scip_copy, objbranchrule, TRUE);
    Utils::configure_scip_instance(scip_copy, true);
    SCIPsetLongintParam(scip_copy, "limits/nodes", nodeLimit);

    SCIPsetIntParam(scip_copy, "display/verblevel", 0);

    if(depth == maxdepth){
        SCIPsetIntParam(scip_copy, "branching/unrealistic/maxdepth", 0);
    }
    // don't use heuristics on recursion levels
    SCIPsetHeuristics(scip_copy, SCIP_PARAMSETTING_OFF, TRUE);
    //SCIPsetRealParam(scip_copy,"limits/time",1e+20);

    SCIP_VAR **vars = new SCIP_VAR*[n];
    memcpy(vars, SCIPgetVars(scip_copy), n*sizeof(SCIP_VAR*));
    sortVarsArray(vars, n);

    // set lower bounds
    for(int i=0; i<n; ++i){
        //SCIPvarSetRemovable(vars[i], FALSE);
        /*if(lb[i] == ub[i]){
            SCIP_Bool fixed, infeasible;
            SCIPfixVar(scip_copy, vars[i], ub[i], &infeasible, &fixed);
            assert(fixed);*/
        //} else {
            SCIPchgVarLbGlobal(scip_copy, vars[i], lb[i]);
            SCIPchgVarUbGlobal(scip_copy, vars[i], ub[i]);
        //}
    }

    assert(firstBrchId < n);
    objbranchrule->setFirstBranch(vars[firstBrchId], left, right);
    SCIPsetObjlimit(scip_copy, objlimit);

    SCIP_Sol* sol;
    SCIPcreateSol(scip_copy, &sol, NULL);
    SCIPsetSolVals(scip_copy, sol, n, vars, bestSolvals);
    SCIP_Bool stored;
    SCIPaddSol(scip_copy, sol, &stored);
    assert(stored);
    SCIPfreeSol(scip_copy, &sol);


    //SCIPsetIntParam(scip_copy, "branching/unrealistic/maxdepth", branchingMaxDepth);
    delete[] vars;

    return scip_copy;
}

void Worker::computeScores(SCIP *scip, SCIP_VAR **lpcands, int nlpcands, std::vector<int> &bestcands, int &bestScore,
                           int depth,
                           int maxdepth, double leafTimeLimit, bool noNodeLimitation, int *varScores) {

    // Check if the worker range can be updated
    /*int flag=0;
    MPI_Iprobe( MPI_ANY_SOURCE, WORKER_RANGE_FLAG, MPI_COMM_WORLD, &flag, &status );
    if(flag){
        getWorkersRange();
    }*/

    SCIP_Real objlimit;
    // get values of the best sol
    SCIP_Sol* bestSol = SCIPgetBestSol(scip);
    //objlimit = SCIPgetObjlimit(scip);

    if(bestSol == NULL) {
        objlimit = SCIPgetObjlimit(scip);
    } else{
        objlimit = SCIPsolGetOrigObj(bestSol);
    }

    //std::cout << rank << " working for " << directorRank << " with " << nWorkers << " workers" << std::endl;

    unsigned minNumberWorker;
    if(depth==maxdepth){
        minNumberWorker = 1;
    } else{
        minNumberWorker = std::max(MIN_NUMBER_WORKER, (int)std::ceil(double(nWorkers)/nlpcands));
    }
    if(nWorkers) {
        int *workerMap = new int[nWorkers]; // workermap[i]=j -> The i th worker works on candidate j
        std::fill(workerMap, workerMap+nWorkers, -1);
        int nActiveWorkers = 0;

        for (int i = 0; i < nlpcands; ++i) {
            unsigned workerRank;
            unsigned n=minNumberWorker;
            assert(n);
            unsigned* subWorkers = findAvailableWorkers(n, i, workerMap);
            int flag=0;
            MPI_Iprobe( MPI_ANY_SOURCE, SCORE_FLAG, MPI_COMM_WORLD, &flag, &status );
            if(flag || n==0) { // an old worker is available or no more available worker
                workerRank = extractScore(lpcands, bestcands, depth, workerMap, bestScore, nullptr);
                updateWork(workerRank, i, workerMap);
                n = -1;
            } else{
                nActiveWorkers++;
                workerRank = subWorkers[0];
            }

            if(n==-1){
                std::cout << rank << " re-sending work to " << workerRank << std::endl;
            } else{
                std::cout << rank << " sending work to " << workerRank << ", n=" << n << "( ";
                for(int h=0; h<n; ++h){
                    std::cout << subWorkers[h] << ", ";
                }
                std::cout << ")" << std::endl;
            }

            short instruction = CONTINUE_WORKING_FLAG;
            MPI_Send(&instruction, 1, MPI_SHORT, workerRank, INSTRUCTION_FLAG, MPI_COMM_WORLD);

            if(n!=-1){
                sendWorkersRange(workerRank, n==1?nullptr:subWorkers+1, n-1);
            } else{
                sendWorkersRange(workerRank, nullptr, -1);
            }

            int nodelimit = -1;
            if(!noNodeLimitation && bestScore < INT_MAX && bestScore != -1){
                nodelimit = bestScore + 1;
            }
            sendNode(scip, workerRank, nodelimit, lpcands[i], depth, maxdepth, objlimit, leafTimeLimit);
            delete[] subWorkers;
        }

        for(int j=0; j<nActiveWorkers; ++j){
            extractScore(lpcands, bestcands, depth, workerMap, bestScore, varScores);
        }

        delete[] workerMap;
    } else{ // no worker available, it's time to do the work by yourself
        for(int i=0; i < nlpcands; ++i){
            int nodelimit = -1;
            if(!noNodeLimitation && bestScore < INT_MAX && bestScore != -1){
                nodelimit = bestScore + 1;
            }
            SCIP* scip_copy = sendNode(scip, rank, nodelimit, lpcands[i], depth, maxdepth,
                                       objlimit, leafTimeLimit);
            SCIPsolve(scip_copy);
            int score = getScore(scip_copy);

            if(!(depth-1))std::cout << "Score of " << SCIPvarGetName(lpcands[i]) <<": " << score << std::endl;

            if (score < bestScore || bestScore == -1) {
                bestcands.clear();
                bestcands.push_back(i);
                bestScore = score;
            } else if (score == bestScore) {
                bestcands.push_back(i);
            }
            if(varScores)varScores[i] = score;
            SCIPfree(&scip_copy);
        }
    }
}

unsigned int
Worker::extractScore(SCIP_VAR *const *lpcands, std::vector<int> &bestcands, int depth, const int *workerMap,
                     int &bestScore, int *varScores) const {
    int score, cand;
    unsigned workerRank;
    MPI_Status status;

    MPI_Recv(&score, 1, MPI_INT, MPI_ANY_SOURCE, SCORE_FLAG, MPI_COMM_WORLD, &status);
    workerRank = status.MPI_SOURCE;

    unsigned workerIndex = findWorkerIndex(workerRank);
    cand = workerMap[workerIndex];

    if(!(depth-1))std::cout << "Score of " << SCIPvarGetName(lpcands[cand]) <<": " << score << std::endl;
    if (score < bestScore || bestScore == -1) {
        bestcands.clear();
        bestcands.push_back(cand);
        bestScore = score;
    } else if (score == bestScore) {
        bestcands.push_back(cand);
    }

    if(varScores)varScores[cand] = score;

    return workerRank;
}

void Worker::broadcastInstance(const char *name) {
    // broadcast the instance name
    int len = strlen(name);

    this->name =  new char[len+1];
    memcpy(this->name, name, len*sizeof(char));
    this->name[len] = '\0';

    MPI_Bcast(&len, 1, MPI_INT, rank, MPI_COMM_WORLD);
    MPI_Bcast((void*)name, len, MPI_CHAR, rank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

SCIP * Worker::sendNode(SCIP *scip, unsigned int workerId, int nodeLimit, SCIP_VAR *varbrch, int depth, int maxdepth,
                        double objlimit, double leafTimeLimit) {
    if(workerId != rank) {
        MPI_Send(&leafTimeLimit, 1, MPI_DOUBLE, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
        MPI_Send(&depth, 1, MPI_INT, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
        MPI_Send(&maxdepth, 1, MPI_INT, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
        MPI_Send(&nodeLimit, 1, MPI_INT, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
    }

    int n = SCIPgetNVars(scip);
    SCIP_VAR **vars = new SCIP_VAR*[n];
    memcpy(vars, SCIPgetVars(scip), n*sizeof(SCIP_VAR*));
    sortVarsArray(vars, n);
    SCIP_Real *lb=new SCIP_Real[n];;
    SCIP_Real *ub=new SCIP_Real[n];;
    SCIP_Real *bestSolVals=new SCIP_Real[n];;
    /*SCIPallocBufferArray(scip, &ub, n);
    SCIPallocBufferArray(scip, &lb, n);
    SCIPallocBufferArray(scip, &bestSolVals, n);*/

    // send lower bounds
    for(int i=0; i<n; ++i){
        SCIP_VAR *var=vars[i];
        lb[i] = SCIPvarGetLbGlobal(var); // TODO: Or SCIPgetVarLbLocal?
        ub[i] = SCIPvarGetUbGlobal(var);
    }
    if(workerId != rank)MPI_Send(lb, n, MPI_DOUBLE, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
    if(workerId != rank)MPI_Send(ub, n, MPI_DOUBLE, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);

    int firstBrchId=-1;
    for(int i=0; i<n; ++i){
        if(vars[i] == varbrch){
            firstBrchId = i;
            break;
        }
    }
    //= SCIPvarGetProbindex(varbrch);
    assert(firstBrchId >= 0);

    SCIP_Sol* bestSol = SCIPgetBestSol(scip);
    SCIPgetSolVals(scip, bestSol, n, vars, bestSolVals);

    double varVal = SCIPvarGetLPSol(varbrch);
    double left = std::floor(varVal);
    double right = std::ceil(varVal);

    int branchingMaxDepth;
    SCIPgetIntParam(scip, "branching/unrealistic/maxdepth", &branchingMaxDepth);

    if(workerId != rank){
        MPI_Send(&firstBrchId, 1, MPI_INT, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
        MPI_Send(&objlimit, 1, MPI_DOUBLE, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
        MPI_Send(bestSolVals, n, MPI_DOUBLE, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
        MPI_Send(&left, 1, MPI_DOUBLE, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
        MPI_Send(&right, 1, MPI_DOUBLE, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
        MPI_Send(&branchingMaxDepth, 1, MPI_INT, workerId, NODE_INFO_FLAG, MPI_COMM_WORLD);
    }


    SCIP* res = nullptr;
    if(workerId == rank){
        res = createScipInstance(leafTimeLimit, depth, maxdepth, nodeLimit, n, lb, ub, firstBrchId, objlimit,
                                 bestSolVals, left, right, branchingMaxDepth);
    }

    /*SCIPfreeBufferArray(scip, &bestSolVals);
    SCIPfreeBufferArray(scip, &lb);
    SCIPfreeBufferArray(scip, &ub);*/
    delete[] bestSolVals;
    delete[] lb;
    delete[] ub;
    delete[] vars;

    return res;
}

void Worker::setWorkersRange(int start, int end) {
    nWorkers = end-start;
    workers = new unsigned[nWorkers];
    for(int i=0; i<nWorkers; i++){
        workers[i] = start+i;
    }
}

void Worker::sendWorkersRange(unsigned int workerRank, unsigned int *workersToSend, unsigned int n) {
    MPI_Send(&n, 1, MPI_UNSIGNED, workerRank, WORKER_RANGE_FLAG, MPI_COMM_WORLD);
    if(n!=0 && n!= -1){
        assert(workersToSend!=nullptr);
        MPI_Send(workersToSend, n, MPI_UNSIGNED, workerRank, WORKER_RANGE_FLAG, MPI_COMM_WORLD);
    }
}

void Worker::getWorkersRange() {
    unsigned n;
    MPI_Recv(&n, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, WORKER_RANGE_FLAG, MPI_COMM_WORLD, &status);
    directorRank = status.MPI_SOURCE;

    if(n != -1){ // if n==-1: The worker keeps its subWorkers
        delete[] workers;
        workers = nullptr;
        nWorkers = n;
        if(n!=0){
            workers = new unsigned[nWorkers];
            MPI_Recv(workers, nWorkers, MPI_UNSIGNED, directorRank, WORKER_RANGE_FLAG, MPI_COMM_WORLD, &status);
        }
    }
}

int Worker::getScore(SCIP *scip) {
    int score = INT_MAX;
    SCIP_STATUS status = SCIPgetStatus(scip);
    switch (status){
        case SCIP_STATUS_NODELIMIT:
        case SCIP_STATUS_OPTIMAL:
        case SCIP_STATUS_INFEASIBLE:
            score = SCIPgetNNodes(scip);
            break;
        case SCIP_STATUS_TIMELIMIT:
            break;
    }
    return score;
}

void Worker::setScipInstance(SCIP *scip) {
    createScipInstance();
    SCIPreadProb(scipmain, SCIPgetProbName(scip), NULL);

    broadcastInstance(SCIPgetProbName(scip));
}

void Worker::broadcastEnd() {
    short instruction = END_SOLVING_FLAG;
    for(int i=0; i<nWorkers; ++i){
        MPI_Send(&instruction, 1, MPI_SHORT, workers[i], INSTRUCTION_FLAG, MPI_COMM_WORLD);
    }
    SCIPfree(&scipmain);
}

bool Worker::isFinished() {
    short instruction;
    MPI_Recv(&instruction, 1, MPI_SHORT, MPI_ANY_SOURCE, INSTRUCTION_FLAG, MPI_COMM_WORLD, &status);
    return instruction==END_SOLVING_FLAG;
}

Worker::~Worker() {
    delete[] workers;
    delete[] name;
}

unsigned *Worker::findAvailableWorkers(unsigned int &n, int task, int *workerMap) {
    unsigned* res = new unsigned[n];
    unsigned nFound=0;
    for(int i=0; i<nWorkers; ++i){
        if(workerMap[i] == -1){
            workerMap[i]=task;
            res[nFound++]=workers[i];
            if(nFound==n)break;
        }
    }

    if(nFound<n){
        n = nFound;
    }
    return res;
}

void Worker::freeWorkers(int task, int *workerMap) {
    for(int i=0; i<nWorkers; ++i){
        if(workerMap[i] == task){
            workerMap[i] = -1;
        }
    }
}

unsigned Worker::findWorkerIndex(unsigned int workerRank) const {
    for(int i=0; i<nWorkers; ++i){
        if(workers[i] == workerRank){
            return i;
        }
    }
    return nWorkers;
}

void Worker::updateWork(unsigned int workerRank, int task, int *workerMap) {
    for(int i=0; i<nWorkers; ++i){
        if(workers[i] == workerRank){
            workerMap[i] = task;
        }
    }
}

void Worker::sortVarsArray(SCIP_VAR **vars, int n) {
    struct {
        bool operator()(SCIP_VAR* a, SCIP_VAR* b) const { return SCIPvarGetIndex(a) < SCIPvarGetIndex(b); }
    } customLess;
    std::sort(vars, vars+n, customLess);
}
