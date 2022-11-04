//
// Created by simon on 26/10/22.
//

#include "../branch_unrealistic.h"
#include "../Utils.h"
#include <scip/scipdefplugins.h>
#include "Worker.h"

Worker* Worker::instance = nullptr;
Worker::Worker(unsigned int rank):
    rank(rank),
    directorRank(rank),
    nWorkers(0){
}

void Worker::setInstance(Worker *node) {
    instance = node;
}

Worker *Worker::getInstance() {
    return instance;
}

const unsigned int Worker::getRank() const {
    return rank;
}

bool Worker::isMaster() {
    return !rank;
}

void Worker::work() {
    createScipInstance();
    retrieveInstance();
    std::cout << "Child init done." << std::endl;

    while (true) {
        getWorkersRange();
        SCIP *scip_copy = retrieveNode();
        SCIPsolve(scip_copy);
        int score = getScore(scip_copy);
        //std::cout << "Score: " << score << " - " << SCIPgetSolOrigObj(scip_copy, SCIPgetBestSol(scip_copy)) << std::endl;

        returnScore(score);
        SCIPfree(&scip_copy);
    }
}

void Worker::returnScore(int score) {
    MPI_Send(&score, 1, MPI_INT, directorRank, 1, MPI_COMM_WORLD);
}

void Worker::retrieveInstance() {
    int n; // number of variables in the problem
    int m; // number of constraints in the problem
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    SCIP_Real *objCoef;
    SCIPallocBlockMemoryArray(scipmain, &objCoef, n);

    // Receive variables data
    VarInfo varInfo;
    for(int i = 0; i < n; ++i){
        SCIP_Var* var;

        MPI_Bcast(&varInfo, sizeof(VarInfo), MPI_BYTE, 0, MPI_COMM_WORLD);

        SCIPcreateVar(scipmain, &var, NULL, varInfo.lb, varInfo.ub, varInfo.obj, varInfo.type, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL);
        SCIPaddVar(scipmain, var);
        SCIPreleaseVar(scipmain, &var);
    }

    SCIP_VAR ** vars = SCIPgetVars(scipmain);
    // Receive constraints data
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

    SCIP_Bool success;
    SCIP_Real *vals;
    SCIP_VAR **consVars;
    int *consVarsIndex;
    SCIPallocBlockMemoryArray(scipmain, &vals, n);
    SCIPallocBlockMemoryArray(scipmain, &consVars, n);
    SCIPallocBlockMemoryArray(scipmain, &consVarsIndex, n);



    ConsInfo consInfo;
    for(int j=0; j<m; ++j){
        MPI_Bcast(&consInfo, sizeof(ConsInfo), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(consVarsIndex, consInfo.nVar, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(vals, consInfo.nVar, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for(int i = 0; i < consInfo.nVar; ++i){
            int varNumber = consVarsIndex[i];
            consVars[i] = vars[varNumber];
        }

        SCIP_CONS *cons;
        SCIPcreateConsBasicLinear(
                scipmain,
                &cons,
                (std::to_string(j)+"cons").c_str(),
                consInfo.nVar,
                consVars,
                vals,
                consInfo.lhs,
                consInfo.rhs
        );

        SCIPaddCons(scipmain, cons);
        SCIPreleaseCons(scipmain, &cons);
    }


    SCIPfreeBlockMemoryArray(scipmain, &vals, n);
    SCIPfreeBlockMemoryArray(scipmain, &consVars, n);
    SCIPfreeBlockMemoryArray(scipmain, &consVarsIndex, n);
    SCIPfreeBlockMemoryArray(scipmain, &objCoef, n);

    MPI_Barrier(MPI_COMM_WORLD);
    /*double objLimit;
    MPI_Recv(&objLimit, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    SCIPsetObjlimit(scipmain, objLimit);

    int varToBranchIndex;
    MPI_Recv(&varToBranchIndex, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    SCIP_VAR* varTobranch = vars[varToBranchIndex];

    // TODO: Slave stay in UNREALISTIC branch rule
    Branch_unrealistic *rule = reinterpret_cast<Branch_unrealistic *>(SCIPfindBranchrule(scipmain, "unrealistic"));
    rule->setFirstBranch(varTobranch);*/

}

void Worker::createScipInstance() {
    Utils::create_scip_instance(&scipmain, true);
    //Utils::configure_slave_scip_instance(scipmain);
    SCIPcreateProbBasic(scipmain, ("slaveProb" + std::to_string(rank)).c_str());
}

SCIP *Worker::retrieveNode() {
    double leafTimeLimit;
    MPI_Recv(&leafTimeLimit, 1, MPI_DOUBLE, directorRank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int depth, maxdepth;
    MPI_Recv(&depth, 1, MPI_INT, directorRank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&maxdepth, 1, MPI_INT, directorRank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    int nodeLimit;
    MPI_Recv(&nodeLimit, 1, MPI_INT, directorRank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    int n = SCIPgetNVars(scipmain);
    SCIP_Real *lb;
    SCIP_Real *ub;
    SCIPallocBlockMemoryArray(scipmain, &lb, n);
    SCIPallocBlockMemoryArray(scipmain, &ub, n);
    // get lower bounds
    MPI_Recv(lb, n, MPI_DOUBLE, directorRank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    // get upper bounds
    MPI_Recv(ub, n, MPI_DOUBLE, directorRank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    int firstBrchId;
    MPI_Recv(&firstBrchId, 1, MPI_INT, directorRank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    SCIP* res = createScipInstance(leafTimeLimit, depth, maxdepth, nodeLimit, n, lb, ub, firstBrchId);

    SCIPfreeBlockMemoryArray(scipmain, &lb, n);
    SCIPfreeBlockMemoryArray(scipmain, &ub, n);
    return res;
}

SCIP* Worker::createScipInstance(double leafTimeLimit, int depth, int maxdepth, int nodeLimit, int n, double* lb, double* ub, int firstBrchId) {
    SCIP *scip_copy;
    SCIP_Bool valid;
    SCIPcreate(&scip_copy);
    SCIP_HashMap *varmap;
    SCIPhashmapCreate(&varmap, SCIPblkmem(scip_copy), SCIPgetNVars(scipmain));
    SCIPcopy(scipmain, scip_copy, varmap, nullptr, "ub_subsolver", FALSE, FALSE, FALSE, FALSE, &valid);
    assert(valid == TRUE);

    SCIPcopyParamSettings(scipmain, scip_copy);

    auto *objbranchrule = new Branch_unrealistic(scip_copy, depth, maxdepth, leafTimeLimit);

    SCIPincludeObjBranchrule(scip_copy, objbranchrule, TRUE);
    Utils::configure_scip_instance(scip_copy, true);
    SCIPsetLongintParam(scip_copy, "limits/nodes", nodeLimit);
    SCIPsetIntParam(scip_copy, "display/verblevel",0);
    // don't use heuristics on recursion levels
    SCIPsetHeuristics(scip_copy, SCIP_PARAMSETTING_OFF, TRUE);
    SCIPsetRealParam(scip_copy,"limits/time",600);

    // TODO: Is usefull?
    SCIPmergeVariableStatistics(scipmain, scip_copy, SCIPgetVars(scipmain), SCIPgetVars(scip_copy), SCIPgetNVars(scipmain));

    SCIP_VAR **vars = SCIPgetVars(scip_copy);
    // set lower bounds
    for(int i=0; i<n; ++i){
        SCIPchgVarLb(scip_copy, vars[i], lb[i]);
        SCIPchgVarUb(scip_copy, vars[i], ub[i]);
    }

    objbranchrule->setFirstBranch(vars[firstBrchId]);

    return scip_copy;
}

void Worker::computeScores(SCIP *scip, SCIP_VAR **lpcands, int nlpcands, std::vector<int> &bestcands, int &bestScore,
                           int depth,
                           int maxdepth) {
    int nUnusedWorkers = nWorkers - nlpcands;
    int nWorkersForEach;
    if(nUnusedWorkers > 0) nWorkersForEach = nUnusedWorkers / nlpcands; //each worker can have  n sub workers
    else nWorkersForEach = 0;

    std::cout << "nWorkers: " << nWorkers << std::endl;

    if(nWorkers) {
        for (int c = 0; c < nlpcands; c += nWorkers) {
            // launch the tasks
            for (int s = 0; s < nWorkers; ++s) {
                int i = c + s - 1; // index of the cand
                if (i == nlpcands)break;

                unsigned workerId = startWorkersRange + s;
                unsigned start = 0, end = 0;
                if (nWorkersForEach > 0) {
                    start = startWorkersRange + nlpcands + s * nWorkersForEach;
                    end = start + nWorkersForEach;
                }

                sendWorkersRange(workerId, start, end);
                sendNode(scip, workerId, bestScore == -1 ? -1 : bestScore + 1, lpcands[i], depth, maxdepth);
            }

            // retrieve result
            int score;
            for (int s = 1; s <= nWorkers; ++s) {
                int i = c + s - 1;
                if (i == nlpcands)break;
                MPI_Recv(&score, 1, MPI_INT, s, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (score < bestScore || bestScore == -1) {
                    bestcands.clear();
                    bestcands.push_back(i);
                    bestScore = score;
                } else if (score == bestScore) {
                    bestcands.push_back(i);
                }
            }
        }
    } else{ // no worker available, it's time to do the work by yourself
        for(int i=0; i < nlpcands; ++i){
            SCIP* scip_copy = sendNode(scip, rank, bestScore == -1 ? -1 : bestScore + 1, lpcands[i], depth, maxdepth);
            SCIPsolve(scip_copy);
            int score = getScore(scip_copy);

            std::cout << "Score: " << score << std::endl;
            /*SCIPdebugMsg(scipmain, (std::string(depth, '\t') + std::to_string(i + 1) + "/" + std::to_string(nlpcands) +
                                " (score: " + std::to_string(score) + ") (var: " + SCIPvarGetName(lpcands[i]) +
                                ")\n").c_str());*/

            if (score < bestScore || bestScore == -1) {
                bestcands.clear();
                bestcands.push_back(i);
                bestScore = score;
            } else if (score == bestScore) {
                bestcands.push_back(i);
            }
        }
    }
}

void Worker::broadcastInstance(SCIP *scip) {
    SCIP_VAR** conssVars = SCIPgetVars(scip);
    int n = SCIPgetNVars(scip);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    VarInfo varInfo;
    for(int i=0; i<n; ++i){
        varInfo.type = SCIPvarGetType(conssVars[i]);
        varInfo.lb = SCIPvarGetLbGlobal(conssVars[i]);
        varInfo.ub = SCIPvarGetUbGlobal(conssVars[i]);
        varInfo.obj = SCIPvarGetObj(conssVars[i]);

        MPI_Bcast(&varInfo, sizeof(VarInfo), MPI_BYTE, 0, MPI_COMM_WORLD);
    }

    int m = SCIPgetNConss(scip);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    SCIP_CONS **conss = SCIPgetConss(scip);

    SCIP_Bool success;
    SCIP_Real *vals;
    int *consVarsIndex;
    SCIP_VAR **consVars;
    SCIPallocBlockMemoryArray(scip, &vals, n);
    SCIPallocBlockMemoryArray(scip, &consVars, n);
    SCIPallocBlockMemoryArray(scip, &consVarsIndex, n);

    // for each constraint
    ConsInfo consInfo;
    for(int j=0; j<m; ++j){
        consInfo.lhs = SCIPconsGetLhs(scip, conss[j], &success);
        consInfo.rhs = SCIPconsGetRhs(scip, conss[j], &success);

        // get the number of vars in the constraint
        SCIPgetConsNVars(scip, conss[j], &(consInfo.nVar), &success);
        MPI_Bcast(&consInfo, sizeof(ConsInfo), MPI_BYTE, 0, MPI_COMM_WORLD);

        SCIPgetConsVars(scip, conss[j], consVars, n, &success); // get the conssVars of the constraint
        SCIPgetConsVals(scip, conss[j], vals, n, &success); // get the coefs of the constraint
        for(int i = 0; i < consInfo.nVar; ++i){
            consVarsIndex[i] =  SCIPvarGetProbindex(consVars[i]);
        }

        MPI_Bcast(consVarsIndex, consInfo.nVar, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(vals, consInfo.nVar, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    SCIPfreeBlockMemoryArray(scip, &vals, n);
    SCIPfreeBlockMemoryArray(scip, &consVars, n);
    SCIPfreeBlockMemoryArray(scip, &consVarsIndex, n);

    MPI_Barrier(MPI_COMM_WORLD);
}

SCIP * Worker::sendNode(SCIP *scip, unsigned int workerId, int nodeLimit, SCIP_VAR *varbrch, int depth, int maxdepth) {
    double leafTimeLimit;
    SCIPgetRealParam(scip, "branching/unrealistic/leaftimelimit", &leafTimeLimit);
    if(workerId != rank) {
        MPI_Send(&leafTimeLimit, 1, MPI_DOUBLE, workerId, 1, MPI_COMM_WORLD);
        MPI_Send(&depth, 1, MPI_INT, workerId, 1, MPI_COMM_WORLD);
        MPI_Send(&maxdepth, 1, MPI_INT, workerId, 1, MPI_COMM_WORLD);
        MPI_Send(&nodeLimit, 1, MPI_INT, workerId, 1, MPI_COMM_WORLD);
    }

    int n = SCIPgetNVars(scip);
    SCIP_VAR **vars = SCIPgetVars(scip);
    SCIP_Real *lb;
    SCIP_Real *ub;
    SCIPallocBlockMemoryArray(scip, &ub, n);
    SCIPallocBlockMemoryArray(scip, &lb, n);

    // send lower bounds
    for(int i=0; i<n; ++i){
        lb[i] = SCIPvarGetLbLocal(vars[i]);
    }
    if(workerId != rank)MPI_Send(lb, n, MPI_DOUBLE, workerId, 1, MPI_COMM_WORLD);

    // send upper bounds
    for(int i=0; i<n; ++i){
        ub[i] = SCIPvarGetUbLocal(vars[i]);
    }
    if(workerId != rank)MPI_Send(ub, n, MPI_DOUBLE, workerId, 1, MPI_COMM_WORLD);

    int firstBrchId = SCIPvarGetProbindex(varbrch);
    if(workerId != rank)MPI_Send(&firstBrchId, 1, MPI_INT, workerId, 1, MPI_COMM_WORLD);

    SCIP* res = nullptr;
    if(workerId == rank){
        res = createScipInstance(leafTimeLimit, depth, maxdepth, nodeLimit, n, lb, ub, firstBrchId);
    }

    SCIPfreeBlockMemoryArray(scip, &lb, n);
    SCIPfreeBlockMemoryArray(scip, &ub, n);

    return res;
}

void Worker::setWorkersRange(int start, int end) {
    startWorkersRange = start;
    endWorkersRange = end;
    nWorkers = endWorkersRange - startWorkersRange + 1;
}

void Worker::sendWorkersRange(unsigned int id, unsigned int start, unsigned int end) {
    unsigned data[] = {start, end};
    MPI_Send(data, 2, MPI_UNSIGNED, id, 1, MPI_COMM_WORLD);
}

void Worker::getWorkersRange() {
    unsigned data[2];
    MPI_Recv(data, 2, MPI_UNSIGNED, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    directorRank = status.MPI_SOURCE;
    setWorkersRange(data[0], data[1]);
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
            SCIPdebugMsg(scip, ("Time limit: " + std::to_string(SCIPgetTotalTime(scip)) + "\n").c_str());
            break;
    }
    return score;
}

void Worker::setScipInstance(SCIP *scip) {
    scipmain = scip;
}