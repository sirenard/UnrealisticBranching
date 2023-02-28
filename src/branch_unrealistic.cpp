//
// Created by simon on 14/08/22.
//
#define SCIP_DEBUG

#include <scip/var.h>
#include <sstream>
#include "branch_unrealistic.h"
#include "Utils.h"
#include "mpi/Worker.h"

#define 	BRANCHRULE_NAME   "unrealistic"
#define 	BRANCHRULE_DESC   "unrealistic branching"
#define 	BRANCHRULE_PRIORITY   200
#define 	BRANCHRULE_MAXDEPTH   -1
#define 	BRANCHRULE_MAXBOUNDDIST   1.0

Branch_unrealistic::Branch_unrealistic(SCIP *scip, int maxdepth, double leafTimeLimit) : ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                                                                                            BRANCHRULE_MAXBOUNDDIST), depth(0), maxdepth(maxdepth), leafTimeLimit(leafTimeLimit), branchingHistory(new std::vector<int>), branchingHistoryValues(new std::vector<double>), branching_count(-1){}


SCIP_RETCODE Branch_unrealistic::branchCopycat(SCIP *scip, SCIP_RESULT *result) {
    SCIP_Var *varbranch = nullptr;
    SCIP_Var** vars = SCIPgetVars(scip);
    int n = SCIPgetNVars(scip);
    for(int i=0;i<n;++i){
        if(SCIPvarGetProbindex(vars[i]) == branchingHistory->at(branching_count)){
            varbranch = vars[i];
            break;
        }
    }

    double value = branchingHistoryValues->at(branching_count);

    if(branching_count == branchingHistory->size()-1){
        branching_count = -1;

        if(depth==maxdepth){
            SCIPsetIntParam(scip, "branching/unrealistic/priority", 0);
            SCIPsetRealParam(scip, "limits/time", leafTimeLimit);
        }
    } else{
        branching_count++;
    }

    *result = SCIP_BRANCHED;
    SCIP_CALL(SCIPbranchVarVal(scip, varbranch, value, NULL, NULL, NULL));

    return SCIP_OKAY;
}

SCIP_RETCODE Branch_unrealistic::branchUnrealistic(SCIP *scip, SCIP_RESULT *result) {
    SCIP_VAR** lpcands;
    SCIP_Real * lpcandsfrac;
    int nlpcands; // number of candidates
    std::vector<int> bestcands;
    int bestScore = -1;

    // get branching candidates
    SCIP_CALL(SCIPgetLPBranchCands(scip, &lpcands, nullptr, &lpcandsfrac, nullptr, &nlpcands, nullptr));
    assert(nlpcands > 0);

    for (int i = 0; i < nlpcands; ++i)lpcandsfrac[i] = SCIPvarGetLPSol(lpcands[i]);

    bool random = false;
    if(dataWriter != nullptr && depth == 0){
        //random = (double) rand() / double(RAND_MAX) < k;
        random=SCIPgetNNodes(scip) <= k;
    }
    int bestcand;
    int *varScores = nullptr; // store every variable's realNnodes

    if(!random) {
        if (dataWriter != nullptr && depth == 0)varScores = new int[nlpcands];
        Worker *worker = Worker::getInstance();
        worker->computeScores(scip, lpcands, nlpcands, bestcands, bestScore, depth + 1, maxdepth, leafTimeLimit,
                              dataWriter != nullptr && depth == 0, varScores);
        bestcand = bestcands[rand() % bestcands.size()];
    } else{
        bestcand = rand()%nlpcands;
    }


    if (!depth) {
        SCIPdebugMsg(scip, ("Var to branch: " + std::to_string(bestcand + 1) + "; " +
                            std::string(SCIPvarGetName(lpcands[bestcand])) + "; Score: " +
                            std::to_string(bestScore) + "\n").c_str());
    }

    branchingHistory->push_back(SCIPvarGetProbindex(lpcands[bestcand]));
    branchingHistoryValues->push_back(lpcandsfrac[bestcand]);
    SCIP_Node* children[2];
    SCIP_CALL(SCIPbranchVar(scip, lpcands[bestcand], &children[0], NULL, &children[1]));
    *result = SCIP_BRANCHED;


    if(dataWriter && depth==0 && !random) {
        // the score must be: how many nodes needed from the current node. Thus remove from each score the current
        // number of nodes
        for(int i=0; i<nlpcands; ++i){
            if(varScores[i] == INT_MAX)continue;
            varScores[i] -= SCIPgetNNodes(scip)-1; // rempve the number of already used nodes
        }
        dataWriter->addNode(scip, children, nlpcands, varScores, lpcands, bestcand, scoreMethod, alpha);
    }

    delete[] varScores;
    return SCIP_OKAY;
}

SCIP_DECL_BRANCHEXECLP(Branch_unrealistic::scip_execlp){
    // an instruction already exists for the first branching (given by the parent)
    if(branching_count >= 0){
        return branchCopycat(scip, result);
    } else{
        return branchUnrealistic(scip, result);
    }
}


DatasetWriter* Branch_unrealistic::dataWriter = nullptr;
void Branch_unrealistic::setDataWriter(DatasetWriter *dataWriter) {
    Branch_unrealistic::dataWriter = dataWriter;
}

int *Branch_unrealistic::getMaxDepthPtr() {
    return &maxdepth;
}

double *Branch_unrealistic::getLeafTimeLimitPtr() {
    return &leafTimeLimit;
}


SCIP_DECL_BRANCHEXIT(Branch_unrealistic::scip_exit){
    Worker *worker = Worker::getInstance();
    if(worker->isMaster() && depth==0){
        worker->broadcastEnd();
    }
    return SCIP_OKAY;
}

void Branch_unrealistic::setDepth(int depth) {
    Branch_unrealistic::depth = depth;
}

Branch_unrealistic::~Branch_unrealistic() {
    delete branchingHistory;
    delete branchingHistoryValues;
}

void Branch_unrealistic::fillBranchHistory(int *history, double *values, int size) {
    branching_count = 0;
    for(int i=0; i<size; ++i){
        branchingHistory->push_back(history[i]);
        branchingHistoryValues->push_back(values[i]);
    }
}

void Branch_unrealistic::setLeafTimeLimit(double leafTimeLimit) {
    Branch_unrealistic::leafTimeLimit = leafTimeLimit;
}

std::vector<int> *Branch_unrealistic::getHistory() {
    return branchingHistory;
}

std::vector<double> *Branch_unrealistic::getBranchingHistoryValues() const {
    return branchingHistoryValues;
}

char *Branch_unrealistic::getScoreMethodPtr() {
    return &scoreMethod;
}

double *Branch_unrealistic::getAlphaPtr() {
    return &alpha;
}

int * Branch_unrealistic::getKPtr() {
    return &k;
}
