//
// Created by simon on 14/08/22.
//
#define SCIP_DEBUG

#include <scip/var.h>
#include <sstream>
#include "branch_unrealistic.h"
#include "Utils.h"
#include "mpi/Worker.h"
#include "EventhdlrUpdateFeatures.h"

#define 	BRANCHRULE_NAME   "unrealistic"
#define 	BRANCHRULE_DESC   "unrealistic branching"
#define 	BRANCHRULE_PRIORITY   200
#define 	BRANCHRULE_MAXDEPTH   -1
#define 	BRANCHRULE_MAXBOUNDDIST   1.0

Branch_unrealistic::Branch_unrealistic(SCIP *scip, int maxdepth, double leafTimeLimit) : ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                                                                                            BRANCHRULE_MAXBOUNDDIST), depth(0), maxdepth(maxdepth), leafTimeLimit(leafTimeLimit), branchingHistory(new BranchingHistory()), branching_count(0){}


SCIP_RETCODE Branch_unrealistic::branchCopycat(SCIP *scip, SCIP_RESULT *result) {
    SCIP_Var *varbranch = nullptr;
    SCIP_Var** vars = SCIPgetVars(scip);
    int n = SCIPgetNVars(scip);

    if(branchingHistory->at(branching_count).varIndex == -1){ // exploration must be repeated
        *result = SCIP_DIDNOTRUN;
    } else {
        for (int i = 0; i < n; ++i) { // search for the good variable
            if (SCIPvarGetProbindex(vars[i]) == branchingHistory->at(branching_count).varIndex) {
                varbranch = vars[i];
                break;
            }
        }
        assert(varbranch != nullptr);
    }


    double value = branchingHistory->at(branching_count).varValue;

    branching_count++;

    if(varbranch) {
        *result = SCIP_BRANCHED;
        SCIP_CALL(SCIPbranchVarVal(scip, varbranch, value, NULL, NULL, NULL));
    }

    return SCIP_OKAY;
}

SCIP_RETCODE Branch_unrealistic::branchUnrealistic(SCIP *scip, SCIP_RESULT *result) {
    SCIP_VAR** lpcands;
    SCIP_Real * lpcandsfrac;
    int nlpcands; // number of candidates
    std::vector<int> bestcands;
    int bestScore = -1;

    // true if current process is generating a dataset
    // false if UB is simply used
    bool generatingData =  dataWriter != nullptr && depth == 0;

    // get branching candidates
    SCIP_CALL(SCIPgetLPBranchCands(scip, &lpcands, nullptr, &lpcandsfrac, nullptr, &nlpcands, nullptr));
    assert(nlpcands > 0);

    for (int i = 0; i < nlpcands; ++i)lpcandsfrac[i] = SCIPvarGetLPSol(lpcands[i]);

    int bestcand;
    int *varScores =new int[nlpcands]; // store every variable's realNnodes


    bool exploration = generatingData && (double) rand() / double(RAND_MAX) < epsilon;

    //if(!exploration) {
        Worker *worker = Worker::getInstance();
        worker->computeScores(scip, lpcands, nlpcands, bestcands, bestScore, depth + 1, maxdepth, leafTimeLimit,
                              dataWriter != nullptr && depth == 0, varScores);
        bestcand = bestcands[rand() % bestcands.size()];
    //}

    // If generating dataset AND with a prob epsilon, exploration by using another scheme
    if(exploration || varScores[bestcand]==INT_MAX){
        *result = SCIP_DIDNOTRUN;
        EventhdlrUpdateFeatures* eventHdlr = dynamic_cast<EventhdlrUpdateFeatures *>(SCIPfindObjEventhdlr(scip, EVENT_HDLR_UPDATE_FEATURES_NAME));
        eventHdlr->informNextOneIsExploration();
        std::cout << "EXPLORE" << std::endl;
    } else { //otherwise, we branch on the best candidate found
        SCIP_CALL(SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL));
        *result = SCIP_BRANCHED;

        if (!depth) {
            SCIPdebugMsg(scip, ("Var to branch: " + std::to_string(bestcand + 1) + "; " +
                                std::string(SCIPvarGetName(lpcands[bestcand])) + "; Score: " +
                                std::to_string(bestScore) + "\n").c_str());
        }
    }


    if(generatingData) {
        // the score must be: how many nodes needed from the current node. Thus remove from each score the current
        // number of nodes
        for (int i = 0; i < nlpcands; ++i) {
            if (varScores[i] == INT_MAX)continue;
            varScores[i] -= SCIPgetNNodes(scip) - 1; // rempve the number of already used nodes
        }
        dataWriter->addNode(scip, nlpcands, varScores, lpcands, scoreMethod, alpha);
    }

    delete[] varScores;
    return SCIP_OKAY;
}

SCIP_DECL_BRANCHEXECLP(Branch_unrealistic::scip_execlp){
    if(branching_count == branchingHistory->size()){ // must turn off copycat branching
        branching_count = -1;

        EventhdlrUpdateFeatures* eventHdlr = dynamic_cast<EventhdlrUpdateFeatures *>(SCIPfindObjEventhdlr(scip, EVENT_HDLR_UPDATE_FEATURES_NAME));
        if(depth==maxdepth){
            Utils::congigure_scip_end_recursion(scip, leafTimeLimit);

            // disable history recording
            eventHdlr->setHistory(nullptr);
        } else{
            eventHdlr->setHistory(branchingHistory);
        }
    }

    // an instruction already exists for the first branching (given by the parent)
    if(branching_count >= 0){
        return branchCopycat(scip, result);
    } else if(depth < maxdepth){
        return branchUnrealistic(scip, result);
    } else{
        *result = SCIP_DIDNOTRUN;
        return SCIP_OKAY;
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
    branchingHistory->clear();
    return SCIP_OKAY;
}

void Branch_unrealistic::setDepth(int depth) {
    Branch_unrealistic::depth = depth;
}

Branch_unrealistic::~Branch_unrealistic() {
    if(branchingHistory)delete branchingHistory;
}

void Branch_unrealistic::setLeafTimeLimit(double leafTimeLimit) {
    Branch_unrealistic::leafTimeLimit = leafTimeLimit;
}

BranchingHistory * Branch_unrealistic::getHistory() {
    return branchingHistory;
}


char *Branch_unrealistic::getScoreMethodPtr() {
    return &scoreMethod;
}

double *Branch_unrealistic::getAlphaPtr() {
    return &alpha;
}

double *Branch_unrealistic::getEpsPtr() {
    return &epsilon;
}

