//
// Created by simon on 14/08/22.
//
#define SCIP_DEBUG

#include <scip/scipdefplugins.h>
#include <scip/var.h>
#include <sstream>
#include <thread>
#include "branch_unrealistic.h"
#include "Utils.h"
#include "mpi/Worker.h"

#define 	BRANCHRULE_NAME   "unrealistic"
#define 	BRANCHRULE_DESC   "unrealistic branching"
#define 	BRANCHRULE_PRIORITY   200
#define 	BRANCHRULE_MAXDEPTH   -1
#define 	BRANCHRULE_MAXBOUNDDIST   1.0

Branch_unrealistic::Branch_unrealistic(SCIP *scip, int depth, int maxdepth, double leafTimeLimit) : ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                                                                   BRANCHRULE_MAXBOUNDDIST), depth(depth), maxdepth(maxdepth), firstBranch(nullptr), leafTimeLimit(leafTimeLimit){}

Branch_unrealistic::Branch_unrealistic(SCIP *scip, int maxdepth, double leafTimeLimit) : ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                                                                                            BRANCHRULE_MAXBOUNDDIST), depth(0), maxdepth(maxdepth), firstBranch(nullptr), leafTimeLimit(leafTimeLimit){}



SCIP_DECL_BRANCHEXECLP(Branch_unrealistic::scip_execlp){
    SCIP_VAR** lpcands;
    SCIP_Real* lpcandsfrac;
    int nlpcands; // number of candidates
    std::vector<int> bestcands;
    int bestScore = -1;
    SCIP_BoundType bestBranchSide; // which child must be cut from the tree (doesn't lead to the optimal solution)

    // get branching candidates
    SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, nullptr, &lpcandsfrac, nullptr, &nlpcands, nullptr) );
    assert(nlpcands > 0);

    for (int i = 0; i < nlpcands; ++i)lpcandsfrac[i] = SCIPvarGetLPSol(lpcands[i]);


    int* varScores = nullptr; // store every variable's realNnodes
    if(dataWriter != nullptr && depth == 0)SCIP_CALL(SCIPallocBufferArray(scip, &varScores, nlpcands));

    // an instruction already exists for the first branching (given by the parent)
    if(firstBranch!=nullptr){
        SCIP_CALL( SCIPbranchVarHole(scip, firstBranch, left, right, nullptr, nullptr) );
        firstBranch=nullptr;

        *result = SCIP_BRANCHED;
        SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
        if(depth>=maxdepth){
            SCIP_CALL( Utils::configure_scip_instance(scip, false) );
            SCIP_CALL( SCIPsetRealParam(scip, "limits/time", leafTimeLimit));
        } else{
            SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
        }
        return SCIP_OKAY;
    }


    Worker* worker = Worker::getInstance();
    worker->computeScores(scip, lpcands, nlpcands, bestcands, bestScore, depth + 1, maxdepth, leafTimeLimit,
                          dataWriter != nullptr && depth == 0, varScores);

    int bestcand = bestcands[rand() % bestcands.size()];
    //int bestcand = bestcands.at(0); // TODO: chose at random
    for(auto k: bestcands){
        if(k<bestcand) bestcand = k;
    }
    if (!depth) {
        SCIPdebugMsg(scip, ("Var to branch: " + std::to_string(bestcand + 1) + "; " +
                            std::string(SCIPvarGetName(lpcands[bestcand])) + "; Score: " +
                            std::to_string(bestScore) + "\n").c_str());
    }

    SCIP_Node* children[2];
    SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], &children[0], NULL, &children[1]) );
    *result = SCIP_BRANCHED;


    if(dataWriter && depth==0) {
        dataWriter->addNode(scip, children, nlpcands, varScores, lpcands, bestcand);
    }

    if(dataWriter != nullptr && depth == 0)SCIPfreeBufferArray(scip, &varScores);
    return SCIP_OKAY;
}


DatasetWriter* Branch_unrealistic::dataWriter = nullptr;
void Branch_unrealistic::setDataWriter(DatasetWriter *dataWriter) {
    Branch_unrealistic::dataWriter = dataWriter;
}

int *Branch_unrealistic::getMaxDepthPtr() {
    return &maxdepth;
}

void Branch_unrealistic::setFirstBranch(SCIP_Var *firstBranch, double left, double right) {
    Branch_unrealistic::firstBranch = firstBranch;
    this->left = left;
    this->right = right;
}

double *Branch_unrealistic::getLeafTimeLimitPtr() {
    return &leafTimeLimit;
}

SCIP_DECL_BRANCHINIT(Branch_unrealistic::scip_init){
    SCIP_CALL( scip::ObjBranchrule::scip_init(scip, branchrule) );
    Worker *worker = Worker::getInstance();
    if(worker->isMaster() && depth==0){
        worker->setScipInstance(scip);
    }
    return SCIP_OKAY;
}

SCIP_DECL_BRANCHEXIT(Branch_unrealistic::scip_exit){
    SCIP_CALL( scip::ObjBranchrule::scip_init(scip, branchrule) );
    Worker *worker = Worker::getInstance();
    if(worker->isMaster() && depth==0){
        worker->broadcastEnd();

        SCIP_Var** vars = SCIPgetVars(scip);
        int n = SCIPgetNVars(scip);
        for(int i=0; i<n; ++i){
            SCIPvarSetRemovable(vars[i], FALSE);
        }
    }
    return SCIP_OKAY;
}