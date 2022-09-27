//
// Created by simon on 14/08/22.
//
#define SCIP_DEBUG

#include <scip/scipdefplugins.h>
#include <scip/var.h>
#include <sstream>
#include "branch_unrealistic.h"
#include "Utils.h"

#define 	BRANCHRULE_NAME   "unrealistic"
#define 	BRANCHRULE_DESC   "unrealistic branching"
#define 	BRANCHRULE_PRIORITY   200
#define 	BRANCHRULE_MAXDEPTH   -1
#define 	BRANCHRULE_MAXBOUNDDIST   1.0

Branch_unrealistic::Branch_unrealistic(SCIP *scip, int depth, int maxdepth) : ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                                                                   BRANCHRULE_MAXBOUNDDIST), depth(depth), maxdepth(maxdepth){}

Branch_unrealistic::Branch_unrealistic(SCIP *scip, int maxdepth) : ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                                                                                            BRANCHRULE_MAXBOUNDDIST), depth(0), maxdepth(maxdepth){}



SCIP_DECL_BRANCHEXECLP(Branch_unrealistic::scip_execlp){
    SCIP_VAR** lpcands;
    SCIP_Real* lpcandsfrac;
    int nlpcands; // number of candidates
    int bestcand = 0;
    int bestScore = -1;
    SCIP_BoundType bestBranchSide; // which child must be cut from the tree (doesn't lead to the optimal solution)

    // get branching candidates
    SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, nullptr, &lpcandsfrac, nullptr, &nlpcands, nullptr) );
    assert(nlpcands > 0);

    for (int i = 0; i < nlpcands; ++i)lpcandsfrac[i] = SCIPvarGetLPSol(lpcands[i]);


    int* varScores; // store every variable's realNnodes
    SCIPallocBlockMemoryArray(scip, &varScores, nlpcands);

    // an instruction already exists for the first branching (given by the parent)
    if(firstBranch!=nullptr){
        //SCIPdebugMsg(scip, ("nnodes: " + std::to_string(SCIPgetNSols(scip) )+ "\n").c_str());
        SCIP_CALL( SCIPbranchVar(scip, firstBranch, nullptr, NULL, nullptr) );

        firstBranch=nullptr;
        *result = SCIP_BRANCHED;
        if(depth>=maxdepth){
            Utils::configure_scip_instance(scip, false);
        }
        return SCIP_OKAY;
    }

    // computing realNnodes for each variable
    for (int i = 0; i < nlpcands; ++i) {
        int score = INT_MAX;
        SCIP_BoundType branchSide;
        SCIP_Real childPrimalBounds[2];
        assert(lpcands[i] != nullptr);

        SCIP_CALL(computeScore(
                scip,
                score,
                childPrimalBounds,
                bestScore,
                lpcandsfrac[i],
                lpcands[i],
                branchSide));

        varScores[i] = score;

        if (bestScore == -1 || score < bestScore) {
            bestScore = score;
            bestBranchSide = branchSide;
            bestcand = i;
        }

        if (!depth)
            SCIPdebugMsg(scip, (std::string(depth, '\t') + std::to_string(i + 1) + "/" + std::to_string(nlpcands) +
                                " (score: " + std::to_string(score) + ") (var: " + SCIPvarGetName(lpcands[i]) +
                                ")\n").c_str());

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
        dataWriter->addNode(scip, children[bestBranchSide], nlpcands, varScores, lpcands,
                            bestScore, bestcand);
    }

    SCIPfreeBlockMemoryArray(scip, &varScores, nlpcands);
    return SCIP_OKAY;
}



SCIP_RETCODE
Branch_unrealistic::computeScore(SCIP *scip, int &score, SCIP_Real *childPrimalBounds, int bestScore,
                                 SCIP_Real fracValue, SCIP_VAR *varbrch, SCIP_BoundType &branchSide) const {

    score = 0;
    SCIP_Real primalBound=SCIP_REAL_MAX;
    int nodeLimit = (dataWriter && depth==0) || bestScore<=0?-1:bestScore; // if realNnodes data are not used, no need to run more than the best realNnodes

    SCIP *scip_copy;
    SCIP_Bool valid;
    SCIP_CALL(SCIPcreate(&scip_copy));
    SCIP_HashMap *varmap;
    SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip_copy), SCIPgetNVars(scip)) );
    SCIP_CALL( SCIPcopy(scip, scip_copy, varmap, nullptr, "ub_subsolver", FALSE, FALSE, FALSE, FALSE, &valid) );
    assert(valid == TRUE);

    SCIPcopyParamSettings(scip, scip_copy);

    auto *objbranchrule = new Branch_unrealistic(scip_copy, depth + 1, maxdepth);
    SCIP_CALL(SCIPincludeObjBranchrule(scip_copy, objbranchrule, TRUE));
    Utils::configure_scip_instance(scip_copy, true);
    SCIPsetLongintParam(scip_copy, "limits/nodes", nodeLimit);
    SCIP_CALL( SCIPsetIntParam(scip_copy, "display/verblevel",0));

    SCIP_VAR* varbrch_copy = (SCIP_VAR*) SCIPhashmapGetImage(varmap, varbrch);
    assert(varbrch != nullptr);

    objbranchrule->setFirstBranch(varbrch_copy);

    setBestSol(scip, scip_copy, varmap);

    SCIPsolve(scip_copy);

    int localScore = SCIPgetNNodes(scip_copy);

    score = localScore;

    // free memory allocated
    SCIPhashmapFree(&varmap);
    SCIPreleaseVar(scip_copy, &varbrch_copy);
    SCIPfree(&scip_copy);

    return SCIP_OKAY;
}

DatasetWriter* Branch_unrealistic::dataWriter = nullptr;
void Branch_unrealistic::setDataWriter(DatasetWriter *dataWriter) {
    Branch_unrealistic::dataWriter = dataWriter;
}

int *Branch_unrealistic::getMaxDepthPtr() {
    return &maxdepth;
}

void Branch_unrealistic::setFirstBranch(SCIP_Var *firstBranch) {
    Branch_unrealistic::firstBranch = firstBranch;
}

const SCIP_Retcode Branch_unrealistic::setBestSol(SCIP *scip, SCIP *scip_copy, SCIP_HashMap *varmap) const{
    // get the initial vars
    int nvars = SCIPgetNVars(scip);
    SCIP_Var **vars = SCIPgetVars(scip);

    // get values of the best sol
    SCIP_Sol* bestSol = SCIPgetBestSol(scip);
    if(bestSol == NULL) {
        return SCIP_OKAY;
    }
    SCIP_Real* vals;
    SCIPallocBlockMemoryArray(scip, &vals, nvars);
    SCIPgetSolVals(scip, bestSol, nvars, vars, vals);

    // copy the solution for scip_copy
    SCIP_Sol* bestSol_copy;
    SCIP_CALL( SCIPcreateSol(scip_copy, &bestSol_copy, NULL) );

    for(int i=0; i<nvars; ++i){
        SCIP_Var* varCopy =(SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);
        SCIPsetSolVal(scip_copy, bestSol_copy, varCopy, vals[i]);
    }

    unsigned int stored;
    SCIPaddSol(scip_copy, bestSol_copy, &stored);
    assert(stored);

    SCIPfreeBlockMemoryArray(scip, &vals, nvars);

}
