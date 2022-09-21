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

    int* varScores; // store every variable's realNnodes
    SCIPallocBlockMemoryArray(scip, &varScores, nlpcands);

    // computing realNnodes for each variable
    for(int i=0; i<nlpcands; ++i){
        int score = INT_MAX;
        SCIP_Real childPrimalBounds[2];
        assert(lpcands[i] != nullptr);

        SCIP_CALL(computeScore(
                scip,
                score,
                childPrimalBounds,
                bestScore,
                lpcandsfrac[i],
                lpcands[i]));

        varScores[i] = score;
        //varScores[i] = 0;
        if(varScores[i] == 0){
            continue;
        }

        if(bestScore == -1 || score < bestScore){
            bestScore = score;
            bestcand = i;
        }

        if(!depth)
            SCIPdebugMsg(scip, (std::string(depth, '\t') + std::to_string(i+1) +"/" + std::to_string(nlpcands) + " (score: " + std::to_string(score) +")\n").c_str());
    }


    if(!depth) {
        SCIPdebugMsg(scip, ("Var to branch: " + std::to_string(bestcand + 1) + "; " +
                            std::string(SCIPvarGetName(lpcands[bestcand])) + "; Score: " +
                            std::to_string(bestScore) + "\n").c_str());
    }

    SCIP_Node* children[2];
    SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], &children[0], NULL, &children[1]) );
    *result = SCIP_BRANCHED;


    if(dataWriter && depth==0){
        //dataWriter->addNode(scip, children[bestBranchSide != SCIP_BOUNDTYPE_UPPER], nlpcands, varScores, lpcands,
        //                    bestScore, bestcand);
    }

    SCIPfreeBlockMemoryArray(scip, &varScores, nlpcands);
    return SCIP_OKAY;
}



SCIP_RETCODE
Branch_unrealistic::computeScore(SCIP *scip, int &score, SCIP_Real *childPrimalBounds, int bestScore,
                                 SCIP_Real fracValue, SCIP_VAR *varbrch) const {

    score = 0;
    //nodeLimit = -1;

    for(auto bound : {SCIP_BOUNDTYPE_UPPER, SCIP_BOUNDTYPE_LOWER}) {
        int nodeLimit = (dataWriter && depth==0) || bestScore<=0?-1:bestScore+1-score; // if realNnodes data are not used, no need to run more than the best realNnodes

        SCIP *scip_copy;
        SCIP_Bool valid;
        SCIP_CALL(SCIPcreate(&scip_copy));
        SCIP_HashMap *varmap;
        SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip_copy), SCIPgetNVars(scip)) );
        SCIP_CALL( SCIPcopy(scip, scip_copy, varmap, nullptr, "ub_subsolver", FALSE, FALSE, FALSE, FALSE, &valid) );

        SCIPcopyParamSettings(scip, scip_copy);

        auto *objbranchrule = new Branch_unrealistic(scip_copy, depth + 1, maxdepth);
        SCIP_CALL(SCIPincludeObjBranchrule(scip_copy, objbranchrule, TRUE));
        Utils::configure_scip_instance(scip_copy, depth+1<maxdepth);
        SCIPsetLongintParam(scip_copy, "limits/nodes", nodeLimit);

        SCIP_CALL( SCIPsetIntParam(scip_copy, "display/verblevel",0));


        assert(valid == TRUE);

        SCIP_VAR* varbrch_copy = (SCIP_VAR*) SCIPhashmapGetImage(varmap, varbrch);
        assert(varbrch != nullptr);

        // add the new constraint
        SCIP_CONS *cons;
        SCIP_Real coef = 1;
        switch(bound){
            case SCIP_BOUNDTYPE_LOWER:
                SCIPcreateConsBasicLinear(scip_copy, &cons, "SIMON_LOWER", 1, &varbrch_copy,
                                          &coef, std::ceil(fracValue),SCIP_REAL_MAX);
                break;
            case SCIP_BOUNDTYPE_UPPER:
                SCIPcreateConsBasicLinear(scip_copy, &cons, "SIMON_UPPER", 1, &varbrch_copy,
                                          &coef, SCIP_REAL_MIN, std::floor(fracValue));
                break;
        }
        assert(cons!=nullptr);
        SCIPaddCons(scip_copy, cons);

        SCIPsolve(scip_copy);
        SCIP_Bool feasible;
        SCIP_CALL( SCIPvalidateSolve(scip, SCIP_UNKNOWN, SCIP_UNKNOWN, SCIPfeastol(scip), TRUE, &feasible, NULL, NULL) );

        SCIP_Real localPrimalBound = SCIPgetPrimalbound(scip_copy);
        childPrimalBounds[bound] = localPrimalBound;

        int localScore = SCIPgetNNodes(scip_copy);

        auto status = SCIPgetStatus(scip_copy);

        if(status != SCIP_STATUS_INFEASIBLE){
            score += localScore;
        } else{
            score += 1;
        }

        if(status == SCIP_STATUS_TIMELIMIT){
            score = INT_MAX;
            break;
        }

        if(score < 25){
            int ii=2;
        }


        SCIP_CALL( SCIPreleaseCons(scip_copy, &cons) );

        SCIPhashmapFree(&varmap);
        SCIPreleaseVar(scip_copy, &varbrch_copy);
        SCIPfree(&scip_copy);

    }

    return SCIP_OKAY;
}

DatasetWriter* Branch_unrealistic::dataWriter = nullptr;
void Branch_unrealistic::setDataWriter(DatasetWriter *dataWriter) {
    Branch_unrealistic::dataWriter = dataWriter;
}
