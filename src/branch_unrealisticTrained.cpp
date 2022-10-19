//
// Created by simon on 18/10/22.
//
#define SCIP_DEBUG


#include "branch_unrealisticTrained.h"

#define 	BRANCHRULE_NAME   "unrealisticTrained"
#define 	BRANCHRULE_DESC   "Compute the score of a variable using a regression model"
#define 	BRANCHRULE_PRIORITY   200
#define 	BRANCHRULE_MAXDEPTH   -1
#define 	BRANCHRULE_MAXBOUNDDIST   1.0


FeaturesCalculator* Branch_unrealisticTrained::featuresCalculator = nullptr;
RegressionModel* Branch_unrealisticTrained::model = nullptr;


Branch_unrealisticTrained::Branch_unrealisticTrained(SCIP *scip):
        ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                      BRANCHRULE_MAXBOUNDDIST){

}

SCIP_DECL_BRANCHEXECLP(Branch_unrealisticTrained::scip_execlp) {
    assert(featuresCalculator);
    assert(model);

    SCIP_VAR** lpcands;
    SCIP_Real* lpcandsfrac;
    int nlpcands; // number of candidates

    int bestcand = -1;
    double bestScore = 1;

    // get branching candidates
    SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, nullptr, &lpcandsfrac, nullptr, &nlpcands, nullptr) );

    featuresCalculator->computeDynamicProblemFeatures(scip, lpcands, nlpcands);
    // estimate a score for each variable
    for (int i = 0; i < nlpcands; ++i) {
        std::vector<double> features = featuresCalculator->getFeatures(lpcands[i]);
        double score = model->predictScore(features);

        /*SCIPdebugMsg(scip, ( std::to_string(i + 1) + "/" + std::to_string(nlpcands) +
                            " (score: " + std::to_string(score) + ") (var: " + SCIPvarGetName(lpcands[i]) +
                            ")\n").c_str());*/

        if(bestcand == -1 || score > bestScore){
            bestScore = score;
            bestcand = i;
        }
    }

    /*SCIPdebugMsg(scip, ("Var to branch: " + std::to_string(bestcand + 1) + "; " +
                        std::string(SCIPvarGetName(lpcands[bestcand])) + "; Score: " +
                        std::to_string(bestScore) + "\n").c_str());*/

    SCIP_Node* children[2];
    SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], &children[0], NULL, &children[1]) );
    *result = SCIP_BRANCHED;

    featuresCalculator->updateBranchCounter(children, lpcands[bestcand]);
    return SCIP_OKAY;
}

void Branch_unrealisticTrained::setFeaturesCalculator(FeaturesCalculator *newCalculator) {
    if(featuresCalculator)delete featuresCalculator;
    featuresCalculator = newCalculator;
}

void Branch_unrealisticTrained::setModel(RegressionModel *newModel) {
    if(model)delete model;
    model = newModel;
}
