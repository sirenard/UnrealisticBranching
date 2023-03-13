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

    featuresCalculator->updateBranching(scip);

    // estimate a score for each variable

    double* nodeFeatures = featuresCalculator->getTreeFeatures(scip, nlpcands);
    int nodeFeaturesSize = featuresCalculator->getTreeFeaturesSize();
    int varFeaturesSize = featuresCalculator->getNoTreeFeaturesSize();

    std::vector<double> features(nodeFeaturesSize + varFeaturesSize);
    for(int k=0;k<nodeFeaturesSize; ++k){
        features.at(k) = nodeFeatures[k];
    }

    for (int i = 0; i < nlpcands; ++i) {
        double* varFeatures = featuresCalculator->getNoTreeFeatures(scip, lpcands[i]);

        for(int k=0;k<varFeaturesSize; ++k){
            features.at(nodeFeaturesSize+k) = varFeatures[k];
        }

        double score = model->predictScore(features);

        if(bestcand == -1 || score > bestScore){
            bestScore = score;
            bestcand = i;
        }

        delete[] varFeatures;
    }


    SCIP_Node* children[2];
    SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], &children[0], NULL, &children[1]) );
    *result = SCIP_BRANCHED;

    delete[] nodeFeatures;
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
