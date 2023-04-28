//
// Created by simon on 18/10/22.
//
#define SCIP_DEBUG


#include "branch_unrealisticTrained.h"

#include <random>

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

    int bestcand = 0;
    double bestScore=0;


    // get branching candidates
    SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, nullptr, &lpcandsfrac, nullptr, &nlpcands, nullptr) );

    std::shuffle(lpcands, lpcands+nlpcands, std::mt19937(std::random_device()()));


    //featuresCalculator->updateBranching(scip);

    // estimate a score for each variable

    int modelInputSize = model->getInputSize();
    featuresCalculator->computeDynamicProblemFeatures(scip, lpcands, nlpcands);
    int featuresSize = featuresCalculator->getNFeatures();
    double *features = new double[modelInputSize];
    // estimate a score for each variable
    for (int i = 0; i < nlpcands; ++i) {
        if(modelInputSize == 2*featuresSize) { // use 1v1 approach
            if(i==0){
                featuresCalculator->getFeatures(lpcands[0], scip, &features[featuresSize]);
                continue;
            } else{
                featuresCalculator->getFeatures(lpcands[i], scip, features);
            }

            double score = model->predictScore(features);
            if(score > alpha){
                bestcand = i;
                memcpy( &features[featuresSize], features, featuresSize*sizeof(double));
            }
        } else{
            featuresCalculator->getFeatures(lpcands[i], scip, features);
            double score = model->predictScore(features);
            if(score > bestScore){
                bestcand = i;
                bestScore = score;
            }
        }





    }

    SCIP_Node *children[2];
    SCIP_CALL(SCIPbranchVar(scip, lpcands[bestcand], &children[0], NULL, &children[1]));
    *result = SCIP_BRANCHED;

    delete[] features;
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

double *Branch_unrealisticTrained::getAlphaPtr() {
    return &alpha;
}
