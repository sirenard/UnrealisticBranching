//
// Created by simon on 11/09/22.
//

#include <scip/scip.h>
#include <iostream>
#include "FeaturesCalculator.h"
#include "EventhdlrUpdateFeatures.h"

#define FOR_BRANCH_DIRS(X) for(SCIP_BRANCHDIR dir:{SCIP_BRANCHDIR_DOWNWARDS, SCIP_BRANCHDIR_UPWARDS}){X}

FeaturesCalculator::FeaturesCalculator(SCIP *scip, int signB, int signC, int signA) :
        nvars(SCIPgetNVars(scip)),
        nbrchs(0){

    EventhdlrUpdateFeatures* eventHdlr = dynamic_cast<EventhdlrUpdateFeatures *>(SCIPgetObjEventhdlr(scip,
                                                                                                     SCIPfindEventhdlr(
                                                                                                             scip,
                                                                                                              EVENT_HDLR_UPDATE_FEATURES_NAME)));
    eventHdlr->setFeatureCalculator(this);
}


const std::vector<double> FeaturesCalculator::getFeatures(SCIP_Var *var, SCIP *scip) {
    std::vector<double> res;

    // TODO
    /*int arraySizes[3] = {
            getNStaticFeatures(),
            getNDynamicFeatures(),
            getNObjectiveIncreaseStatics()
    };
    const double* features[3] = {
            getStaticFeatures(var, scip),
            getDynamicProblemFeatures(var),
            getDynamicOptimizationFeatures(var)
    };

    for(int i=0; i<3; ++i){
        //line += "start " + std::to_string(i) + ";";
        for(int k=0; k<arraySizes[i]; ++k){
            res.push_back(features[i][k]);
        }
    }*/

    return res;
}

double FeaturesCalculator::varScore(double s_i, double s_avg) {
    return 1.0 - (1.0/ (1.0 + s_i/std::max(s_avg, 0.1)));
}

double FeaturesCalculator::gNormMax(double x) {
    return std::max(x/(x+1.0), 0.1);
}

double FeaturesCalculator::relDist(double x, double y) {
    if(x<0){
        return 0;
    } else{
        return std::abs(x-y)/(std::max(std::abs(x), std::max(std::abs(y), std::pow(10,-10))));
    }
}

double FeaturesCalculator::relPos(double z, double x, double y) {
    return std::abs(x-z)/std::abs(x-y);
}

double *FeaturesCalculator::getNoTreeFeatures(SCIP *scip, SCIP_Var *var) {
    double* features = new double[nNotreeFeatures];
    int idx = 0;

    double lpsol = SCIPvarGetLPSol(var);
    features[idx++] = lpsol;
    features[idx++] = SCIPvarGetAvgSol(var);

    FOR_BRANCH_DIRS(features[idx++] = 1.0 - (double)(1+SCIPvarGetAvgBranchdepthCurrentRun(var, dir)) / (double)(1+SCIPgetMaxDepth(scip));)

    features[idx++] = varScore(SCIPgetVarConflictScore(scip, var), SCIPgetAvgConflictScore(scip));
    features[idx++] = varScore(SCIPgetVarConflictlengthScore(scip, var), SCIPgetAvgConflictlengthScore(scip));
    features[idx++] = varScore(SCIPgetVarAvgInferenceScore(scip, var), SCIPgetAvgInferenceScore(scip));
    features[idx++] = varScore(SCIPgetVarAvgCutoffScore(scip, var), SCIPgetAvgCutoffScore(scip));
    features[idx++] = varScore(SCIPgetVarPseudocostScore(scip, var, lpsol), SCIPgetAvgPseudocostScore(scip));

    FOR_BRANCH_DIRS(features[idx++] = (1.0+SCIPgetVarPseudocostCountCurrentRun(scip, var, dir))/ (1.0+SCIPgetPseudocostCount(scip, dir, 1));)
    FOR_BRANCH_DIRS(features[idx++] = (1.0+SCIPgetVarPseudocostCountCurrentRun(scip, var, dir))/ (1.0+SCIPvarGetNBranchingsCurrentRun(var, dir));)
    FOR_BRANCH_DIRS(features[idx++] = (1.0+SCIPgetVarPseudocostCountCurrentRun(scip, var, dir))/ (1.0+nbrchs);)

    for(int i:{0,1})features[idx++] = SCIPvarGetNImpls(var, i);

    for(int i:{0,1})features[idx++] = (1.0+SCIPvarGetNCliques(var, i)) / (1.0+SCIPgetNCliques(scip));
    FOR_BRANCH_DIRS(features[idx++] = gNormMax(SCIPgetVarAvgCutoffsCurrentRun(scip, var, dir));)
    FOR_BRANCH_DIRS(features[idx++] = gNormMax(SCIPgetVarAvgConflictlengthCurrentRun(scip, var, dir));)
    FOR_BRANCH_DIRS(features[idx++] = gNormMax(SCIPgetVarAvgInferencesCurrentRun(scip, var, dir));)

    return features;
}

double *FeaturesCalculator::getTreeFeatures(SCIP *scip, int nlpcands) {
    double* features = new double[nTreeFeatures];
    int idx = 0;

    SCIP_Node* node = SCIPgetCurrentNode(scip);
    features[idx++] = (1.0+SCIPnodeGetDepth(node)) / (1.0+SCIPgetMaxDepth(scip));
    features[idx++] = (1.0+SCIPgetPlungeDepth(scip)) / (1.0+SCIPgetMaxDepth(scip));
    features[idx++] = relDist(SCIPgetLowerbound(scip), SCIPgetLPObjval(scip));
    features[idx++] = relDist(SCIPgetLowerboundRoot(scip), SCIPgetLPObjval(scip));
    features[idx++] = relDist(SCIPgetUpperbound(scip), SCIPgetLPObjval(scip));
    features[idx++] = relPos(SCIPgetLPObjval(scip), SCIPgetUpperbound(scip), SCIPgetLowerbound(scip));
    features[idx++] = (1.0+nlpcands)/(1.0+nvars-nlpcands);

    int nleaves = SCIPgetNLeaves(scip);
    int ninternalnodes = SCIPgetNNodes(scip) - nleaves;
    int ncreatednode = SCIPgetNNodes(scip); // NOT SURE
    int nactivatednodes = 0; // tode
    int ndeactivatednodes = SCIPgetNNodes(scip) - nactivatednodes;
    features[idx++] = (1.0+SCIPgetNObjlimLeaves(scip))/(1.0+nleaves);
    features[idx++] = (1.0+SCIPgetNInfeasibleLeaves(scip))/(1.0+nleaves);
    features[idx++] = (1.0+SCIPgetNFeasibleLeaves(scip))/(1.0+nleaves);
    features[idx++] = (double) (SCIPgetNInfeasibleLeaves(scip)+1)/(double)(SCIPgetNObjlimLeaves(scip)+1);
    features[idx++] = (double) SCIPgetNNodesLeft(scip)/(double)SCIPgetNNodes(scip);
    features[idx++] = (double) ninternalnodes/(double)SCIPgetNNodes(scip);
    features[idx++] = (double) SCIPgetNNodes(scip)/(double)ncreatednode;

    features[idx++] = (double) nactivatednodes/(double) SCIPgetNNodes(scip);
    features[idx++] = (double) ndeactivatednodes/(double) SCIPgetNNodes(scip);
    features[idx++] = (1.0+SCIPgetPlungeDepth(scip))/(1.0+SCIPgetMaxDepth(scip));
    features[idx++] = (1.0+SCIPgetNBacktracks(scip))/(1.0+SCIPgetNNodes(scip));

    features[idx++] = std::log((double) SCIPgetNLPIterations(scip)/(double)SCIPgetNNodes(scip));
    features[idx++] = std::log((double) SCIPgetNLPs(scip)/(double)SCIPgetNNodes(scip));
    features[idx++] = (double)SCIPgetNNodes(scip)/(double) SCIPgetNLPs(scip);
    features[idx++] = (double)SCIPgetNNodeLPs(scip)/(double) SCIPgetNLPs(scip);

    //TODO: Gap(4)
    int primaldualintegral = 1; //TODO
    double lastsolgap = 1; // TODO
    double firstsolgap = 1;
    features[idx++] = std::log(primaldualintegral);
    features[idx++] = SCIPgetGap(scip) / lastsolgap;
    features[idx++] = SCIPgetGap(scip) / firstsolgap;


    features[idx++] = relDist(SCIPgetLowerboundRoot(scip), SCIPgetLowerbound(scip));
    features[idx++] = relDist(SCIPgetLowerboundRoot(scip), SCIPgetAvgLowerbound(scip));
    features[idx++] = relDist(SCIPgetUpperbound(scip), SCIPgetLowerbound(scip));
    features[idx++] = SCIPisPrimalboundSol(scip);
    int nnodesbeforefirst = 1;
    features[idx++] = (double)nnodesbeforefirst / SCIPgetNNodes(scip);

    features[idx++] = gNormMax(SCIPgetAvgConflictScore(scip));
    features[idx++] = gNormMax(SCIPgetAvgConflictlengthScore(scip));
    features[idx++] = gNormMax(SCIPgetAvgInferenceScore(scip));
    features[idx++] = gNormMax(SCIPgetAvgCutoffScore(scip));
    features[idx++] = gNormMax(SCIPgetAvgPseudocostScore(scip));
    features[idx++] = gNormMax(SCIPgetAvgCutoffScore(scip));
    FOR_BRANCH_DIRS(features[idx++] = gNormMax(SCIPgetAvgCutoffs(scip, dir));)
    FOR_BRANCH_DIRS(features[idx++] = gNormMax(SCIPgetAvgInferences(scip, dir));)
    FOR_BRANCH_DIRS(features[idx++] = gNormMax(SCIPgetPseudocostVariance(scip, dir, 1));)
    FOR_BRANCH_DIRS(features[idx++] = gNormMax(SCIPgetNConflictConssApplied(scip));)

    //TODO: Open nodes bounds (12), Open nodes depths (4)

    features[idx++] = 0; // TODO
    features[idx++] = relDist(SCIPgetLowerbound(scip), open_lbs_max);
    features[idx++] = relDist(open_lbs_min, open_lbs_max);
    features[idx++] = relDist(open_lbs_min, SCIPgetUpperbound(scip));
    features[idx++] = relPos(open_lbs_mean, SCIPgetUpperbound(scip), SCIPgetLowerbound(scip));
    features[idx++] = relPos(open_lbs_min, SCIPgetUpperbound(scip), SCIPgetLowerbound(scip));
    features[idx++] = relPos(open_lbs_max, SCIPgetUpperbound(scip), SCIPgetLowerbound(scip));
    features[idx++] = 0; // TODO
    features[idx++] = open_lbs_std / open_lbs_mean;
    features[idx++] = 0; // TODO

    features[idx++] = (double)open_ds_mean/ (1.0+SCIPgetMaxDepth(scip));
    features[idx++] = 0; // TODO
    features[idx++] = open_ds_std / open_ds_mean;
    features[idx++] = 0; // TODO

    return features;
}

void FeaturesCalculator::updateBranching(SCIP *scip) {
    std::cout << "Inform" << std::endl;
    double lb = SCIPgetLocalLowerbound(scip);
    int depth = SCIPgetDepth(scip)+1;
    nbrchs++;
    if(nbrchs == 1){
        open_lbs_max = lb;
        open_lbs_min = lb;
        open_lbs_mean = lb;
        open_lbs_std = 0;
        open_lbs_squared_sum = lb*lb;

        open_ds_max = depth;
        open_ds_min = depth;
        open_ds_mean = depth;
        open_ds_std = 0;
        open_ds_squared_sum = depth*depth;
    } else{
        open_lbs_squared_sum += lb*lb;
        open_ds_squared_sum += depth*depth;

        if(lb > open_lbs_max)open_lbs_max = lb;
        if(lb < open_lbs_min)open_lbs_min = lb;
        double open_lbs_sum = open_lbs_mean*(nbrchs-1);
        open_lbs_mean = (open_lbs_sum + lb)/nbrchs ;
        open_lbs_std = std::sqrt((nbrchs * open_lbs_squared_sum - open_lbs_sum*open_lbs_sum)/(nbrchs * (nbrchs+1)));


        if(lb > open_ds_max)open_ds_max = depth;
        if(lb < open_ds_min)open_ds_min = depth;
        double open_ds_sum = open_ds_mean*(nbrchs-1);
        open_ds_mean = (open_ds_sum + lb)/nbrchs ;
        open_ds_std = std::sqrt((nbrchs * open_ds_squared_sum - open_ds_sum*open_ds_sum)/(nbrchs * (nbrchs+1)));
    }
}

int FeaturesCalculator::getNoTreeFeaturesSize() {
    return nNotreeFeatures;
}

int FeaturesCalculator::getTreeFeaturesSize() {
    return nTreeFeatures;
}

