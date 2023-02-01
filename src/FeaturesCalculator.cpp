//
// Created by simon on 11/09/22.
//

#include <scip/scip.h>
#include <iostream>
#include "FeaturesCalculator.h"

FeaturesCalculator::FeaturesCalculator(SCIP *scip, int signB, int signC, int signA) :
        nStaticFeatures(1 + (signC==0) + 2*(1+(signB==0)) + 2*(1+(signC==0)) + 2*(1+3*(signA==0))),
        varnameIndexMap(),
        staticFeaturesMap(),
        dynamicFeaturesMap(),
        numberBrchMap(),
        nvars(SCIPgetNVars(scip)),
        nbrchs(0){
    SCIP_Var **vars = SCIPgetVars(scip);

    sumObjCoefs[0]=0;
    sumObjCoefs[1]=0;
    for(int i=0; i<nvars; ++i){
        varnameIndexMap[SCIPvarGetName(vars[i])] = i;

        double varObjCoef = SCIPvarGetObj(vars[i]);
        if(varObjCoef > 0){
            sumObjCoefs[0] += varObjCoef;
        } else{
            sumObjCoefs[1] += -varObjCoef;
        }
    }

    // compute features for every variable i
    for(int i=0; i<nvars; ++i){
        computStaticFeatures(scip, signB, signC, signA, i);
    }
}

void FeaturesCalculator::computStaticFeatures(SCIP *scip, int signB, int signC, int signA, int i) {
    SCIP_Var **vars = SCIPgetVars(scip);
    int nconss = SCIPgetNOrigConss(scip);
    SCIP_Cons **conss = SCIPgetOrigConss(scip);

    SCIP_Real *consCoefTemp = new SCIP_Real[nvars];
    SCIP_Real *consCoef =new SCIP_Real[nvars];
    SCIP_Var **consVars = new SCIP_Var*[nvars];

    unsigned featureIndex = 0;

    auto* features = new double[nStaticFeatures];
    std::fill(features, features + nStaticFeatures, 0);
    //SCIPallocBlockMemoryArray(scipmain, &features, nfeatures);

    //staticFeaturesMap.insert(std::pair<const char*, const double*>(SCIPvarGetName(vars[i]), features));
    std::string varName = "t_" + std::string(SCIPvarGetName(vars[i]));
    staticFeaturesMap[varName] = features;
    dynamicFeaturesMap[varName] = new double[nDynamicFeatures];
    objectiveIncreaseStaticsMap[varName] = new double[nObjectiveIncreaseStatics];
    std::fill(objectiveIncreaseStaticsMap[varName], objectiveIncreaseStaticsMap[varName] + nObjectiveIncreaseStatics, 0);

    numberBrchMap[varName] = 0;

    double objCoefVari= SCIPvarGetObj(vars[i]);

    computeSet1StaticFeatures(signC, features, featureIndex, objCoefVari);

    // compute second set of features
//init features array to compute min/max
    int sign = 1;
    for(int numMinMaxFeature=0; numMinMaxFeature < nStaticFeatures - featureIndex; ++numMinMaxFeature){
        features[featureIndex + numMinMaxFeature] = sign*DBL_MAX;
        sign *= -1;
    }

    SCIP_Bool success;
    for(int j=0; j<nconss; ++j){
        featureIndex = 1 + (signC==0);
        std::fill(consCoef, consCoef + nvars, 0);
        std::fill(consCoefTemp, consCoefTemp + nvars, 0);
        std::fill(consVars, consVars + nvars, nullptr);
        SCIPgetConsVals(scip, conss[j], consCoefTemp, nvars, &success);
        SCIPgetConsVars(scip, conss[j], consVars, nvars, &success);
        assert(success);

        for(int k=0; k < nvars; ++k){
            if(consVars[k]){
                int index = varnameIndexMap[SCIPvarGetName(consVars[k])];
                consCoef[index] = consCoefTemp[k];
            }
        }

        double bj = SCIPconsGetRhs(scip, conss[j], &success);
        bool canonicalForm = bj < 1e+20; // default upper bound
        if(!canonicalForm){ // constraint not in canonical form, must change it
            continue; // must be a constraint of type x_k >= 0, don't care because must be canonical form
        }
        int localIndexOffset;
        double val;
        computeM1Features(signB, consCoef, i, featureIndex, features, bj);

        int indexOffset = 2*(1+(signB==0));
        computeM2Features(signC, consCoef, i, featureIndex, features, objCoefVari, val,
                          indexOffset);

        //SCIPdebugMessagePrint(scipmain, (varName + ", cons: " + std::to_string(j) + " -> " + std::to_string(features[3] == 1.7976931348623157e+308) + "\n").c_str());

        indexOffset += 2*(1+(signC==0));
        computeM3Features(signA, consCoef, i, featureIndex, features, indexOffset, localIndexOffset, val);
    }

    // Every unmodified field is set to 0
    sign = 1;
    for(int numMinMaxFeature=0; numMinMaxFeature < nStaticFeatures - featureIndex; ++numMinMaxFeature){
        if(features[featureIndex + numMinMaxFeature] == sign*DBL_MAX){
            features[featureIndex + numMinMaxFeature] = 0;
        }
        sign *= -1;
    }

    delete[] consVars;
    delete[] consCoef;
    delete[] consCoefTemp;
}

void FeaturesCalculator::computeSet1StaticFeatures(int signC, double *features, unsigned int &featureIndex,
                                                   double &objCoefVari) const {
    for(int k: {0, 1}) {
        if(signC!=0 && k==1)continue;
        int idx = (signC>=0)?k:1;
        features[featureIndex] = std::abs(objCoefVari) / sumObjCoefs[idx];
        ++featureIndex;
    }
}

void FeaturesCalculator::computeM3Features(int signA, const double *consCoef, int i, unsigned int featureIndex,
                                           double *features, int indexOffset, int &localIndexOffset,
                                           double &val) const {
    // compute M3 -> M3++ M3+- M3-+ M3-- (8 features)
    if(signA==0)
        localIndexOffset = 2 * (consCoef[i] < 0);
    else
        localIndexOffset = 0;

    double ajkSum[2] = {0,0};
    for(int k=0; k < nvars; ++k){
        double ajk = consCoef[k];
        ajkSum[ajk<0] += std::abs(ajk);
    }

    for(auto ajkIndex: {0, 1}){
        // if every element of A has the same sign, we don't need of the part using sum of A's element of the other sign
        if((ajkIndex==0 && signA==-1) || (ajkIndex==1 && signA==1))
            continue;
        //if(signA==0)
        //    ajkIndex = 0;

        if(ajkSum[ajkIndex] != 0) {
            val = std::abs(consCoef[i]) / ajkSum[ajkIndex];
            if (val < features[featureIndex + indexOffset + 4 * ajkIndex + localIndexOffset]) { // keep the min
                features[featureIndex + indexOffset + 4 * ajkIndex + localIndexOffset] = val;
            }
            if (val > features[featureIndex + indexOffset + 4 * ajkIndex + localIndexOffset + 1]) { //keep the max
                features[featureIndex + indexOffset + 4 * ajkIndex + localIndexOffset + 1] = val;
            }
        }
    }
    ////

}

void FeaturesCalculator::computeM2Features(int signC, const double *consCoef, int i, unsigned int featureIndex,
                                           double *features, double objCoefVari, double val, int indexOffset) const {
    // compute M2 (4 features) |c_i|/A_ji
    val = std::abs(objCoefVari) / consCoef[i];
    int localIndexOffset;

    if (signC == 0)
        localIndexOffset = 2 * (objCoefVari < 0);
    else
        localIndexOffset = 0;

    if (val < features[featureIndex + indexOffset + localIndexOffset]) { // keep the min
        features[featureIndex + indexOffset + localIndexOffset] = val;
    }
    if (val > features[featureIndex + indexOffset + localIndexOffset + 1]) { //keep the max
        features[featureIndex + indexOffset + localIndexOffset + 1] = val;
    }
}

void FeaturesCalculator::computeM1Features(int signB, const double *consCoef, int i, unsigned int featureIndex,
                                           double *features, double bj) const {// compute M1 (4 features)
// TODO: M1- utils if c_j >= 0 for all j ???? Maybe less features for this types of problems. Same reflexion for other features set
    int localIndexOffset=0;
    if(signB==0)
    localIndexOffset = 2 * (bj<0);

    double val = consCoef[i] / std::abs(bj);

    if (val < features[featureIndex + localIndexOffset]) { // keep the min
        features[featureIndex + localIndexOffset] = val;
    }
    if (val > features[featureIndex + localIndexOffset + 1]) {//keep the max
        features[featureIndex + localIndexOffset + 1] = val;
    }
}

FeaturesCalculator::~FeaturesCalculator() {
    std::map<std::string, double*>::iterator it;
    for (it = staticFeaturesMap.begin(); it != staticFeaturesMap.end(); it++){
        delete[] it->second;
    }

    for (it = dynamicFeaturesMap.begin(); it != dynamicFeaturesMap.end(); it++){
        delete[] it->second;
    }

    for (it = objectiveIncreaseStaticsMap.begin(); it != objectiveIncreaseStaticsMap.end(); it++){
        delete[] it->second;
    }
}

void FeaturesCalculator::updateBranchCounter(SCIP_NODE **nodes, SCIP_VAR *var) {
    double increase = 0;

    double parentObj = SCIPnodeGetLowerbound(SCIPnodeGetParent(nodes[0]));
    for(auto i:{0,1}){
        double tmp = SCIPnodeGetLowerbound(nodes[i]) - parentObj;
        if(tmp > increase) increase = tmp;
    }

    increase /= parentObj;

    std::string key = std::string(SCIPvarGetName(var));
    int n = numberBrchMap[key] + 1;
    numberBrchMap[key] = n ;
    double* statistics = objectiveIncreaseStaticsMap[key];

    // update min/max
    if(increase < statistics[0] || n == 1){
        statistics[0] = increase;
    }
    if(increase > statistics[1] || n == 1){
        statistics[1] = increase;
    }

    if(n>1){
        // update of the mean
        double oldMean = statistics[2];
        statistics[2] = 1.0/n * (increase + (double)(n-1)*oldMean);

        // update of std
        statistics[3] = ((n-2)*statistics[3] + (n-1)*std::pow(oldMean-statistics[2], 2) + std::pow(increase-statistics[2], 2))/(n-1);
    } else{
        statistics[2] = increase;
        statistics[3] = 0;
    }

    nbrchs++;
    statistics[4] = (double)numberBrchMap[key] / nbrchs;
}


void FeaturesCalculator::computeDynamicProblemFeatures(SCIP *scip, SCIP_Var **vars, int varsSize) {
    SCIP_Var **allVars = SCIPgetVars(scip);

    double* lb = new double[varsSize];
    double* ub = new double[varsSize];
    computeSensitivity(scip, lb, ub, vars, varsSize);

    // count number of fixed variables
    int nFixedVar = 0;
    for(int i=0; i<nvars; ++i){
        if(SCIPvarGetUbLocal(allVars[i]) == SCIPvarGetLbLocal(allVars[i]))
            ++nFixedVar;
    }

    for(int i=0; i<varsSize; ++i){
        std::string key = std::string(SCIPvarGetName(vars[i]));
        auto* features = dynamicFeaturesMap[key];
        features[0] = (double)nFixedVar/nvars;

        // up and down fractionalities of the variable
        double varObj = SCIPvarGetLPSol(vars[i]);
        features[1] = std::ceil(varObj) - varObj;
        features[2] = varObj - std::floor(varObj);

        double varObjCoef = std::abs(SCIPvarGetObj(vars[i]));
        if(varObjCoef!=0){
            features[3] = lb[i]/ varObjCoef;
            features[4] = ub[i]/ varObjCoef;
        } else{
            features[3] = lb[i];
            features[4] = ub[i];
        }
    }

    delete[] lb;
    delete[] ub;
}

const double *FeaturesCalculator::getStaticFeatures(SCIP_VAR *var) {
    return staticFeaturesMap[SCIPvarGetName(var)];
}

const double *FeaturesCalculator::getDynamicProblemFeatures(SCIP_VAR *var) {
    return dynamicFeaturesMap[SCIPvarGetName(var)];
}

const double *FeaturesCalculator::getDynamicOptimizationFeatures(SCIP_VAR *var) {
    return objectiveIncreaseStaticsMap[SCIPvarGetName(var)];
}

const int FeaturesCalculator::getNStaticFeatures() const {
    return nStaticFeatures;
}

const int FeaturesCalculator::getNDynamicFeatures() const {
    return nDynamicFeatures;
}

const int FeaturesCalculator::getNObjectiveIncreaseStatics() const {
    return nObjectiveIncreaseStatics;
}

const std::vector<double> FeaturesCalculator::getFeatures(SCIP_Var *var) {
    std::vector<double> res;

    int arraySizes[3] = {
            getNStaticFeatures(),
            getNDynamicFeatures(),
            getNObjectiveIncreaseStatics()
    };
    const double* features[3] = {
            getStaticFeatures(var),
            getDynamicProblemFeatures(var),
            getDynamicOptimizationFeatures(var)
    };

    for(int i=0; i<3; ++i){
        //line += "start " + std::to_string(i) + ";";
        for(int k=0; k<arraySizes[i]; ++k){
            res.push_back(features[i][k]);
        }
    }

    return res;
}

void FeaturesCalculator::computeSensitivity(SCIP *scip, double *lb, double *ub, SCIP_Var **vars, int varsSize) {
    int nrows=SCIPgetNLPRows(scip);
    int ncols=SCIPgetNLPCols(scip);
    SCIP_Var** allvars = SCIPgetVars(scip);

    std::fill(lb, lb+varsSize, SCIP_REAL_MIN);
    std::fill(ub, ub+varsSize, SCIP_REAL_MAX);

    double* row = new double[ncols];
    double* redcost = new double[ncols];
    for(int i=0; i<ncols; ++i){
        redcost[i] = SCIPgetVarRedcost(scip, allvars[i]);
    }

    int count=nrows-1;
    double eps = SCIPepsilon(scip);
    for(int i=0; i<varsSize; ++i){
        if(redcost[i] != 0)continue;
        double varobj = SCIPvarGetObj(vars[i]);
        int k = count--;
        SCIPgetLPBInvARow(scip, k, NULL, row, NULL, NULL);
        for(int j=0; j<ncols; ++j){
            if(redcost[j] >= -eps && redcost[j] <= eps)continue; // equals 0
            double val = varobj + redcost[j]/row[j];
            if(row[j] < -eps && val > lb[i]){
                lb[i] = val;
            } else if(row[j] > eps && val < ub[i]){
                ub[i] = val;
            }
        }
    }


    delete[] row;
    delete[] redcost;
}

