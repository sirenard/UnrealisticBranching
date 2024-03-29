//
// Created by simon on 11/09/22.
//

#include <scip/scip.h>
#include <iostream>
#include "FeaturesCalculator.h"
#include <scip/struct_lp.h>

FeaturesCalculator::FeaturesCalculator(SCIP *scip, int signB, int signC, int signA) :
        nStaticFeatures(1 + (signC==0) + 2*(1+(signB==0)) + 2*(1+(signC==0)) + 2*(1+3*(signA==0))),
        signA(signA),
        signB(signB),
        signC(signC),
        varnameIndexMap(),
        dynamicFeatures(),
        numberBrchs(),
        nvars(SCIPgetNVars(scip)),
        nbrchs(0){
    SCIP_Var **vars = SCIPgetVars(scip);

    staticFeatures = new double*[nvars];
    dynamicFeatures = new double*[nvars];
    objectiveIncreaseStaticsMap = new double*[nvars];
    numberBrchs = new int[nvars];
    std::fill(staticFeatures, staticFeatures + nvars, nullptr);
    std::fill(dynamicFeatures, dynamicFeatures + nvars, nullptr);
    std::fill(numberBrchs, numberBrchs + nvars, 0);

    sumObjCoefs[0]=0;
    sumObjCoefs[1]=0;
    for(int i=0; i<nvars; ++i){
        std::string key = SCIPvarGetName(vars[i]);
        if(key.find("t_")==0){
            key.erase(0, 2);
        }
        varnameIndexMap[key] = i;
        dynamicFeatures[i] = new double[nDynamicFeatures];
        objectiveIncreaseStaticsMap[i] = new double[nObjectiveIncreaseStatics];
        std::fill(objectiveIncreaseStaticsMap[i], objectiveIncreaseStaticsMap[i] + nObjectiveIncreaseStatics, 0);

        double varObjCoef = SCIPvarGetObj(vars[i]);
        if(varObjCoef > 0){
            sumObjCoefs[0] += varObjCoef;
        } else{
            sumObjCoefs[1] += -varObjCoef;
        }
    }
}

void FeaturesCalculator::computStaticFeatures(SCIP *scip, SCIP_Var *var) {
    int nconss = SCIPgetNOrigConss(scip);
    SCIP_Cons **conss = SCIPgetOrigConss(scip);

    SCIP_Real *consCoefTemp = new SCIP_Real[nvars];
    SCIP_Real *consCoef =new SCIP_Real[nvars];
    SCIP_Var **consVars = new SCIP_Var*[nvars];

    unsigned featureIndex = 0;

    auto* features = new double[nStaticFeatures];
    std::fill(features, features + nStaticFeatures, 0);

    int index = getVarKey(var);
    staticFeatures[index] = features;

    numberBrchs[index] = 0;

    double objCoefVari= SCIPvarGetObj(var);

    computeSet1StaticFeatures(features, featureIndex, objCoefVari);

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

        int nVarsInCons;
        SCIPgetConsNVars(scip, conss[j], &nVarsInCons, &success);
        for(int k=0; k < nVarsInCons; ++k){
            int index = getVarKey(consVars[k]);
            consCoef[index] = consCoefTemp[k];
        }

        double bj = SCIPconsGetRhs(scip, conss[j], &success);
        bool canonicalForm = bj < 1e+20; // default upper bound
        if(!canonicalForm){ // constraint not in canonical form, must change it
            continue; // must be a constraint of type x_k >= 0, don't care because must be canonical form
        }
        int localIndexOffset;
        double val;
        computeM1Features(consCoef, index, featureIndex, features, bj);

        int indexOffset = 2*(1+(signB==0));
        computeM2Features(consCoef, index, featureIndex, features, objCoefVari, val,
                          indexOffset);

        //SCIPdebugMessagePrint(scipmain, (index + ", cons: " + std::to_string(j) + " -> " + std::to_string(features[3] == 1.7976931348623157e+308) + "\n").c_str());

        indexOffset += 2*(1+(signC==0));
        computeM3Features(consCoef, index, featureIndex, features, indexOffset, localIndexOffset, val);
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

void FeaturesCalculator::computeSet1StaticFeatures(double *features, unsigned int &featureIndex,
                                                   double &objCoefVari) const {
    for(int k: {0, 1}) {
        if(signC!=0 && k==1)continue;
        int idx = (signC>=0)?k:1;
        features[featureIndex] = std::abs(objCoefVari) / sumObjCoefs[idx];
        ++featureIndex;
    }
}

void FeaturesCalculator::computeM3Features(const double *consCoef, int i, unsigned int featureIndex, double *features,
                                           int indexOffset, int &localIndexOffset, double &val) const {
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

void FeaturesCalculator::computeM2Features(const double *consCoef, int i, unsigned int featureIndex,
                                           double *features, double objCoefVari, double val, int indexOffset) const {
    // compute M2 (4 features) |c_i|/A_ji
    val = std::abs(objCoefVari) / consCoef[i];
    if(consCoef[i]==0){
        val=0;
    }
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

void FeaturesCalculator::computeM1Features(const double *consCoef, int i, unsigned int featureIndex, double *features,
                                           double bj) const {// compute M1 (4 features)
// TODO: M1- utils if c_j >= 0 for all j ???? Maybe less features for this types of problems. Same reflexion for other features set
    int localIndexOffset=0;
    if(signB==0)
        localIndexOffset = 2 * (bj<0);

    double val = consCoef[i] / std::abs(bj);
    if(bj == 0){
        val=0;
    }

    if (val < features[featureIndex + localIndexOffset]) { // keep the min
        features[featureIndex + localIndexOffset] = val;
    }
    if (val > features[featureIndex + localIndexOffset + 1]) {//keep the max
        features[featureIndex + localIndexOffset + 1] = val;
    }
}

FeaturesCalculator::~FeaturesCalculator() {
    //std::map<std::string, double*>::iterator it;
    for (int i=0; i<nvars; ++i){
        delete[] staticFeatures[i];
        delete[] dynamicFeatures[i];
        delete[] objectiveIncreaseStaticsMap[i];
    }
    delete[] staticFeatures;
    delete[] dynamicFeatures;
    delete[] objectiveIncreaseStaticsMap;

    /*for (it = dynamicFeatures.begin(); it != dynamicFeatures.end(); it++){
        delete[] it->second;
    }
    for (it = objectiveIncreaseStaticsMap.begin(); it != objectiveIncreaseStaticsMap.end(); it++){
        delete[] it->second;
    }*/
}

void FeaturesCalculator::updateBranchCounter(SCIP_NODE **nodes, SCIP_VAR *var) {
    double increase = 1;

    double parentObj = SCIPnodeGetLowerbound(SCIPnodeGetParent(nodes[0]));
    for(auto i:{0,1}){
        double lb = SCIPnodeGetLowerbound(nodes[i]);
        if (lb==1e20){
            lb=1;
        }
        double tmp = lb - parentObj;
        increase *= std::max(tmp, 10e-3);
    }

    increase /= parentObj;

    int index = getVarKey(var);
    int n = numberBrchs[index] + 1;
    numberBrchs[index] = n ;
    double* statistics = objectiveIncreaseStaticsMap[index];

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
    statistics[4] = (double)numberBrchs[index] / nbrchs;
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
        int index = getVarKey(vars[i]);
        auto* features = dynamicFeatures[index];
        features[0] = (double)nFixedVar/nvars;

        // up and down fractionalities of the variable
        double varObj = SCIPvarGetLPSol(vars[i]);
        features[1] = std::ceil(varObj) - varObj;
        features[2] = varObj - std::floor(varObj);

        double varObjCoef = std::abs(SCIPvarGetObj(vars[i]));
        if(varObjCoef==0)varObjCoef =1;
        /*std::cout << SCIPvarGetName(vars[i]) << "="<< SCIPvarGetLPSol(vars[i]) << std::endl;
        std::cout << "LB: " << lb[i] << std::endl;
        std::cout << "UB: " << ub[i] << std::endl;*/

        if(lb[i] != SCIP_REAL_MAX && lb[i] != SCIP_REAL_MIN)
            features[3] = lb[i]/ varObjCoef;
        else
            features[3] = 0;
        if(ub[i] != SCIP_REAL_MAX && ub[i] != SCIP_REAL_MIN)
            features[4] = ub[i]/ varObjCoef;
        else
            features[4] = 0;

    }

    delete[] lb;
    delete[] ub;
}

const double *FeaturesCalculator::getStaticFeatures(SCIP_VAR *var, SCIP *scip) {
    int index = getVarKey(var);
    if(staticFeatures[index] == nullptr){
        computStaticFeatures(scip, var);
    }

    return staticFeatures[index];
}

const double *FeaturesCalculator::getDynamicProblemFeatures(SCIP_VAR *var) {
    int index = getVarKey(var);
    return dynamicFeatures[index];
}

const double *FeaturesCalculator::getDynamicOptimizationFeatures(SCIP_VAR *var) {
    int index = getVarKey(var);
    return objectiveIncreaseStaticsMap[index];
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

void FeaturesCalculator::getFeatures(SCIP_Var *var, SCIP *scip, double *features) {
    int idx=0;
    memcpy(features, getStaticFeatures(var, scip), getNStaticFeatures()*sizeof(double));
    idx += getNStaticFeatures();
    memcpy(&features[idx], getDynamicProblemFeatures(var), getNDynamicFeatures()*sizeof(double));
    idx += getNDynamicFeatures();
    memcpy(&features[idx], getDynamicOptimizationFeatures(var), getNObjectiveIncreaseStatics()*sizeof(double));
}

void FeaturesCalculator::computeSensitivity(SCIP *scip, double *lb, double *ub, SCIP_Var **vars, int varsSize) {
    int nrows=SCIPgetNLPRows(scip);
    int ncols=SCIPgetNLPCols(scip);
    SCIP_Var** allvars = SCIPgetVars(scip);

    std::fill(lb, lb+varsSize, SCIP_REAL_MIN);
    std::fill(ub, ub+varsSize, SCIP_REAL_MAX);

    double* row = new double[ncols];
    double* col = new double[nrows];
    double* redcost = new double[ncols];
    for(int i=0; i<ncols; ++i){
        redcost[i] = SCIPgetVarRedcost(scip, allvars[i]);
    }

    int count=nrows-1;
    double eps = SCIPepsilon(scip);
    for(int i=0; i<varsSize; ++i){
        double rc = SCIPgetVarRedcost(scip, vars[i]);
        if(rc != 0)continue;
        double varobj = SCIPvarGetObj(vars[i]);
//        int k = count--;
        SCIPgetLPBInvACol(scip, SCIPvarGetCol(vars[i])->index, col, NULL, NULL);
        int k;
        for(k=0; k<nrows; ++k){
            if(col[k]>=1.0-eps)break;
        }

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
    delete[] col;
    delete[] redcost;
}

int FeaturesCalculator::getVarKey(SCIP_Var *var) {
    //int colNum = SCIPcolGetLPPos(SCIPvarGetCol(var));
    std::string key = SCIPvarGetName(var);
    if(key.find("t_")==0){
        key.erase(0, 2);
    }
    return varnameIndexMap[key];
}

const int FeaturesCalculator::getNFeatures() const {
    return getNDynamicFeatures() + getNObjectiveIncreaseStatics() + getNStaticFeatures();
}
