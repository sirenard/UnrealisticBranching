//
// Created by simon on 11/09/22.
//

#include <scip/scip.h>
#include <iostream>
#include "FeaturesCalculator.h"

FeaturesCalculator::FeaturesCalculator(SCIP *scip, int signB, int signC, int signA) :
        nStaticFeatures(1 + (signC==0) + 2*(1+(signB==0)) + 2*(1+(signC==0)) + 2*(1+3*(signA==0))),
        nfeatures(nStaticFeatures),
        varnameIndexMap(),
        staticFeaturesMap(),
        dynamicFeaturesMap(),
        numberBrchMap(),
        nvars(SCIPgetNVars(scip)),
        nbrchs(0){
    SCIP_Var **vars = SCIPgetVars(scip);
    int nconss = SCIPgetNConss(scip);
    SCIP_Cons **conss = SCIPgetConss(scip);

    SCIP_Real *consCoefTemp;
    SCIPallocBlockMemoryArray(scip, &consCoefTemp, nvars);
    SCIP_Real *consCoef;
    SCIPallocBlockMemoryArray(scip, &consCoef, nvars);
    SCIP_Var **consVars;
    SCIPallocBlockMemoryArray(scip, &consVars, nvars);

    for(int i=0; i<nvars; ++i){
        varnameIndexMap[SCIPvarGetName(vars[i])] = i;
    }

    // compute features for evey variable i
    for(int i=0; i<nvars; ++i){
        unsigned featureIndex = 0;

        auto* features = new double[nStaticFeatures];
        std::fill(features, features+nStaticFeatures, 0);
        //SCIPallocBlockMemoryArray(scip, &features, nfeatures);

        //staticFeaturesMap.insert(std::pair<const char*, const double*>(SCIPvarGetName(vars[i]), features));
        std::string varName = "t_" + std::string(SCIPvarGetName(vars[i]));
        staticFeaturesMap[varName] = features;
        dynamicFeaturesMap[varName] = new double[nDynamicFeatures];
        objectiveIncreaseStaticsMap[varName] = new double[nObjectiveIncreaseStatics];
        std::fill(objectiveIncreaseStaticsMap[varName], objectiveIncreaseStaticsMap[varName]+nObjectiveIncreaseStatics, 0);

        numberBrchMap[varName] = 0;

        double objCoefVari= SCIPvarGetObj(vars[i]);// compute first set of static features |c_i| / sum_j(c_j) (2 features)

        computeSet1StaticFeatures(signC, vars, i, features, featureIndex, objCoefVari);

        // compute second set of features
        //init features array to compute min/max
        int sign = 1;
        for(int numMinMaxFeature=0; numMinMaxFeature<nStaticFeatures-featureIndex; ++numMinMaxFeature){
            features[featureIndex + numMinMaxFeature] = sign*DBL_MAX;
            sign *= -1;
        }

        SCIP_Bool success;
        for(int j=0; j<nconss; ++j){
            featureIndex = 1 + (signC==0);
            std::fill(consCoef, consCoef+nvars, 0);
            std::fill(consCoefTemp, consCoefTemp+nvars, 0);
            std::fill(consVars, consVars+nvars, nullptr);
            SCIPgetConsVals(scip, conss[j], consCoefTemp, nvars, &success);
            SCIPgetConsVars(scip, conss[j], consVars, nvars, &success);
            assert(success);

            for(int k=0; k<nvars; ++k){
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

            //SCIPdebugMessagePrint(scip, (varName + ", cons: " + std::to_string(j) + " -> " + std::to_string(features[3] == 1.7976931348623157e+308) + "\n").c_str());

            indexOffset += 2*(1+(signC==0));
            computeM3Features(signA, consCoef, i, featureIndex, features, indexOffset, localIndexOffset, val);
        }

        // Every unmodified field is set to 0
        sign = 1;
        for(int numMinMaxFeature=0; numMinMaxFeature<nStaticFeatures-featureIndex; ++numMinMaxFeature){
            if(features[featureIndex + numMinMaxFeature] == sign*DBL_MAX){
                features[featureIndex + numMinMaxFeature] = 0;
            }
            sign *= -1;
        }
    }

    SCIPfreeBlockMemoryArray(scip, &consCoef, nvars);
    SCIPfreeBlockMemoryArray(scip, &consCoefTemp, nvars);
    SCIPfreeBlockMemoryArray(scip, &consVars, nvars);
}

void FeaturesCalculator::computeSet1StaticFeatures(int signC, SCIP_Var *const *vars, int i, double *features,
                                                   unsigned int &featureIndex, double &objCoefVari) const {
    for(int j=0; j < nvars; ++j){
        double objCoef = SCIPvarGetObj(vars[j]);
        features[featureIndex + (objCoef<0 && signC==0)] += std::abs(objCoef);
    }
    for(int k: {0, 1}) {
        if(signC!=0 && k==1)continue;
        if(features[featureIndex])
            features[featureIndex] = std::abs(objCoefVari) / features[featureIndex];
        else
            features[featureIndex] = 0;
        ++featureIndex;
    }
    /////

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
    if(consCoef[i]!=0) {
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
    ////

}

void FeaturesCalculator::computeM1Features(int signB, const double *consCoef, int i, unsigned int featureIndex,
                                           double *features, double bj) const {// compute M1 (4 features)
// TODO: M1- utils if c_j >= 0 for all j ???? Maybe less features for this types of problems. Same reflexion for other features set
    int localIndexOffset=0;
    if(signB==0)
        localIndexOffset = 2 * (bj<0);

    if(bj!=0) {
        double val = consCoef[i] / std::abs(bj);

        if (val < features[featureIndex + localIndexOffset]) { // keep the min
            features[featureIndex + localIndexOffset] = val;
        }
        if (val > features[featureIndex + localIndexOffset + 1]) {//keep the max
            features[featureIndex + localIndexOffset + 1] = val;
        }
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

void FeaturesCalculator::updateBranchCounter(SCIP_NODE *node, SCIP_VAR *var) {
    double obj = SCIPnodeGetEstimate(node);
    double parentObj = SCIPnodeGetEstimate(SCIPnodeGetParent(node));
    double increase = (parentObj - obj)/obj;
    std::string key = std::string(SCIPvarGetName(var));
    numberBrchMap[key] = numberBrchMap[key] + 1 ;
    int n = numberBrchMap[SCIPvarGetName(var)];
    double* statistics = objectiveIncreaseStaticsMap[SCIPvarGetName(var)];

    if(increase < statistics[0] || n == 1){
        statistics[0] = increase;
    }
    if(increase > statistics[1] || n == 1){
        statistics[1] = increase;
    }

    if(n>1){
        // update of the mean
        double oldMean = statistics[2];
        statistics[2] = 1.0/n * (increase + 1.0/(n-1)*oldMean);

        // update of std
        statistics[3] = ((n-2)*statistics[3] + (n-1)*std::pow(oldMean-statistics[2], 2) + std::pow(increase-statistics[2], 2))/(n-1);
    } else{
        statistics[2] = increase;
        statistics[3] = 0;
    }
    nbrchs++;
}

void FeaturesCalculator::computeDynamicProblemFeatures(SCIP *scip) {
    SCIP_Var **vars = SCIPgetVars(scip);

    int nFixedVar = 0;
    for(int i=0; i<nvars; ++i){
        if(SCIPvarGetUbGlobal(vars[i]) == SCIPvarGetLbGlobal(vars[i]))
            ++nFixedVar;
    }

    for(int i=0; i<nvars; ++i){
        std::string key = std::string(SCIPvarGetName(vars[i]));
        auto* features = dynamicFeaturesMap[key];
        features[0] = (double)nFixedVar/nvars;
        double varObj = SCIPvarGetLPSol(vars[i]);
        features[1] = std::ceil(varObj) - varObj;
        if(nbrchs)
            features[2] = (double)numberBrchMap[key] / nbrchs;
        else
            features[2] = 0;
    }
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

