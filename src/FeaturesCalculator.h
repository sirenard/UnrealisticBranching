//
// Created by simon on 11/09/22.
//

#ifndef LEARNING_FEATURESCALCULATOR_H
#define LEARNING_FEATURESCALCULATOR_H

#include <string>
#include <map>
#include <vector>
#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

class FeaturesCalculator {
    const int nStaticFeatures;
    const int nDynamicFeatures = 4;
    const int nObjectiveIncreaseStatics = 4;
    const int nfeatures;

    std::map<std::string , int> varnameIndexMap;
    std::map<std::string , double*> staticFeaturesMap;
    std::map<std::string , double*> dynamicFeaturesMap;
    std::map<std::string , int> numberBrchMap;

    // one entry is composed of:
    // min, max, mean, std
    std::map<std::string , double*> objectiveIncreaseStaticsMap;
    int nvars;
    int nbrchs; // number of branching

    void
    computeM1Features(int signB, const double *consCoef, int i, unsigned int featureIndex, double *features,
                      double bj) const;

    void computeM2Features(int signC, const double *consCoef, int i, unsigned int featureIndex,
                           double *features, double objCoefVari, double val, int indexOffset) const;

    void computeM3Features(int signA, const double *consCoef, int i, unsigned int featureIndex, double *features,
                           int indexOffset, int &localIndexOffset, double &val) const;

    void
    computeSet1StaticFeatures(int signC, SCIP_Var *const *vars, int i, double *features, unsigned int &featureIndex,
                              double &objCoefVari) const;
public:
    /**
     * Build Feature calculator and compute static features
     * @param scip
     * @param signB sign of elements in b vector. 1: all positive, -1 all negative, 0: positive and positive
     * @param signC sign of elements in v vector. 1: all positive, -1 all negative, 0: positive and positive
     * @param signA sign of elements in A matrix. 1: all positive, -1 all negative, 0: positive and positive
     */
    FeaturesCalculator(SCIP *scip, int signB, int signC, int signA);
    void updateBranchCounter(SCIP_NODE **nodes, SCIP_VAR *var);
    void computeDynamicProblemFeatures(SCIP *scip, SCIP_Var** vars, int varsSize);
    void computeDynamicProblemFeatures(SCIP *scip);
    const double* getStaticFeatures(SCIP_VAR *var);
    const double *getDynamicProblemFeatures(SCIP_VAR *var);
    const double* getDynamicOptimizationFeatures(SCIP_VAR *var);

    const int getNStaticFeatures() const;

    const int getNDynamicFeatures() const;

    const int getNObjectiveIncreaseStatics() const;

    const std::vector<double> getFeatures(SCIP_Var *var);
    virtual ~FeaturesCalculator();
};


#endif //LEARNING_FEATURESCALCULATOR_H
