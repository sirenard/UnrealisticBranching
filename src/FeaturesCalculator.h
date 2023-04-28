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
#include <dlib/matrix.h>

class FeaturesCalculator {
    const int nStaticFeatures;
    const int nDynamicFeatures = 5;
    const int nObjectiveIncreaseStatics = 5;

    double sumObjCoefs[2]; // sum of positive c_i and absolute sum if negative c_i
    int signA, signB, signC;

    std::map<std::string, int> varnameIndexMap; // map varindex to position
    double** staticFeatures;
    double** dynamicFeatures;
    int* numberBrchs;

    // one entry is composed of:
    // min, max, mean, std
    double** objectiveIncreaseStaticsMap;
    int nvars;
    int nbrchs; // number of branching

    void
    computeM1Features(const double *consCoef, int i, unsigned int featureIndex, double *features,
                      double bj) const;

    void computeM2Features(const double *consCoef, int i, unsigned int featureIndex,
                           double *features, double objCoefVari, double val, int indexOffset) const;

    void computeM3Features(const double *consCoef, int i, unsigned int featureIndex, double *features,
                           int indexOffset, int &localIndexOffset, double &val) const;

    void
    computeSet1StaticFeatures(double *features, unsigned int &featureIndex,
                              double &objCoefVari) const;

    int getVarKey(SCIP_Var *var);
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

    const double *getStaticFeatures(SCIP_VAR *var, SCIP *scip);
    const double *getDynamicProblemFeatures(SCIP_VAR *var);
    const double* getDynamicOptimizationFeatures(SCIP_VAR *var);

    const int getNStaticFeatures() const;

    const int getNDynamicFeatures() const;

    const int getNObjectiveIncreaseStatics() const;

    const int getNFeatures() const;

    void getFeatures(SCIP_Var *var, SCIP *scip, double *features);
    virtual ~FeaturesCalculator();

    /**
     * Compute the sensitivity range of the variables in vars. All var in vars must be in the basis of the simplex
     * @param scip
     * @param lb lb[i] = lower bound for the coefficient of variable vars[i]. Must be of size at least varsSize
     * @param ub ub[i] = upper bound for the coefficient of variable vars[i]. Must be of size at least varsSize
     * @param vars Array of vars
     * @param varsSize Size of vars
     */
    void computeSensitivity(SCIP *scip, double *lb, double *ub, SCIP_Var **vars, int varsSize);

    void computStaticFeatures(SCIP *scip, SCIP_Var *var);
};


#endif //LEARNING_FEATURESCALCULATOR_H
