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
    const int nNotreeFeatures = 25;
    const int nTreeFeatures = 58;

    double open_lbs_min, open_lbs_max, open_lbs_mean, open_lbs_std, open_lbs_squared_sum;
    double open_ds_min, open_ds_max, open_ds_mean, open_ds_std, open_ds_squared_sum;

    int nvars;
    int nbrchs;

    double varScore(double s_i, double s_avg);
    double gNormMax(double x);
    double relDist(double x, double y);
    double relPos(double z, double x, double y);

public:
    /**
     * Build Feature calculator and compute static features
     * @param scip
     * @param signB sign of elements in b vector. 1: all positive, -1 all negative, 0: positive and positive
     * @param signC sign of elements in v vector. 1: all positive, -1 all negative, 0: positive and positive
     * @param signA sign of elements in A matrix. 1: all positive, -1 all negative, 0: positive and positive
     */
    FeaturesCalculator(SCIP *scip, int signB, int signC, int signA);
    const std::vector<double> getFeatures(SCIP_Var *var, SCIP *scip);

    double *getNoTreeFeatures(SCIP *scip, SCIP_Var *var);
    double *getTreeFeatures(SCIP *scip, int nlpcands);
    void updateBranching(SCIP *scip);

    int getNoTreeFeaturesSize();

    int getTreeFeaturesSize();
};


#endif //LEARNING_FEATURESCALCULATOR_H
