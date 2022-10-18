//
// Created by simon on 18/10/22.
//

#ifndef LEARNING_BRANCH_UNREALISTICTRAINED_H
#define LEARNING_BRANCH_UNREALISTICTRAINED_H


#include <string>
#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#include "objscip/objscip.h"
#include "DatasetWriter.h"
#include "RegressionModel.h"

class Branch_unrealisticTrained: public scip::ObjBranchrule {
private:
    static FeaturesCalculator *featuresCalculator;
    static RegressionModel *model;

    SCIP_DECL_BRANCHEXECLP(scip_execlp) override;

public:
    explicit Branch_unrealisticTrained(SCIP *scip);

    static void setFeaturesCalculator(FeaturesCalculator *newCalculator);

    static void setModel(RegressionModel *newModel);
};

#endif //LEARNING_BRANCH_UNREALISTICTRAINED_H
