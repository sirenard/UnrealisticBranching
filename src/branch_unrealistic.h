//
// Created by simon on 14/08/22.
//

#ifndef LEARNING_BRANCH_UNREALISTIC_H
#define LEARNING_BRANCH_UNREALISTIC_H



#include <string>
#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#include "objscip/objscip.h"
#include "DatasetWriter.h"

class Branch_unrealistic: public scip::ObjBranchrule {
private:
    static DatasetWriter *dataWriter;
    int depth, maxdepth;
    SCIP_Var* firstBranch;
    double leafTimeLimit;
    double left, right;

    SCIP_DECL_BRANCHEXECLP(scip_execlp) override;

    SCIP_DECL_BRANCHINIT(scip_init) override;


public:
    explicit Branch_unrealistic(SCIP *scip, int maxdepth=1, double leafTimeLimit=-1);
    int* getMaxDepthPtr();

    void setFirstBranch(SCIP_Var *firstBranch, double d, double d1);

    static void setDataWriter(DatasetWriter *dataWriter);

    double *getLeafTimeLimitPtr();

    Branch_unrealistic(SCIP *scip, int depth, int maxdepth, double leafTimeLimit);
};


#endif //LEARNING_BRANCH_UNREALISTIC_H
