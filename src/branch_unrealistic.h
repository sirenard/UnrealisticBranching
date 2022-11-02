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

    SCIP_DECL_BRANCHEXECLP(scip_execlp) override;

    SCIP_DECL_BRANCHINIT(scip_init) override;


    /**
     * Compute the realNnodes i.e. the minimum number of nodes needed to solve the 2 children problem if the branching
     * is done on a given variable
     * @param scip
     * @param score will be filled by the realNnodes of the variable
     * @param bestScore current best realNnodes seen so far, -1 if none
     * @param fracValue fractional value of the tested variable
     * @param varbrch variable to branch on
     * @return
     */
    SCIP_RETCODE computeScore(SCIP *scip, int &score, SCIP_Real *childPrimalBounds, int bestScore,
                              SCIP_Real fracValue, SCIP_VAR *varbrch, SCIP_BoundType &branchSide) const;

    /**
     * translate the best sol of scipmain to scip_copy
     */
    const SCIP_Retcode setBestSol(SCIP *scip, SCIP *scip_copy) const;

public:
    explicit Branch_unrealistic(SCIP *scip, int maxdepth=1, double leafTimeLimit=-1);
    int* getMaxDepthPtr();

    void setFirstBranch(SCIP_Var *firstBranch);

    static void setDataWriter(DatasetWriter *dataWriter);

    double *getLeafTimeLimitPtr();

    Branch_unrealistic(SCIP *scip, int depth, int maxdepth, double leafTimeLimit);
};


#endif //LEARNING_BRANCH_UNREALISTIC_H
