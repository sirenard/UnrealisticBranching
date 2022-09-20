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
    SCIP_DECL_BRANCHEXECLP(scip_execlp) override;

    /**
     * Compute the realNnodes i.e. the minimum number of nodes needed to solve the problem if the branching
     * is done on a given variable
     * @param scip
     * @param score will be filled by the realNnodes of the variable
     * @param branchSide will be filled by 0 (1) if the optimal solution is get by exploring the left (adding an upper bound) child
     * when branching on the given variable
     * @param bestScore current best realNnodes seen so far, -1 if none
     * @param fracValue fractional value of the tested variable
     * @param varbrch variable to branch on
     * @return
     */
    SCIP_RETCODE computeScore(SCIP *scip, int *score, SCIP_Real *childPrimalBounds, int bestScore,
                              SCIP_Real fracValue, SCIP_VAR *varbrch, SCIP_BoundType &branchSide) const;

    Branch_unrealistic(SCIP *scip, int depth, int maxdepth);

public:
    explicit Branch_unrealistic(SCIP *scip, int maxdepth=1);

    static void setDataWriter(DatasetWriter *dataWriter);
};


#endif //LEARNING_BRANCH_UNREALISTIC_H
