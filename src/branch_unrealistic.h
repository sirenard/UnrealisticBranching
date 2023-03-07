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
    int depth;
    int maxdepth;
    double leafTimeLimit;

    std::vector<int>* branchingHistory;
    std::vector<double>* branchingHistoryValues;
    int branching_count;

    char scoreMethod='c';
    double alpha=0;
    double epsilon=0; // diversification factor

    SCIP_DECL_BRANCHEXECLP(scip_execlp) override;

    /**
     * Branch by following the branchingHistory vector. Used to make a copy of a tree.
     */
    SCIP_RETCODE branchCopycat(SCIP *scip, SCIP_RESULT *result);

    /**
     * Apply the unrealistic branching. Branch on the variable with less nodes.
     * If a datawriter is given, use it to store new features.
     */
    SCIP_RETCODE branchUnrealistic(SCIP *scip, SCIP_RESULT *result);


    SCIP_DECL_BRANCHEXIT(scip_exit) override;
public:
    explicit Branch_unrealistic(SCIP *scip, int maxdepth=1, double leafTimeLimit=-1);
    int* getMaxDepthPtr();

    static void setDataWriter(DatasetWriter *dataWriter);

    double *getLeafTimeLimitPtr();

    char* getScoreMethodPtr();

    double* getAlphaPtr();

    double * getEpsPtr();

    void fillBranchHistory(int *history, double *values, int size);

    void setLeafTimeLimit(double leafTimeLimit);

    void setDepth(int depth);

    ~Branch_unrealistic() override;

    std::vector<int> *getHistory();

    std::vector<double> *getBranchingHistoryValues() const;
};


#endif //LEARNING_BRANCH_UNREALISTIC_H
