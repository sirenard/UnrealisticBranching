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

    SCIP_DECL_BRANCHEXECLP(scip_execlp) override;
    SCIP_DECL_BRANCHEXIT(scip_exit) override;
public:
    explicit Branch_unrealistic(SCIP *scip, int maxdepth=1, double leafTimeLimit=-1);
    int* getMaxDepthPtr();

    static void setDataWriter(DatasetWriter *dataWriter);

    double *getLeafTimeLimitPtr();

    void fillBranchHistory(int *history, double *values, int size);

    void setLeafTimeLimit(double leafTimeLimit);

    void setDepth(int depth);

    ~Branch_unrealistic() override;

    std::vector<int> *getHistory();

    std::vector<double> *getBranchingHistoryValues() const;
};


#endif //LEARNING_BRANCH_UNREALISTIC_H
