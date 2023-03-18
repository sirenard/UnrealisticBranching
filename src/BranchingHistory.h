//
// Created by simon on 17/03/23.
//

#ifndef LEARNING_BRANCHINGHISTORY_H
#define LEARNING_BRANCHINGHISTORY_H

#include <vector>
#include <scip/scip_var.h>

struct BranchingItem{
    int varIndex;
    SCIP_Real varValue;
};

class BranchingHistory: public std::vector<BranchingItem>{
public:
    explicit BranchingHistory();
    void fill(BranchingItem* items, int size);
    void addElement(SCIP_VAR* var);
};


#endif //LEARNING_BRANCHINGHISTORY_H
