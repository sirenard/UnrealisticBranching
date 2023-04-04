//
// Created by simon on 17/03/23.
//

#include "BranchingHistory.h"
#include <scip/scip.h>


BranchingHistory::BranchingHistory() : std::vector<BranchingItem>() {}

void BranchingHistory::addElement(SCIP_VAR *var) {
    int index = SCIPvarGetProbindex(var);
    double value = SCIPvarGetLPSol(var);
    push_back({index, value});
}

void BranchingHistory::fill(BranchingItem *items, int size) {
    clear();
    for(int i=0; i<size; i++) {
        push_back(items[i]);
    }
}
