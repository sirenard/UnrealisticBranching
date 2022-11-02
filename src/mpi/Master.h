//
// Created by simon on 26/10/22.
//

#ifndef LEARNING_MASTER_H
#define LEARNING_MASTER_H

#include <vector>
#include "Node.h"

// Singleton master
class Master: public Node {
    int nslaves;
public:
    Master(int nslaves);

    void run() override;

    void computeScores(SCIP *scip, SCIP_VAR **lpcands, int nlpcands, std::vector<int> &bestcands, int &bestScore);

    void broadcastInstance(SCIP *scip);

    void sendNode(SCIP *scip, int slaveId, int nodeLimit, SCIP_VAR *varbrch);
};


#endif //LEARNING_MASTER_H
