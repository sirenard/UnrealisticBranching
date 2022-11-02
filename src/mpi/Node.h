//
// Created by simon on 26/10/22.
//

#ifndef LEARNING_NODE_H
#define LEARNING_NODE_H


#include <mpi.h>
#include <scip/scip.h>

class Node {
    static Node* instance;
protected:
    struct VarInfo{
        SCIP_Vartype type;
        double lb;
        double ub;
        double obj;
    };

    struct ConsInfo{
        double lhs;
        double rhs;
        int nVar;
    };
    MPI_Status status;
    SCIP* scipmain;

public:
    explicit Node(unsigned rank);
    static void setInstance(Node* node);
    static Node* getInstance();
    bool isMaster();

    const unsigned int getRank() const;

    virtual void run() = 0;

    unsigned rank;
};


#endif //LEARNING_NODE_H
