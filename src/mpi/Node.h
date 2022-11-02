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
    const unsigned rank;
    MPI_Status status;
    SCIP* scip;

public:
    explicit Node(unsigned rank);
    static void setInstance(Node* node);
    static Node* getInstance();
    bool isMaster();

    const unsigned int getRank() const;

    virtual void run() = 0;
};


#endif //LEARNING_NODE_H
