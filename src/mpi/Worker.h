//
// Created by simon on 26/10/22.
//

#ifndef LEARNING_WORKER_H
#define LEARNING_WORKER_H


#include <vector>
#include <mpi.h>
#include <scip/scip.h>

class Worker {
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

    static Worker* instance;

    void returnScore(int score);

    void retrieveInstance();

    // a worker can use the workers in the range [startWorkersRange, endWorkersRange]
    unsigned startWorkersRange;
    unsigned endWorkersRange;
    unsigned nWorkers;
    unsigned directorRank;
    const unsigned rank;

    MPI_Status status;
    SCIP* scipmain;

public:
    explicit Worker(unsigned rank);
    static void setInstance(Worker* node);
    static Worker* getInstance();
    bool isMaster();

    const unsigned int getRank() const;

    void work();

    void setWorkersRange(int start, int end);

    void createScipInstance();

    SCIP *retrieveNode();

    void
    computeScores(SCIP *scip, SCIP_VAR **lpcands, int nlpcands, std::vector<int> &bestcands, int &bestScore, int depth,
                  int maxdepth);

    void broadcastInstance(SCIP *scip);

    SCIP * sendNode(SCIP *scip, unsigned int workerId, int nodeLimit, SCIP_VAR *varbrch, int depth, int maxdepth);

    void sendWorkersRange(unsigned int id, unsigned int start, unsigned int end);

    void getWorkersRange();

    SCIP *
    createScipInstance(double leafTimeLimit, int depth, int maxdepth, int nodeLimit, int n, double *lb, double *ub,
                       int firstBrchId);

    int getScore(SCIP *scip);

    void setScipInstance(SCIP *scip);
};


#endif //LEARNING_WORKER_H