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
        int nameLength;
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
                  int maxdepth, double leafTimeLimit);

    void broadcastInstance(const char *name);

    SCIP *
    sendNode(SCIP *scip, unsigned int workerId, int nodeLimit, SCIP_VAR *varbrch, int depth, int maxdepth, double objlimit, double leafTimeLimit);

    void sendWorkersRange(unsigned int workerRank, unsigned int start, unsigned int end);

    void getWorkersRange();

    SCIP *
    createScipInstance(double leafTimeLimit, int depth, int maxdepth, int nodeLimit, int n, double *lb, double *ub,
                       int firstBrchId, double objlimit);

    int getScore(SCIP *scip);

    void setScipInstance(SCIP *scip);

    unsigned int
    extractScore(SCIP_VAR *const *lpcands, std::vector<int> &bestcands, int depth, const int *workerMap,
                 int &bestScore) const;
};


#endif //LEARNING_WORKER_H
