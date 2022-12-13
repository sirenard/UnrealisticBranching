//
// Created by simon on 26/10/22.
//

#ifndef LEARNING_WORKER_H
#define LEARNING_WORKER_H

#define CONTINUE_WORKING_FLAG 0
#define END_SOLVING_FLAG 1
#define INSTRUCTION_FLAG 2
#define SCORE_FLAG 3
#define NODE_INFO_FLAG 4
#define WORKER_RANGE_FLAG 5

#define MIN_NUMBER_WORKER 10

#include <vector>
#include <mpi.h>
#include <scip/scip.h>

class Worker {
    static Worker* instance;

    void returnScore(int score);

    // a worker can use the workers in the range [startWorkersRange, endWorkersRange]
    unsigned nWorkers;
    unsigned* workers;
    unsigned directorRank;
    const unsigned rank;

    MPI_Status status;

    bool isFinished();

    /**
     * return an array of at most n free workers to work on task j
     * @param n
     * @param task
     * @return
     */
    unsigned* findAvailableWorkers(unsigned& n, int task, int *workerMap);
    void updateWork(unsigned workerRank, int task, int* workerMap);

    void transferWorkers(int fromTask, int *workerMap, int n);

    /**
     * Find the index i s.t. workers[i] = workerRank
     * @param workerRank
     * @return
     */
    unsigned findWorkerIndex(unsigned workerRank) const;

    static void sortVarsArray(SCIP_VAR** vars, int n);
    static int findVar(SCIP_VAR *var, SCIP_VAR **vars, int n);

public:
    explicit Worker(unsigned rank);
    static void setInstance(Worker* node);
    static Worker* getInstance();
    bool isMaster();

    void work();

    virtual ~Worker();

    void setWorkersRange(int start, int end);

    SCIP *retrieveNode();

    void
    computeScores(SCIP *scip, SCIP_VAR **lpcands, int nlpcands, std::vector<int> &bestcands, int &bestScore, int depth,
                  int maxdepth, double leafTimeLimit, bool noNodeLimitation, int *varScores);


    SCIP *
    sendNode(SCIP *scip, unsigned int workerId, int nodeLimit, SCIP_VAR *varbrch, int depth, int maxdepth, double objlimit, double leafTimeLimit);

    void sendWorkersRange(unsigned int workerRank, unsigned int *workersToSend, unsigned int n);

    void getWorkersRange();

    SCIP *
    createScipInstance(double leafTimeLimit, int depth, int maxdepth, int nodeLimit, int firstBrchId, double left,
                       double right, double objlimit, int branchingMaxDepth, double *bestSolVals,
                       const char *filename);

    int getScore(SCIP *scip);

    unsigned int
    extractScore(SCIP_VAR *const *lpcands, std::vector<int> &bestcands, int depth, const int *workerMap,
                 int &bestScore, int *varScores) const;

    void broadcastEnd();
};


#endif //LEARNING_WORKER_H
