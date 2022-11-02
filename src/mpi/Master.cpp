//
// Created by simon on 26/10/22.
//

#include "Master.h"


Master::Master(int nslaves):
        Node(0),
        nslaves(nslaves){
}

void Master::run(){

}

void Master::computeScores(SCIP *scip, SCIP_VAR **lpcands, int nlpcands, std::vector<int> &bestcands, int &bestScore) {
    for(int c = 0; c < nlpcands; c += nslaves){
        // launch the tasks
        for(int s=1; s <= nslaves; ++s){
            int i = c + s - 1;
            if(i == nlpcands)break;
            sendNode(scip, s, bestScore==-1?-1:bestScore+1, lpcands[i]);
        }

        // retrieve result
        int score;
        for(int s=1; s <= nslaves; ++s){
            int i = c + s - 1;
            if(i == nlpcands)break;
            MPI_Recv(&score, 1, MPI_INT, s, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(score < bestScore || bestScore==-1){
                bestcands.clear();
                bestcands.push_back(i);
                bestScore = score;
            } else if(score == bestScore){
                bestcands.push_back(i);
            }
        }
    }
}

void Master::broadcastInstance(SCIP *scip) {
    SCIP_VAR** conssVars = SCIPgetVars(scip);
    int n = SCIPgetNVars(scip);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    VarInfo varInfo;
    for(int i=0; i<n; ++i){
        varInfo.type = SCIPvarGetType(conssVars[i]);
        varInfo.lb = SCIPvarGetLbGlobal(conssVars[i]);
        varInfo.ub = SCIPvarGetUbGlobal(conssVars[i]);
        varInfo.obj = SCIPvarGetObj(conssVars[i]);

        MPI_Bcast(&varInfo, sizeof(VarInfo), MPI_BYTE, 0, MPI_COMM_WORLD);
    }

    int m = SCIPgetNConss(scip);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    SCIP_CONS **conss = SCIPgetConss(scip);

    SCIP_Bool success;
    SCIP_Real *vals;
    int *consVarsIndex;
    SCIP_VAR **consVars;
    SCIPallocBlockMemoryArray(scip, &vals, n);
    SCIPallocBlockMemoryArray(scip, &consVars, n);
    SCIPallocBlockMemoryArray(scip, &consVarsIndex, n);

    // for each constraint
    ConsInfo consInfo;
    for(int j=0; j<m; ++j){
        consInfo.lhs = SCIPconsGetLhs(scip, conss[j], &success);
        consInfo.rhs = SCIPconsGetRhs(scip, conss[j], &success);

        // get the number of vars in the constraint
        SCIPgetConsNVars(scip, conss[j], &(consInfo.nVar), &success);
        MPI_Bcast(&consInfo, sizeof(ConsInfo), MPI_BYTE, 0, MPI_COMM_WORLD);

        SCIPgetConsVars(scip, conss[j], consVars, n, &success); // get the conssVars of the constraint
        SCIPgetConsVals(scip, conss[j], vals, n, &success); // get the coefs of the constraint
        for(int i = 0; i < consInfo.nVar; ++i){
            consVarsIndex[i] =  SCIPvarGetProbindex(consVars[i]);
        }

        MPI_Bcast(consVarsIndex, consInfo.nVar, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(vals, consInfo.nVar, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    SCIPfreeBlockMemoryArray(scip, &vals, n);
    SCIPfreeBlockMemoryArray(scip, &consVars, n);
    SCIPfreeBlockMemoryArray(scip, &consVarsIndex, n);

    MPI_Barrier(MPI_COMM_WORLD);
}

void Master::sendNode(SCIP *scip, int slaveId, int nodeLimit, SCIP_VAR *varbrch) {
    double leafTimeLimit;
    SCIPgetRealParam(scip, "branching/unrealistic/leaftimelimit", &leafTimeLimit);
    MPI_Send(&leafTimeLimit, 1, MPI_DOUBLE, slaveId, 1, MPI_COMM_WORLD);
    MPI_Send(&nodeLimit, 1, MPI_INT, slaveId, 1, MPI_COMM_WORLD);

    int n = SCIPgetNVars(scip);
    SCIP_VAR **vars = SCIPgetVars(scip);
    SCIP_Real *vals;
    SCIPallocBlockMemoryArray(scip, &vals, n);

    // send lower bounds
    for(int i=0; i<n; ++i){
        vals[i] = SCIPvarGetLbLocal(vars[i]);
    }
    MPI_Send(vals, n, MPI_DOUBLE, slaveId, 1, MPI_COMM_WORLD);

    // send upper bounds
    for(int i=0; i<n; ++i){
        vals[i] = SCIPvarGetUbLocal(vars[i]);
    }
    MPI_Send(vals, n, MPI_DOUBLE, slaveId, 1, MPI_COMM_WORLD);

    int id = SCIPvarGetProbindex(varbrch);
    MPI_Send(&id, 1, MPI_INT, slaveId, 1, MPI_COMM_WORLD);


    SCIPfreeBlockMemoryArray(scip, &vals, n);
}