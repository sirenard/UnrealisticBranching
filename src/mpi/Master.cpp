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
            //broadcastInstance(s, scip, lpcands[i]);
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

        SCIPgetConsVars(scip, conss[j], conssVars, n, &success); // get the conssVars of the constraint
        SCIPgetConsVals(scip, conss[j], vals, n, &success); // get the coefs of the constraint
        for(int i = 0; i < consInfo.nVar; ++i){
            consVarsIndex[i] =  SCIPvarGetProbindex(conssVars[i]);
        }

        MPI_Bcast(consVarsIndex, consInfo.nVar, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(vals, consInfo.nVar, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    }
    SCIPfreeBlockMemoryArray(scip, &vals, n);
    SCIPfreeBlockMemoryArray(scip, &consVars, n);
    SCIPfreeBlockMemoryArray(scip, &consVarsIndex, n);

    MPI_Barrier(MPI_COMM_WORLD);
    /*SCIP_Real objLimit;
    // get values of the best sol
    SCIP_Sol* bestSol = SCIPgetBestSol(scip);
    if(bestSol == NULL) {
        objLimit = SCIPgetObjlimit(scip);
    } else{
        objLimit = SCIPsolGetOrigObj(bestSol);
    }
    MPI_Send(&objLimit, 1, MPI_DOUBLE, slaveNumber, 1, MPI_COMM_WORLD);

    int varToBranchIndex = SCIPvarGetProbindex(varToBranch);
    MPI_Send(&varToBranchIndex, 1, MPI_INT, slaveNumber, 1, MPI_COMM_WORLD);*/
}
