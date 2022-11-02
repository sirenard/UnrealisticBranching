//
// Created by simon on 26/10/22.
//

#include <scip/scipdefplugins.h>
#include "Slave.h"
#include "../Utils.h"
#include "../branch_unrealistic.h"

void Slave::run() {
    createScipInstance();
    retrieveInstance();
    SCIPsolve(scip);

    while (true) {}

    /*while (true) {
        retrieveInstance();
        SCIPsolve(scip);
        int score = SCIPgetNNodes(scip);
        std::cout << score << " - " << SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip)) << std::endl;
        //SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, 0);

        returnScore(score);
        SCIPfree(&scip);
    }*/
}

void Slave::returnScore(int score) {
    MPI_Send(&score, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
}

void Slave::retrieveInstance() {
    int n; // number of variables in the problem
    int m; // number of constraints in the problem
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    SCIP_Real *objCoef;
    SCIPallocBlockMemoryArray(scip, &objCoef, n);

    // Receive variables data
    VarInfo varInfo;
    for(int i = 0; i < n; ++i){
        SCIP_Var* var;

        MPI_Bcast(&varInfo, sizeof(VarInfo), MPI_BYTE, 0, MPI_COMM_WORLD);

        SCIPcreateVar(scip, &var, NULL, varInfo.lb, varInfo.ub, varInfo.obj, varInfo.type, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL);
        SCIPaddVar(scip, var);
        SCIPreleaseVar(scip, &var);
    }

    SCIP_VAR ** vars = SCIPgetVars(scip);
    // Receive constraints data
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

    SCIP_Bool success;
    SCIP_Real *vals;
    SCIP_VAR **consVars;
    int *consVarsIndex;
    SCIPallocBlockMemoryArray(scip, &vals, n);
    SCIPallocBlockMemoryArray(scip, &consVars, n);
    SCIPallocBlockMemoryArray(scip, &consVarsIndex, n);



    ConsInfo consInfo;
    for(int j=0; j<m; ++j){
        MPI_Bcast(&consInfo, sizeof(ConsInfo), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(consVarsIndex, consInfo.nVar, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(vals, consInfo.nVar, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for(int i = 0; i < consInfo.nVar; ++i){
            int varNumber = consVarsIndex[i];
            consVars[i] = vars[varNumber];
        }

        SCIP_CONS *cons;
        SCIPcreateConsBasicLinear(
            scip,
            &cons,
            (std::to_string(j)+"cons").c_str(),
            consInfo.nVar,
            consVars,
            vals,
            consInfo.lhs,
            consInfo.rhs
        );

        SCIPaddCons(scip, cons);
        SCIPreleaseCons(scip, &cons);
    }

    SCIPfreeBlockMemoryArray(scip, &vals, n);
    SCIPfreeBlockMemoryArray(scip, &consVars, n);
    SCIPfreeBlockMemoryArray(scip, &consVarsIndex, n);
    SCIPfreeBlockMemoryArray(scip, &objCoef, n);

    /*double objLimit;
    MPI_Recv(&objLimit, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    SCIPsetObjlimit(scip, objLimit);

    int varToBranchIndex;
    MPI_Recv(&varToBranchIndex, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    SCIP_VAR* varTobranch = vars[varToBranchIndex];

    // TODO: Slave stay in UNREALISTIC branch rule
    Branch_unrealistic *rule = reinterpret_cast<Branch_unrealistic *>(SCIPfindBranchrule(scip, "unrealistic"));
    rule->setFirstBranch(varTobranch);*/

}

Slave::Slave(unsigned int rank) : Node(rank) {
}

void Slave::createScipInstance() {
    Utils::create_scip_instance(&scip, false);
    //Utils::configure_slave_scip_instance(scip);
    SCIPcreateProbBasic(scip, ("slaveProb" + std::to_string(rank)).c_str());
}
