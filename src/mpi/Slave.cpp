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
    std::cout << "Child init done." << std::endl;

    while (true) {
        SCIP *scip_copy = retrieveNode();
        SCIPsolve(scip_copy);
        int score = SCIPgetNNodes(scip_copy);
        std::cout << "Score: " << score << " - " << SCIPgetSolOrigObj(scip_copy, SCIPgetBestSol(scip_copy)) << std::endl;

        returnScore(score);
        SCIPfree(&scip_copy);
    }
}

void Slave::returnScore(int score) {
    MPI_Send(&score, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
}

void Slave::retrieveInstance() {
    int n; // number of variables in the problem
    int m; // number of constraints in the problem
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    SCIP_Real *objCoef;
    SCIPallocBlockMemoryArray(scipmain, &objCoef, n);

    // Receive variables data
    VarInfo varInfo;
    for(int i = 0; i < n; ++i){
        SCIP_Var* var;

        MPI_Bcast(&varInfo, sizeof(VarInfo), MPI_BYTE, 0, MPI_COMM_WORLD);

        SCIPcreateVar(scipmain, &var, NULL, varInfo.lb, varInfo.ub, varInfo.obj, varInfo.type, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL);
        SCIPaddVar(scipmain, var);
        SCIPreleaseVar(scipmain, &var);
    }

    SCIP_VAR ** vars = SCIPgetVars(scipmain);
    // Receive constraints data
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

    SCIP_Bool success;
    SCIP_Real *vals;
    SCIP_VAR **consVars;
    int *consVarsIndex;
    SCIPallocBlockMemoryArray(scipmain, &vals, n);
    SCIPallocBlockMemoryArray(scipmain, &consVars, n);
    SCIPallocBlockMemoryArray(scipmain, &consVarsIndex, n);



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
                scipmain,
                &cons,
                (std::to_string(j)+"cons").c_str(),
                consInfo.nVar,
                consVars,
                vals,
                consInfo.lhs,
                consInfo.rhs
        );

        SCIPaddCons(scipmain, cons);
        SCIPreleaseCons(scipmain, &cons);
    }


    SCIPfreeBlockMemoryArray(scipmain, &vals, n);
    SCIPfreeBlockMemoryArray(scipmain, &consVars, n);
    SCIPfreeBlockMemoryArray(scipmain, &consVarsIndex, n);
    SCIPfreeBlockMemoryArray(scipmain, &objCoef, n);

    MPI_Barrier(MPI_COMM_WORLD);
    /*double objLimit;
    MPI_Recv(&objLimit, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    SCIPsetObjlimit(scipmain, objLimit);

    int varToBranchIndex;
    MPI_Recv(&varToBranchIndex, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    SCIP_VAR* varTobranch = vars[varToBranchIndex];

    // TODO: Slave stay in UNREALISTIC branch rule
    Branch_unrealistic *rule = reinterpret_cast<Branch_unrealistic *>(SCIPfindBranchrule(scipmain, "unrealistic"));
    rule->setFirstBranch(varTobranch);*/

}

Slave::Slave(unsigned int rank) : Node(rank) {
}

void Slave::createScipInstance() {
    Utils::create_scip_instance(&scipmain, true);
    //Utils::configure_slave_scip_instance(scipmain);
    SCIPcreateProbBasic(scipmain, ("slaveProb" + std::to_string(rank)).c_str());
}

SCIP *Slave::retrieveNode() {
    SCIP *scip_copy;
    SCIP_Bool valid;
    SCIPcreate(&scip_copy);
    SCIP_HashMap *varmap;
    SCIPhashmapCreate(&varmap, SCIPblkmem(scip_copy), SCIPgetNVars(scipmain));
    SCIPcopy(scipmain, scip_copy, varmap, nullptr, "ub_subsolver", FALSE, FALSE, FALSE, FALSE, &valid);
    assert(valid == TRUE);

    SCIPcopyParamSettings(scipmain, scip_copy);

    double leafTimeLimit;

    MPI_Recv(&leafTimeLimit, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    SCIPsetRealParam(scipmain, "branching/unrealistic/leaftimelimit", leafTimeLimit);
    auto *objbranchrule = new Branch_unrealistic(scip_copy, 1, 1, leafTimeLimit);

    int nodeLimit;
    MPI_Recv(&nodeLimit, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    SCIPincludeObjBranchrule(scip_copy, objbranchrule, TRUE);
    Utils::configure_scip_instance(scip_copy, true);
    SCIPsetLongintParam(scip_copy, "limits/nodes", nodeLimit);
    SCIPsetIntParam(scip_copy, "display/verblevel",0);
    // don't use heuristics on recursion levels
    SCIPsetHeuristics(scip_copy, SCIP_PARAMSETTING_OFF, TRUE);
    SCIPsetRealParam(scip_copy,"limits/time",600);

    // TODO: Is usefull?
    SCIPmergeVariableStatistics(scipmain, scip_copy, SCIPgetVars(scipmain), SCIPgetVars(scip_copy), SCIPgetNVars(scipmain));

    int n = SCIPgetNVars(scip_copy);
    SCIP_VAR **vars = SCIPgetVars(scip_copy);
    SCIP_Real *vals;
    SCIPallocBlockMemoryArray(scip_copy, &vals, n);
    // get lower bounds
    MPI_Recv(vals, n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    for(int i=0; i<n; ++i){
        SCIPchgVarLb(scip_copy, vars[i], vals[i]);
    }


    // get upper bounds
    MPI_Recv(vals, n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    for(int i=0; i<n; ++i){
        SCIPchgVarUb(scip_copy, vars[i], vals[i]);
    }

    int id;
    MPI_Recv(&id, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    //SCIP_VAR* varbrch_copy = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[id]);
    //assert(varbrch != nullptr);
    objbranchrule->setFirstBranch(vars[id]);

    SCIPfreeBlockMemoryArray(scip_copy, &vals, n);
    return scip_copy;
}
