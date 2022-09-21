//
// Created by simon on 14/08/22.
//

#include "Utils.h"
#include "objscip/objscipdefplugins.h"
#include "nodesel_random.h"
#include "branch_unrealistic.h"

using namespace scip;

SCIP_Retcode Utils::create_scip_instance(SCIP** scipp, bool addBranchScheme) {
    *scipp = nullptr;

    SCIP_CALL(SCIPcreate(scipp));
    SCIP* scip = *scipp;

    /* we explicitly enable the use of a debug solution for this main SCIP instance */
    SCIPenableDebugSol(*scipp);

    /* include default SCIP plugins */
    SCIP_CALL( SCIPincludeDefaultPlugins(*scipp) );
    SCIP_CALL(SCIPincludeObjBranchrule(scip, new Branch_unrealistic(scip), TRUE));
    //SCIP_CALL( SCIPincludeObjNodesel(scip, new NodeSelRandom(scip), TRUE) );

    configure_scip_instance(scip, addBranchScheme);

    return SCIP_OKAY;
}

SCIP_Retcode Utils::configure_scip_instance(SCIP *scip, bool addBranchScheme) {
    if(addBranchScheme) {
        SCIP_CALL( SCIPsetIntParam(scip,"branching/unrealistic/priority",536870911) );
        SCIP_CALL( SCIPsetIntParam(scip,"branching/vanillafullstrong/priority",0) );
        SCIP_CALL( SCIPsetIntParam(scip,"branching/pscost/priority",0) );

    } else{
        SCIP_CALL( SCIPsetIntParam(scip,"branching/pscost/priority",536870911) );
        SCIP_CALL( SCIPsetIntParam(scip,"branching/unrealistic/priority",0) );
        SCIP_CALL( SCIPsetRealParam(scip,"limits/time",10) );
    }

    /* for column generation instances, disable restarts */
    SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );
    SCIP_CALL( SCIPsetIntParam(scip,"limits/restarts",0) );

    /* turn off all separation algorithms */
    SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
    SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
    SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );

    //SCIP_CALL( SCIPsetIntParam(scip,"nodeselection/random/stdpriority",900000) );
    //SCIP_CALL( SCIPsetIntParam(scip,"nodeselection/dfs/stdpriority",900000) );
    SCIP_CALL( SCIPsetIntParam(scip,"propagating/maxroundsroot",0) );
    //SCIP_CALL( SCIPsetIntParam(scip,"propagating/maxrounds",0) );
    SCIP_CALL( SCIPsetIntParam(scip,"display/freq",1) );
    SCIP_CALL( SCIPsetIntParam(scip,"lp/disablecutoff",1) );

    SCIP_CALL( SCIPsetIntParam(scip,"lp/solvefreq",1) );
    SCIP_CALL( SCIPsetIntParam(scip,"lp/threads",1) );
    //SCIP_CALL( SCIPsetLongintParam(scip,"lp/iterlim",1) );
    //SCIP_CALL( SCIPsetLongintParam(scip,"lp/rootiterlim",1) );

    SCIP_CALL( SCIPsetIntParam(scip,"separating/maxruns",0) );
    SCIP_CALL( SCIPsetIntParam(scip,"pricing/maxvars",1) );
    SCIP_CALL( SCIPsetIntParam(scip,"pricing/maxvarsroot",1) );

    SCIP_CALL( SCIPsetBoolParam(scip,"benders/copybenders",FALSE) );
    SCIP_CALL( SCIPsetBoolParam(scip,"benders/cutlpsol",FALSE) );
    return SCIP_OKAY;
}
