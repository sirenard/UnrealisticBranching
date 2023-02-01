//
// Created by simon on 14/08/22.
//

#include "Utils.h"
#include "objscip/objscipdefplugins.h"
#include "nodesel_random.h"
#include "branch_unrealistic.h"
#include "dialog_generateDataset.h"
#include "branch_unrealisticTrained.h"
#include "dialog_loadModel.h"
#include "dialog_trainmodel.h"

using namespace scip;

SCIP_Retcode Utils::create_scip_instance(SCIP** scipp, bool addBranchScheme) {
    *scipp = nullptr;

    SCIP_CALL(SCIPcreate(scipp));
    SCIP* scip = *scipp;

    /* we explicitly enable the use of a debug solution for this main SCIP instance */
    SCIPenableDebugSol(*scipp);

    /* include default SCIP plugins */
    SCIP_CALL( SCIPincludeDefaultPlugins(*scipp) );

    Branch_unrealistic *objbranchrule = new Branch_unrealistic(scip);
    SCIP_CALL(SCIPincludeObjBranchrule(scip, objbranchrule, TRUE));

    Branch_unrealisticTrained *objbranchruleTrained = new Branch_unrealisticTrained(scip);
    SCIP_CALL(SCIPincludeObjBranchrule(scip, objbranchruleTrained, TRUE));

    SCIP_CALL(SCIPincludeObjDialog(scip, new DialogGenerateDataset(scip), TRUE));
    SCIP_CALL(SCIPincludeObjDialog(scip, new DialogLoadModel(scip), TRUE));
    SCIP_CALL(SCIPincludeObjDialog(scip, new DialogTrainModel(scip), TRUE));
    SCIP_CALL(SCIPaddIntParam(
            scip,
            "branching/unrealistic/recursiondepth",
            "How depth the unrealistic branching is performed with the recursion, -1 for unlimited depth",
            objbranchrule->getMaxDepthPtr(),
            FALSE,
            1,
            -1,
            INT_MAX,
            NULL,
            NULL
    ));

    SCIP_CALL(SCIPaddRealParam(
            scip,
            "branching/unrealistic/leaftimelimit",
            "Time limit allowed to the leaf of the recursion tree (not solved using the unrealistic branching). -1 for no limit ",
            objbranchrule->getLeafTimeLimitPtr(),
            FALSE,
            120,
            0,
            1e+20,
            NULL,
            NULL
    ));

    configure_scip_instance(scip, addBranchScheme);

    return SCIP_OKAY;
}

SCIP_Retcode Utils::configure_scip_instance(SCIP *scip, bool addBranchScheme) {
    if(addBranchScheme) {
        SCIP_CALL( SCIPsetIntParam(scip,"branching/unrealistic/priority",536870911) );
        SCIP_CALL( SCIPsetIntParam(scip,"branching/vanillafullstrong/priority",536870900) );
        SCIP_CALL( SCIPsetIntParam(scip,"display/freq",1) );
        SCIP_CALL( SCIPsetRealParam(scip,"limits/time",1e+20) );

    } else{
        SCIP_CALL( SCIPsetIntParam(scip,"branching/vanillafullstrong/priority",536870911) );
        SCIP_CALL( SCIPsetIntParam(scip,"branching/unrealistic/priority",0) );
        SCIP_CALL( SCIPsetIntParam(scip,"display/freq",100) );
        SCIP_CALL(     SCIPsetIntParam(scip, "display/verblevel", 0));
        //SCIP_CALL( SCIPsetRealParam(scipmain,"limits/time",1) );
    }

    /* for column generation instances, disable restarts */
    //SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrestarts", 0));
    //SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", 0));
    //SCIP_CALL(SCIPsetIntParam(scip, "limits/restarts", 0));
    //SCIP_CALL( SCIPsetBoolParam(scip, "lp/presolving", FALSE) );
    //SCIP_CALL(SCIPsetIntParam(scip, "lp/solutionpolishing", 0));
    //SCIP_CALL(SCIPsetIntParam(scip, "lp/scaling", 0));



    /* turn off all separation algorithms */
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE));

    SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );

    //SCIP_CALL( SCIPsetIntParam(scip,"nodeselection/random/stdpriority",900000) );
    //SCIP_CALL( SCIPsetIntParam(scip,"nodeselection/dfs/stdpriority",9000000) );
    //SCIP_CALL(SCIPsetIntParam(scip, "propagating/maxroundsroot", 0));
    //SCIP_CALL( SCIPsetIntParam(scip,"propagating/maxrounds",0) );
    //SCIP_CALL(SCIPsetIntParam(scip, "lp/disablecutoff", 1));

    //SCIP_CALL(SCIPsetIntParam(scip, "lp/solvefreq", 1));
    //SCIP_CALL(SCIPsetIntParam(scip, "lp/threads", 1));
    //SCIP_CALL( SCIPsetLongintParam(scip,"lp/iterlim",1) );
    //SCIP_CALL( SCIPsetLongintParam(scip,"lp/rootiterlim",1) );

    //SCIP_CALL(SCIPsetIntParam(scip, "separating/maxruns", 0));
    //SCIP_CALL(SCIPsetIntParam(scip, "pricing/maxvars", 1));
    //SCIP_CALL(SCIPsetIntParam(scip, "pricing/maxvarsroot", 1));

    //SCIP_CALL(SCIPsetBoolParam(scip, "benders/copybenders", FALSE));
    //SCIP_CALL(SCIPsetBoolParam(scip, "benders/cutlpsol", FALSE));

    /*SCIP_CALL(SCIPsetBoolParam(scip, "randomization/permuteconss", FALSE));
    SCIP_CALL(SCIPsetIntParam(scip, "randomization/lpseed", 12));
    SCIP_CALL(SCIPsetIntParam(scip, "randomization/permutationseed", 12));
    SCIP_CALL(SCIPsetIntParam(scip, "randomization/randomseedshift", 12));
    SCIP_CALL(SCIPsetIntParam(scip, "branching/random/seed", 12));*/

    //SCIP_CALL( SCIPsetBoolParam(scip, "branching/vanillafullstrong/idempotent", TRUE) );
    //SCIP_CALL(SCIPsetCharParam(scip, "nodeselection/childsel", 'd'));

    //srand(12);

    return SCIP_OKAY;
}

SCIP_Retcode Utils::remove_handlers(SCIP *scip) {
    /*for(const std::string hdlrName:{"linear", "integral"}){
        SCIP_CALL(SCIPsetIntParam(scip, (std::string("constraints/")+hdlrName+"/propfreq").c_str(), -1));
        SCIP_CALL(SCIPsetIntParam(scip, (std::string("constraints/")+hdlrName+"/sepafreq").c_str(), -1));
    }*/
    return SCIP_OKAY;
}
