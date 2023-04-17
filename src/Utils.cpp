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
#include "EventhdlrUpdateFeatures.h"

using namespace scip;

SCIP_Retcode Utils::create_scip_instance(SCIP **scipp) {
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

    SCIP_CALL(SCIPincludeObjEventhdlr(scip, new EventhdlrUpdateFeatures(scip), TRUE) );

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

    char possibilities[] = {'c', 'a', 'n'};
    SCIP_CALL(SCIPaddCharParam(
            scip,
            "dataset/scoremethod",
            "How the score must be computed. 'c'ontinuous mapping between 0 and 1. Or use of 'a'lpha: the score is 1 "
            "if variable SB score is <= (1+alpha)bestUbScore. Or 'n'one: no score computed for each features vector"
            "the current SB score, the minimum and maximum SB score are reported",
            objbranchrule->getScoreMethodPtr(),
            FALSE,
            'c',
            possibilities,
            NULL,
            NULL
    ));


    SCIP_CALL(SCIPaddRealParam(
            scip,
            "dataset/alpha",
            "If the scoremethod is 'a'",
            objbranchrule->getAlphaPtr(),
            FALSE,
            0.2,
            0,
            SCIP_REAL_MAX,
            NULL,
            NULL
    ));

    SCIP_CALL(SCIPaddRealParam(
            scip,
            "dataset/epsilon",
            "Random branching is performed with probability epsilon at each step during dataset generation",
            objbranchrule->getEpsPtr(),
            FALSE,
            0.3,
            0,
            1,
            NULL,
            NULL
    ));

    SCIP_CALL(SCIPaddRealParam(
            scip,
            "branching/unrealisticTrained/threshold",
            "Above which tthreshold does the predicted score is relevant",
            objbranchruleTrained->getAlphaPtr(),
            FALSE,
            0,
            0,
            1,
            NULL,
            NULL
    ));

    configure_scip_instance(scip);

    return SCIP_OKAY;
}

SCIP_Retcode Utils::configure_scip_instance(SCIP *scip) {
    SCIP_CALL( SCIPsetIntParam(scip,"branching/unrealistic/priority",100000) );
    SCIP_CALL( SCIPsetIntParam(scip,"branching/vanillafullstrong/priority",50000) );
    SCIP_CALL( SCIPsetIntParam(scip,"display/freq",1) );
    SCIP_CALL( SCIPsetRealParam(scip,"limits/time",1e+20) );


    /* turn off all separation algorithms */
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE));

    SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );

    SCIP_CALL(SCIPsetIntParam(scip, "lp/threads", 1));
    SCIP_CALL(SCIPsetCharParam(scip, "estimation/restarts/restartpolicy", 'n'));

    SCIP_CALL(SCIPsetIntParam(scip, "randomization/lpseed", 12));
    SCIP_CALL(SCIPsetIntParam(scip, "randomization/permutationseed", 12));
    SCIP_CALL(SCIPsetIntParam(scip, "randomization/randomseedshift", 12));
    SCIP_CALL(SCIPsetIntParam(scip, "randomization/randomseedshift", 12));
    SCIPinitializeRandomSeed(scip, 12);

    return SCIP_OKAY;
}

SCIP_Retcode Utils::congigure_scip_end_recursion(SCIP *scip, double leafTimeLimit) {
    // set all branch rules
    int nBranchingRules = SCIPgetNBranchrules(scip);
    SCIP_Branchrule** branchrules = SCIPgetBranchrules(scip);

    for(int i=0; i<nBranchingRules; ++i){
        int priority = SCIPbranchruleGetPriority(branchrules[i]);
        if(priority>50000) {
            const char *name = SCIPbranchruleGetName(branchrules[i]);
            std::string param = "branching/" + std::string(name) + "/priority";
            SCIPsetIntParam(scip, param.c_str(), 0);
        }
    }

    SCIPsetRealParam(scip, "limits/time", leafTimeLimit);

    return SCIP_OKAY;
}

SCIP_Branchrule ** Utils::get_branching_rules(SCIP *scip) {
    SCIP_Branchrule** branchrules = SCIPgetBranchrules(scip);
    int nBranchingRule = SCIPgetNBranchrules(scip);
    SCIP_Branchrule** res = new SCIP_Branchrule*[nBranchingRule];
    memcpy(res, branchrules, nBranchingRule* sizeof(SCIP_Branchrule*));

    std::sort(res, res+nBranchingRule,
         [](SCIP_Branchrule* & a, SCIP_Branchrule* & b) -> bool
         {
             return SCIPbranchruleCompName(a, b)>0;
         });

    return res;
}
