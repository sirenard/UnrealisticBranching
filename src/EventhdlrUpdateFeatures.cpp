//
// Created by simon on 17/03/23.
//

#include <iostream>
#include "EventhdlrUpdateFeatures.h"

EventhdlrUpdateFeatures::EventhdlrUpdateFeatures(SCIP *scip):
ObjEventhdlr(scip, EVENT_HDLR_UPDATE_FEATURES_NAME,"event handler for new focused node, update if it exists, the data for features computation"),
featureCalculator(nullptr),
history(nullptr),
nextOneIsExploration(false)
{
}

SCIP_DECL_EVENTINITSOL(EventhdlrUpdateFeatures::scip_initsol){
    SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL) );
    SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_NODEBRANCHED, eventhdlr, NULL, NULL) );
    return SCIP_OKAY;
}


SCIP_DECL_EVENTEXEC(EventhdlrUpdateFeatures::scip_exec){
    if(SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEBRANCHED && (history || featureCalculator)) {
        int nchildren = SCIPgetNChildren(scip);
        assert(nchildren == 2);
        SCIP_NODE **children;
        SCIPgetChildren(scip, &children, NULL);

        int depth = 1;

        SCIP_Var **vars = new SCIP_Var *[depth];
        SCIP_Real *branchbounds = new SCIP_Real[depth];
        SCIP_BOUNDTYPE *boundtypes = new SCIP_BOUNDTYPE[depth];
        int n;

        SCIPnodeGetAncestorBranchings(children[0], vars, branchbounds, boundtypes, &n, depth);

        if(featureCalculator)
            featureCalculator->updateBranchCounter(children, vars[0]);

        if(history)
            history->addElement(vars[0]);

        delete[] vars;
        delete[] branchbounds;
        delete[] boundtypes;
    }

    if(featureCalculator && SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEFOCUSED) {
        // ?
    }
    return SCIP_OKAY;
}



void EventhdlrUpdateFeatures::setFeatureCalculator(FeaturesCalculator *featureCalculator) {
    EventhdlrUpdateFeatures::featureCalculator = featureCalculator;
}

void EventhdlrUpdateFeatures::setHistory(BranchingHistory *history) {
    EventhdlrUpdateFeatures::history = history;
}

void EventhdlrUpdateFeatures::informNextOneIsExploration() {
    nextOneIsExploration = true;
}
