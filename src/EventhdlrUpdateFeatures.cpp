//
// Created by simon on 17/03/23.
//

#include <iostream>
#include "EventhdlrUpdateFeatures.h"

EventhdlrUpdateFeatures::EventhdlrUpdateFeatures(SCIP *scip):
ObjEventhdlr(scip, EVENT_HDLR_UPDATE_FEATURES_NAME,"event handler for new focused node, update if it exists, the data for features computation"),
featureCalculator(nullptr)
{
}

SCIP_DECL_EVENTINITSOL(EventhdlrUpdateFeatures::scip_initsol){
    SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL) );
    return SCIP_OKAY;
}


SCIP_DECL_EVENTEXEC(EventhdlrUpdateFeatures::scip_exec){
    if(featureCalculator){
        featureCalculator->updateBranching(scip);
    }
    return SCIP_OKAY;
}



void EventhdlrUpdateFeatures::setFeatureCalculator(FeaturesCalculator *featureCalculator) {
    EventhdlrUpdateFeatures::featureCalculator = featureCalculator;
}
