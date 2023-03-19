//
// Created by simon on 17/03/23.
//

#ifndef LEARNING_EVENTHDLRUPDATEFEATURES_H
#define LEARNING_EVENTHDLRUPDATEFEATURES_H

#define EVENT_HDLR_UPDATE_FEATURES_NAME "hdlr_updatefeatures"
#include "objscip/objscip.h"
#include "FeaturesCalculator.h"
#include "BranchingHistory.h"


class EventhdlrUpdateFeatures: public scip::ObjEventhdlr {
    FeaturesCalculator* featureCalculator;
    BranchingHistory* history;
    bool nextOneIsExploration;
    public:
    /** default constructor */
    EventhdlrUpdateFeatures(
            SCIP* scip
    );

    /** destructor */
    virtual ~EventhdlrUpdateFeatures()
    {
    }

    virtual SCIP_DECL_EVENTINITSOL(scip_initsol);

    virtual SCIP_DECL_EVENTEXEC(scip_exec);

    void setFeatureCalculator(FeaturesCalculator *featureCalculator);

    void setHistory(BranchingHistory *history);

    void informNextOneIsExploration();
};


#endif //LEARNING_EVENTHDLRUPDATEFEATURES_H
