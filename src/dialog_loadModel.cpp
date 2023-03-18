//
// Created by simon on 18/10/22.
//

#include "dialog_loadModel.h"
#include "FeaturesCalculator.h"
#include "RegressionModel.h"
#include "branch_unrealisticTrained.h"
#include "EventhdlrUpdateFeatures.h"

#define DIALOG_NAME "loadmodel"
#define DIALOG_DESC "Load a regression model that the unrealistic trained branching will use"
#define DIALOG_ISSUBMENU 0

#define TEST_END_DIALOG(ptr) if(ptr[0]=='\0'){*nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);return SCIP_OKAY;}

DialogLoadModel::DialogLoadModel(SCIP *scip)
        : ObjDialog(scip, DIALOG_NAME, DIALOG_DESC, DIALOG_ISSUBMENU) {}

SCIP_DECL_DIALOGEXEC(DialogLoadModel::scip_exec) {
    char* modelPath;
    char* tmp;
    int signB, signC, signA;
    SCIP_Bool endoffile;


    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter the path to the model: ", &modelPath, &endoffile) );
    if( endoffile){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(modelPath);



    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter input problems signA: ", &tmp, &endoffile) );
    if( endoffile ){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(tmp);
    signA = atoi(tmp);
    //SCIPfreeBufferArray(scipmain, &tmp);

    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter input problems signB: ", &tmp, &endoffile) );
    if( endoffile ){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(tmp);
    signB = atoi(tmp);
    //SCIPfreeBufferArray(scipmain, &tmp);


    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter input problems signC: ", &tmp, &endoffile) );
    if( endoffile ){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(tmp);
    signC = atoi(tmp);

    auto *featuresCalculator = new FeaturesCalculator(scip, signB, signC, signA);
    auto *model = new RegressionModel(modelPath);

    Branch_unrealisticTrained::setFeaturesCalculator(featuresCalculator);
    Branch_unrealisticTrained::setModel(model);

    EventhdlrUpdateFeatures* eventHdlr = dynamic_cast<EventhdlrUpdateFeatures *>(SCIPfindObjEventhdlr(scip, EVENT_HDLR_UPDATE_FEATURES_NAME));
    eventHdlr->setFeatureCalculator(featuresCalculator);

    *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
    return SCIP_OKAY;
}
