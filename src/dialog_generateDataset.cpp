//
// Created by simon on 23/09/22.
//

#include "dialog_generateDataset.h"
#include "DatasetWriter.h"
#include "branch_unrealistic.h"
#include "Utils.h"
#include "EventhdlrUpdateFeatures.h"

#define TEST_END_DIALOG(ptr) if(ptr[0]=='\0'){*nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);return SCIP_OKAY;}

DialogGenerateDataset::DialogGenerateDataset(SCIP *scip)
        : ObjDialog(scip, DIALOG_NAME, DIALOG_DESC, DIALOG_ISSUBMENU) {}

SCIP_DECL_DIALOGEXEC(DialogGenerateDataset::scip_exec){
    char* inputfilename;
    char* outputfilename;
    char* tmp;
    int signB, signC, signA;
    SCIP_Bool endoffile;

    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter input problem file: ", &inputfilename, &endoffile) );
    if( endoffile){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(inputfilename);

    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter output filename: ", &outputfilename, &endoffile) );
    if( endoffile ){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(outputfilename);



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
    //SCIPfreeBufferArray(scipmain, &tmp);

    SCIP_CALL(SCIPreadProb(
            scip,
            inputfilename,
            "lp"
    )
    );

    DatasetWriter writer("node.csv", outputfilename);
    EventhdlrUpdateFeatures* eventHdlr = dynamic_cast<EventhdlrUpdateFeatures *>(SCIPfindObjEventhdlr(scip, EVENT_HDLR_UPDATE_FEATURES_NAME));
    FeaturesCalculator featuresCalculator(scip, signB, signC, signA);
    writer.setFeaturesCalculator(&featuresCalculator);

    Branch_unrealistic::setDataWriter(&writer);

    BranchingHistory* history = new BranchingHistory();
    Branch_unrealistic *objbranchrule = (Branch_unrealistic*)SCIPfindObjBranchrule(scip, "unrealistic");
    eventHdlr->setHistory(history);

    objbranchrule->setBranchingHistory(history);
    objbranchrule->disableCopyCatBranching();




    SCIP_CALL( SCIPsolve(scip) );

    //SCIP_Longint score = SCIPgetNNodes(scipmain);
    //SCIPinfoMessage(scipmain, NULL, ("Solved in " + std::to_string(score) + " nodes\n").c_str());

    //SCIPfreeBufferArray(scipmain, &inputfilename);
    //SCIPfreeBufferArray(scipmain, &outputfilename);

    writer.setFeaturesCalculator(nullptr);


    *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
    return SCIP_OKAY;
}

