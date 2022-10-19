//
// Created by simon on 19/10/22.
//

#include "dialog_trainmodel.h"
#include "RegressionModel.h"

#define DIALOG_NAME "trainmodel"
#define DIALOG_DESC "Train a random forest"
#define DIALOG_ISSUBMENU 0

#define TEST_END_DIALOG(ptr) if(ptr[0]=='\0'){*nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);return SCIP_OKAY;}

DialogTrainModel::DialogTrainModel(SCIP *scip)
        : ObjDialog(scip, DIALOG_NAME, DIALOG_DESC, DIALOG_ISSUBMENU) {}

SCIP_DECL_DIALOGEXEC(DialogTrainModel::scip_exec) {
    char* modelPath;
    char* datasetPath;
    char* tmp;
    int nTrees;
    SCIP_Bool endoffile;


    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Enter the path to store the model: ", &modelPath, &endoffile) );
    if( endoffile){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(modelPath);



    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Enter the path of the dataset: ", &datasetPath, &endoffile) );
    if( endoffile ){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(datasetPath);

    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Enter the number of trees for the random forest (if not a integer, 100 by default) : ", &tmp, &endoffile) );
    if( endoffile ){
        *nextdialog = NULL;
        return SCIP_OKAY;
    }
    TEST_END_DIALOG(tmp);

    char* pEnd;
    nTrees = strtol (tmp,&pEnd,10);
    if(nTrees<=0)nTrees = 100;

    RegressionModel model(nTrees);
    try{
        double mse = model.train(datasetPath);
        SCIPinfoMessage(scip, NULL, ("MSE: " + std::to_string(mse) + "\n").c_str());
    } catch(std::invalid_argument &e) {
        SCIPerrorMessage((std::string(e.what()) + "\n").c_str());
    }


    std::fstream stream;
    stream.open(modelPath, std::fstream::out);
    if(!stream.is_open()){
        SCIPinfoMessage(scip, NULL, ("Cannot open file " + std::string(modelPath) +"\n").c_str());
    } else {
        model.save(stream);
    }
    stream.close();


    *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
    return SCIP_OKAY;
}
