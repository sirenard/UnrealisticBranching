//
// Created by simon on 23/09/22.
//

#ifndef LEARNING_DIALOG_GENERATEDATASET_H
#define LEARNING_DIALOG_GENERATEDATASET_H

#define DIALOG_NAME "generateDataset"
#define DIALOG_DESC "Generate a dataset from a given problem"
#define DIALOG_ISSUBMENU 0

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#include "objscip/objscip.h"

class DialogGenerateDataset: public scip::ObjDialog {
public:
    DialogGenerateDataset(SCIP *scip);
    SCIP_DECL_DIALOGEXEC(scip_exec) override;
};


#endif //LEARNING_DIALOG_GENERATEDATASET_H
