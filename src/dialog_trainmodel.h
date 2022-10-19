//
// Created by simon on 19/10/22.
//

#ifndef LEARNING_DIALOG_TRAINMODEL_H
#define LEARNING_DIALOG_TRAINMODEL_H

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#include "objscip/objscip.h"

class DialogTrainModel: public scip::ObjDialog {
public:
    DialogTrainModel(SCIP *scip);
    SCIP_DECL_DIALOGEXEC(scip_exec) override;
};


#endif //LEARNING_DIALOG_TRAINMODEL_H
