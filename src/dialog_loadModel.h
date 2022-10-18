//
// Created by simon on 18/10/22.
//

#ifndef LEARNING_DIALOG_LOADMODEL_H
#define LEARNING_DIALOG_LOADMODEL_H

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#include "objscip/objscip.h"

class DialogLoadModel: public scip::ObjDialog {
public:
    DialogLoadModel(SCIP *scip);
    SCIP_DECL_DIALOGEXEC(scip_exec) override;
};


#endif //LEARNING_DIALOG_LOADMODEL_H
