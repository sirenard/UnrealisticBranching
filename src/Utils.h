//
// Created by simon on 14/08/22.
//

#ifndef LEARNING_UTILS_H
#define LEARNING_UTILS_H


#include <scip/scip.h>


class Utils {
public:
    static SCIP_Retcode create_scip_instance(SCIP** scipp, bool addBranchScheme=true);
    static SCIP_Retcode configure_scip_instance(SCIP* scip, bool addBranchScheme);
    static SCIP_Retcode configure_slave_scip_instance(SCIP* scip);
};


#endif //LEARNING_UTILS_H
