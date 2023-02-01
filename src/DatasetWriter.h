//
// Created by simon on 9/09/22.
//

#ifndef LEARNING_DATASETWRITER_H
#define LEARNING_DATASETWRITER_H

#include <string>
#include <fstream>
#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#include "objscip/objscip.h"
#include "FeaturesCalculator.h"

class DatasetWriter {
    std::fstream nodeSelectionStream, branchStream;
    FeaturesCalculator* featuresCalculator;
public:
    DatasetWriter(const char* nodeSelectionFilename, const char* branchFilename);
    ~DatasetWriter();
    SCIP_RETCODE addNode(SCIP *scip, SCIP_NODE **node, int nlpcands, int *varScores, SCIP_VAR **lpcands, int bestCand);

    void writeLine(SCIP_VAR *var, double score, SCIP *scip);

    void setFeaturesCalculator(FeaturesCalculator *featuresCalculator);
};


#endif //LEARNING_DATASETWRITER_H
