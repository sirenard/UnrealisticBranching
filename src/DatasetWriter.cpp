//
// Created by simon on 9/09/22.
//

#include "DatasetWriter.h"

DatasetWriter::DatasetWriter(const char *nodeSelectionFilename, const char *branchFilename):
nodeSelectionStream(),
branchStream(){
    branchStream.open(branchFilename, std::fstream::out);
}

SCIP_RETCODE
DatasetWriter::addNode(SCIP *scip, SCIP_NODE *node, int nlpcands, int *varScores, SCIP_VAR **lpcands, int bestScore,
                       int bestCand) {
    // compute the features
    featuresCalculator->computeDynamicProblemFeatures(scip);
    int maxScore=0;
    for(int i=0; i<nlpcands; ++i){
        if(varScores[i] > maxScore)maxScore=varScores[i];
    }

    for(int i=0; i<nlpcands; ++i){
        if(varScores[i] == 0)continue;
        double score = 1.0 - ((double)varScores[i]/maxScore);
        SCIP_VAR* var = lpcands[i];
        writeLine(var, score);
    }

    featuresCalculator->updateBranchCounter(node, lpcands[bestCand]);
    return SCIP_OKAY;
}

DatasetWriter::~DatasetWriter() {
    branchStream.close();
    //nodeSelectionStream.close();
}

void DatasetWriter::writeLine(SCIP_VAR *var, double score) {
    int arraySizes[3] = {
        featuresCalculator->getNStaticFeatures(),
        featuresCalculator->getNDynamicFeatures(),
        featuresCalculator->getNObjectiveIncreaseStatics()
    };
    const double* features[3] = {
            featuresCalculator->getStaticFeatures(var),
            featuresCalculator->getDynamicProblemFeatures(var),
            featuresCalculator->getDynamicOptimizationFeatures(var)
    };

    std::string line;

    for(int i=0; i<3; ++i){
        //line += "start " + std::to_string(i) + ";";
        for(int k=0; k<arraySizes[i]; ++k){
            line += std::to_string(features[i][k]) + ";";
        }
    }

    line += std::to_string(score) +  "\n";

    branchStream.write(line.c_str(), line.length());
}

void DatasetWriter::setFeaturesCalculator(FeaturesCalculator *featuresCalculator) {
    DatasetWriter::featuresCalculator = featuresCalculator;
}
