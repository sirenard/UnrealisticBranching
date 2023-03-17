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
DatasetWriter::addNode(SCIP *scip, int nlpcands, int *varScores, SCIP_VAR **lpcands, char scoreMethod, double alpha) {
    // compute the features
    int maxScore=0;
    int minScore=INT_MAX;
    for(int i=0; i<nlpcands; ++i){
        if(varScores[i] == INT_MAX)continue;
        if(varScores[i] > maxScore)maxScore=varScores[i];
        if(varScores[i] < minScore)minScore=varScores[i];
    }

    double* nodeFeatures = featuresCalculator->getTreeFeatures(scip, nlpcands);
    int nodeFeaturesSize = featuresCalculator->getTreeFeaturesSize();

    for(int i=0; i<nlpcands; ++i){
        if(varScores[i] == INT_MAX)continue;
        if(maxScore == minScore)continue;

        double score;
        switch (scoreMethod) {
            case('c'):
                //score = (double)(maxScore-varScores[i])/(maxScore-minScore);
                score = 1.0 - (double)varScores[i]/maxScore;
                break;
            case('a'):
                score = varScores[i] <= (1.0+alpha)*minScore;
                break;
            default:
                score=-1;
        }

        SCIP_VAR* var = lpcands[i];
        writeLine(var, score, scip, varScores[i], minScore, maxScore, nodeFeatures, nodeFeaturesSize);
    }

    delete[] nodeFeatures;
    return SCIP_OKAY;
}

DatasetWriter::~DatasetWriter() {
    branchStream.close();
    //nodeSelectionStream.close();
}

void
DatasetWriter::writeLine(SCIP_VAR *var, double score, SCIP *scip, int ubScore, int smallestUbScore, int biggestUbScore,
                         double *nodeFeatures, int nodeFeaturesSize) {
    // TODO: use getFeatures instead
    double* varFeatures = featuresCalculator->getNoTreeFeatures(scip, var);
    int varFeaturesSize = featuresCalculator->getNoTreeFeaturesSize();


    std::string line;

    for(int k=0; k<nodeFeaturesSize; ++k){
        line += std::to_string(nodeFeatures[k]) + ",";
    }
    for(int k=0; k<varFeaturesSize; ++k){
        line += std::to_string(varFeatures[k]) + ",";
    }



    if(score>=0){
        line += std::to_string(score)+"\n";
    } else {
        line += std::to_string(ubScore) + "," + std::to_string(smallestUbScore) + "," +
                std::to_string(biggestUbScore) + "\n";
    }

    branchStream.write(line.c_str(), line.length());
    delete[] varFeatures;
}

void DatasetWriter::setFeaturesCalculator(FeaturesCalculator *featuresCalculator) {
    DatasetWriter::featuresCalculator = featuresCalculator;
}
