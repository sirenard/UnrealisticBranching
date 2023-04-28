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
    featuresCalculator->computeDynamicProblemFeatures(scip, lpcands, nlpcands);
    for(int i=0; i<nlpcands; ++i){
        if(varScores[i] == INT_MAX)continue;
        if(varScores[i] > maxScore)maxScore=varScores[i];
        if(varScores[i] < minScore)minScore=varScores[i];
    }

    if(scoreMethod == '1'){
        for(int i=0; i<nlpcands; ++i) {
            if (maxScore == minScore)continue;
            for (int j = 0; j < nlpcands; ++j) {
                if (i == j || varScores[i] == INT_MAX && varScores[j] == INT_MAX)continue;
                double score = varScores[i] < varScores[j];

                SCIP_Var* vars[] = {lpcands[i], lpcands[j]};
                writeLine(vars, 2, score, scip, varScores[i], minScore, maxScore);
            }
        }
    } else for(int i=0; i<nlpcands; ++i){
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
        writeLine(&var, 1, score, scip, varScores[i], minScore, maxScore);
    }

    return SCIP_OKAY;
}

DatasetWriter::~DatasetWriter() {
    branchStream.close();
    //nodeSelectionStream.close();
}

void
DatasetWriter::writeLine(SCIP_VAR **vars, int nVars, double score, SCIP *scip, int ubScore, int smallestUbScore,
                         int biggestUbScore) {
    // TODO: use getFeatures instead
    std::string line;

    int arraySizes[3] = {
            featuresCalculator->getNStaticFeatures(),
            featuresCalculator->getNDynamicFeatures(),
            featuresCalculator->getNObjectiveIncreaseStatics()
    };

    for(int s=0;s<nVars;s++) {
        SCIP_Var *var = vars[s];
        const double *features[3] = {
                featuresCalculator->getStaticFeatures(var, scip),
                featuresCalculator->getDynamicProblemFeatures(var),
                featuresCalculator->getDynamicOptimizationFeatures(var)
        };


        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < arraySizes[i]; ++k) {
                line += std::to_string(features[i][k]) + ",";
            }
        }
    }

    if(score>=0){
        line += std::to_string(score)+"\n";
    } else {
        line += std::to_string(ubScore) + "," + std::to_string(smallestUbScore) + "," +
                std::to_string(biggestUbScore) + "\n";
    }

    branchStream.write(line.c_str(), line.length());
}

void DatasetWriter::setFeaturesCalculator(FeaturesCalculator *featuresCalculator) {
    DatasetWriter::featuresCalculator = featuresCalculator;
}
