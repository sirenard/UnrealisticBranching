//
// Created by simon on 18/10/22.
//

#include "RegressionModel.h"

RegressionModel::RegressionModel(std::string path) {
    std::fstream f(path);
    deserialize(rf, f);
}

void RegressionModel::train(std::string csvPath) {
    std::vector<sample_type> x;
    std::vector <double> y;

    readDataSet(x, y, csvPath);
    dlib::random_forest_regression_trainer<dlib::dense_feature_extractor> trainer;
    trainer.be_verbose();
    rf = trainer.train(x, y);
    //std::cout << x(9, 8) << std::endl;

}

void RegressionModel::readDataSet(std::vector<sample_type> &x, std::vector<double> &y, std::string path) {
    std::ifstream data(path);
    std::string line;
    std::string cell;

    bool firstLine = true;
    int nCol=0;
    std::getline(data, line);
    std::stringstream lineStream(line);
    while (std::getline(lineStream, cell, ';')) {
        ++nCol;
    }

    --nCol;

    data.close();
    data.open(path);

    //int nCol = 22;


    int lineNumber = 0;
    while (std::getline(data, line)) {
        int columnNumber = 0;
        std::stringstream lineStream(line);
        std::string cell;
        std::vector <double> parsedRow;

        double val;
        bool first = false;
        sample_type features;
        features.set_size(nCol);
        while (std::getline(lineStream, cell, ';')) {
            if(first){
                features(columnNumber++) = val;
            }
            first = true;
            val = std::stod(cell);
        }
        lineNumber++;
        y.push_back(val);
        x.push_back(features);
    }
    data.close();
}

void RegressionModel::save(std::fstream out) {
    serialize(rf, out);
}

RegressionModel::RegressionModel() {

}
