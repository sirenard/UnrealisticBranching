//
// Created by simon on 18/10/22.
//

#ifndef LEARNING_REGRESSIONMODEL_H
#define LEARNING_REGRESSIONMODEL_H

#include <dlib/random_forest.h>
typedef dlib::matrix<double,0 ,1> sample_type;

class RegressionModel {
    dlib::random_forest_regression_function<> rf;
    void readDataSet( std::vector<sample_type> &x, std::vector<double> &y, std::string path);
public:
    RegressionModel();
    RegressionModel(std::string path);
    void train(std::string csvPath);
    void save(std::fstream out);
    double predictScore(const std::vector<double> &features);
};


#endif //LEARNING_REGRESSIONMODEL_H
