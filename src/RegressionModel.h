//
// Created by simon on 18/10/22.
//

#ifndef LEARNING_REGRESSIONMODEL_H
#define LEARNING_REGRESSIONMODEL_H

#include <dlib/random_forest.h>
typedef dlib::matrix<double,0 ,1> sample_type;

class RegressionModel {
    dlib::random_forest_regression_function<> rf;
    int nTrees;
    static void readDataSet( std::vector<sample_type> &x, std::vector<double> &y, std::string path);
public:
    RegressionModel(int nTrees=100);
    RegressionModel(std::string path);

    /**
     * Train a random forest and output the mean square error of the training using a leave one out approach.
     * @param csvPath path to the dataset
     * @return mse
     */
    double train(std::string csvPath);
    void save(std::fstream &out);
    void setNTrees(int newNTrees);
    double predictScore(const std::vector<double> &features);
};


#endif //LEARNING_REGRESSIONMODEL_H
