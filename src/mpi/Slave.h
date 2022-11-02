//
// Created by simon on 26/10/22.
//

#ifndef LEARNING_SLAVE_H
#define LEARNING_SLAVE_H


#include "Node.h"

class Slave: public Node { // MPI slave
    void retrieveInstance();
    void returnScore(int score);
public:
    explicit Slave(unsigned rank);

    void run() override;

    void createScipInstance();
};


#endif //LEARNING_SLAVE_H
