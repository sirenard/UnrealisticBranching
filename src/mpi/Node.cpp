//
// Created by simon on 26/10/22.
//

#include "Node.h"

Node* Node::instance = nullptr;
Node::Node(unsigned int rank):
    rank(rank) {
}

void Node::setInstance(Node *node) {
    instance = node;
}

Node *Node::getInstance() {
    return instance;
}

const unsigned int Node::getRank() const {
    return rank;
}

bool Node::isMaster() {
    return !rank;
}
