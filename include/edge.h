//
// Created by Lee James O'Riordan on 23/07/15.
//

#ifndef LATTICEGRAPH_EDGE_H
#define LATTICEGRAPH_EDGE_H

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include "node.h"

namespace LatticeGraph {

    class Node;
    class Edge {

    private:
	    static unsigned int suid; //Incremented id for new

	    Node *n1, *n2; //Points to the connected nodes

	    int direction; //1 n1->n2, 0 undirected, -1 n2->n1

	    double weight; // 0.0 if not relevant, otherwise +/- reals

    public:
	    unsigned int uid;

	    Edge();
	    ~Edge();

	    Edge(Node &n1, Node &n2);
	    Edge(Node &n1, Node &n2, int dir, double weight);

	    unsigned int getUid();

	    int getDirection(); //Return nodes in order of flow

	    double getWeight();

	    Node *getNode(int idx);

	    void setDirection(int direction);

	    void setWeight(double weight);

	    void updateNode(int idx, Node *n_new);
    };

}
#endif //LATTICEGRAPH_EDGE_H
