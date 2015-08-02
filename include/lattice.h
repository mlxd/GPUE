//
// Created by Lee James O'Riordan on 23/07/15.
//

#ifndef LATTICEGRAPH_LATTICE_H
#define LATTICEGRAPH_LATTICE_H

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include "node.h"
#include "edge.h"
#include "tracker.h"

namespace LatticeGraph {

    class Lattice {

    private:
	    std::vector <Node> vortices; //The nodes

    public:
	    Lattice();
	    ~Lattice();

	    Node &getElement(int idx);

	    int getElementUidIdx(unsigned int uid); //get element index by uid

	    Node &getElementUid(unsigned int uid); // gets element by uid. Assumes it exists.

	    std::vector <Node> &getVortices(); //get vortices vector

	    double getNodeDistance(Node &n1, Node &n2); //compare distances between nodes;

	    void setElement(int idx, Node &vtx);

	    void addNode(Node& vtx);

	    void removeNode(Node &n);

	    void removeEdges(Node &n1);

	    void removeEdge(Node &n1, Node &n2);

	    void removeNode(int idx);

	    void createEdges(); //Check all nodes and see if an edge can be created between pairs

	    void genAdjMat(unsigned int *mat);  //generate adjacency matrix

	    void addEdge(Edge &e, Node &n1, Node &n2);

	    void adjMatMtca(unsigned int *mat);
    };
}
#endif //LATTICEGRAPH_LATTICE_H
