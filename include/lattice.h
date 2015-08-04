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
	    //std::vector <Node> vortices; //The nodes
	    std::vector <std::shared_ptr <Node> > vortices; //The nodes

    public:
	    Lattice();
	    ~Lattice();

	    std::shared_ptr<Node> getElement(int idx);

	    int getElementUidIdx(unsigned int uid); //get element index by uid

	    std::shared_ptr<Node> getElementUid(unsigned int uid); // gets element by uid. Assumes it exists.

	    std::vector <std::shared_ptr <Node> > &getVortices(); //get vortices vector

	    double getNodeDistance(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2); //compare distances between nodes;

	    void setElement(int idx, std::shared_ptr<Node> n);

	    void addNode(std::shared_ptr<Node> n);

	    void removeNode(std::shared_ptr<Node> n);

	    void removeEdges(std::shared_ptr<Node> n1);

	    void removeEdge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);

	    void removeNode(int idx);

	    void createEdges(unsigned int radius); //Check all nodes and see if an edge can be created between pairs

	    void genAdjMat(unsigned int *mat);  //generate adjacency matrix

	    void addEdge(std::shared_ptr<Edge> e, std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);

	    void adjMatMtca(unsigned int *mat);
    };
}
#endif //LATTICEGRAPH_LATTICE_H
