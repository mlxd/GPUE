//
// Created by Lee James O'Riordan on 23/07/15.
//

#ifndef LATTICEGRAPH_NODE_H
#define LATTICEGRAPH_NODE_H

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include "edge.h"
#include "tracker.h"

namespace LatticeGraph {
	class Edge;
    class Node {

    private:
	    static unsigned int suid;

	    Tracker::Vortex data;

	    std::vector<Edge *> edges; //all connected edges

    public:
	    unsigned int uid;

	    Node();
	    ~Node();

	    Node(Tracker::Vortex &data);

	    unsigned int getUid();

	    Tracker::Vortex &getData();

	    std::vector<Edge *> &getEdges(); //returns all connected edges

	    Edge *getEdge(int idx); //returns edge at index idx

	    void setData(Tracker::Vortex &data);

	    void addEdge(Node &n, int dir, double weight);

	    void addEdge(Edge &e);

	    void removeEdge(Node &n); //remove edge connecting this to Node n

	    void removeEdge(unsigned int uid); //remove edge at index idx

	    void removeEdge(Edge &e); //remove edge at index idx

	    void removeEdges(); //remove all edges

	    Node *getConnectedNode(Edge *e); //Return the node on the other side of the edge

	    void getConnectedNodes(unsigned int &nodes); //get all connected nodes to this
    };

}
#endif //LATTICEGRAPH_NODE_H
