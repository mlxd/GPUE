//
// Created by Lee James O'Riordan on 23/07/15.
//

#ifndef LATTICEGRAPH_NODE_H
#define LATTICEGRAPH_NODE_H

#include <cstdlib>
#include <cmath>
#include <memory>
#include <vector>
#include "edge.h"
#include "tracker.h"

namespace LatticeGraph {
	class Edge;
    class Node {

    private:
	    static unsigned int suid;

	    Tracker::Vortex data;

	    std::vector<std::shared_ptr <Edge> > edges; //all connected edges

    public:
	    unsigned int uid;

	    Node();
	    ~Node();

	    Node(Tracker::Vortex &data);

	    unsigned int getUid();

	    Tracker::Vortex &getData();

	    std::vector<std::shared_ptr <Edge> > &getEdges(); //returns all connected edges

	    std::shared_ptr<Edge> getEdge(int idx); //returns edge at index idx

	    void setData(Tracker::Vortex &data);

	    void addEdge(std::shared_ptr<Node> n, int dir, double weight);

	    void addEdge(std::shared_ptr<Edge> e);

	    void removeEdge(std::shared_ptr<Node> n); //remove edge connecting this to Node n

	    void removeEdge(unsigned int uid); //remove edge at index idx

	    void removeEdge(std::shared_ptr<Edge> e); //remove edge at index idx

	    void removeEdges(); //remove all edges

	    std::shared_ptr<Node> getConnectedNode(std::shared_ptr<Edge> e); //Return the node on the other side of the edge

	    void getConnectedNodes(unsigned int &nodes); //get all connected nodes to this
    };

}
#endif //LATTICEGRAPH_NODE_H
