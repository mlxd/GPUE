///@endcond
//##############################################################################
/**
 *  @file    node.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Allow vortex to be treated as node in a graph
 *
 *  @section DESCRIPTION
 *  As part of the lattice graph generation, this class allows a vortex to be
 *	treated as a node. Edges can be created or destroyed, connected nodes can
 *	be returned, UID can be determine.
 */
 //##############################################################################

#ifndef LATTICEGRAPH_NODE_H
#define LATTICEGRAPH_NODE_H

#include <cstdlib>
#include <cmath>
#include <memory>
#include <vector>
#include "edge.h"
#include "tracker.h"

//##############################################################################

namespace LatticeGraph {
	class Edge;
    class Node {

    private:
	    static unsigned int suid;
	    Vtx::Vortex data;
	    std::vector<std::weak_ptr <Edge> > edges; //all connected edges

    public:
	    unsigned int uid;

	    Node();
	    ~Node();

	    Node(Vtx::Vortex &data);

//##############################################################################

		/**
		* @brief	Get vortex (node) UID
		* @ingroup	graph
		* @return	Vortex UID
		*/
	    unsigned int getUid();
		/**
		* @brief	Get vortex (node) static UID for new UID generation
		* @ingroup	graph
		* @return	Static class UID
		*/
	    unsigned int& getSuid();
		/**
		* @brief	Get vortex (node) data struct
		* @ingroup	graph
		* @return	Vortex data struct
		*/
	    Vtx::Vortex &getData();
		/**
		* @brief	Get all connected edges to vortex (node)
		* @ingroup	graph
		* @return	Vector of weak pointers to the connected edges
		*/
	    std::vector<std::weak_ptr <Edge> > &getEdges();
		/**
		* @brief	Get edge at index idx. Assumes indices exist.
		* @ingroup	graph
		* @param	idx Index of the requested edge
		* @return	Vector of weak pointers to the requested edge at index idx
		*/
	    std::weak_ptr<Edge> getEdge(int idx);
		/**
		* @brief	Get the node on the other side of the edge e
		* @ingroup	graph
		* @param	e Edge sharing connection with required node
		* @return	Share pointer to connected vortex (node)
		*/
	    std::shared_ptr<Node> getConnectedNode(std::shared_ptr<Edge> e);
		/**
		* @brief	Get all connected nodes to the current vortex. PassByRef.
		* @ingroup	graph
		* @param	&nodes Pass by reference location for nodes result
		*/
	    void getConnectedNodes(unsigned int &nodes);

//##############################################################################

		/**
		* @brief	Set the vortex data (in node)
		* @ingroup	graph
		* @param	&data Reference to vortex struct to set.
		*/
	    void setData(Vtx::Vortex &data);

//##############################################################################
	    //void addEdge(std::shared_ptr<Node> n, int dir, double weight);

		/**
		* @brief	Add edge e to the current vortex (node)
		* @ingroup	graph
		* @param	e Weak pointer to the edge to add
		*/
	    void addEdge(std::weak_ptr<Edge> e);

//##############################################################################

		/**
		* @brief	Remove edge connecting this to Node n
		* @ingroup	graph
		* @param	n Shared pointer to vortex (node) edge connected with.
		*/
	    void removeEdge(std::shared_ptr<Node> n);
		/**
		* @brief	Remove edge with UID uid
		* @ingroup	graph
		* @param	uid UID of requested edge to remove
		*/
	    void removeEdgeUid(unsigned int uid);
		/**
		* @brief	Remove edge with index idx
		* @ingroup	graph
		* @param	idx Index of the requested edge to remove
		*/
	    void removeEdgeIdx(unsigned int idx);
		/**
		* @brief	Remove edge e directly
		* @ingroup	graph
		* @param	e Shared pointer to edge for removal
		*/
	    void removeEdge(std::weak_ptr<Edge> e);
		/**
		* @brief	Remove all connected edges
		* @ingroup	graph
		*/
	    void removeEdges();

//##############################################################################
    };
}
#endif //LATTICEGRAPH_NODE_H
