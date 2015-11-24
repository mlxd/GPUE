///@cond LICENSE
/*** node.h - GPUE: Split Operator based GPU solver for Nonlinear
Schrodinger Equation, Copyright (C) 2011-2015, Lee J. O'Riordan
<loriordan@gmail.com>, Tadhg Morgan, Neil Crowley.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
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
