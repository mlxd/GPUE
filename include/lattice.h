///@endcond
//##############################################################################
/**
 *  @file    lattice.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Class for treating vortices as a graph
 *
 *  @section DESCRIPTION
 *  Vortices are passed in as nodes, with edges created between them for
 *  generating graphs. Edges and nodes can be created as necessary. An adjacency
 *  matrix can also be output, with Mathematica input syntax in mind.
 */
//##############################################################################

#ifndef LATTICEGRAPH_LATTICE_H
#define LATTICEGRAPH_LATTICE_H

#include <cstdlib>
#include <cmath>
//#include <algorithm>
#include <vector>
#include "node.h"
#include "edge.h"
#include "tracker.h"

namespace LatticeGraph {

    class Lattice {

    private:
	    std::vector <std::shared_ptr <Node> > vortices; //The vortices (nodes)
	    std::vector <std::shared_ptr <Edge> > edges; //The edges

    public:
		/**
		* @brief	Makes stuff exist
		* @ingroup	graph
		*/
	    Lattice();
		/**
		* @brief	Makes stuff anti-exist
		* @ingroup	graph
		*/
	    ~Lattice();

//##############################################################################

		/**
		* @brief	Returns the vectors for vortices
		* @ingroup	graph
		* @return	Vector of shared_ptr for vortices (Nodes)
		*/
	    std::vector <std::shared_ptr <Node> > &getVortices();
		/**
		* @brief	Returns the edges for edges
		* @ingroup	graph
		* @return	Vector of shared_ptr for edges (Nodes)
		*/
	    std::vector <std::shared_ptr <Edge> > &getEdges();

//##############################################################################

		/**
		* @brief	Returns vortex specified at index idx
		* @ingroup	graph
		* @param	idx Index of vortex position
		* @return	Shared pointer to vortex (Node) at index idx
		*/
	    std::shared_ptr<Node> getVortexIdx(unsigned int idx);
		/**
		* @brief	Returns edge specified at index idx
		* @ingroup	graph
		* @param	idx Index of vortex position
		* @return	Shared pointer to edge at index idx
		*/
	    std::shared_ptr<Edge> getEdgeIdx(unsigned int idx);

//##############################################################################

		/**
		* @brief	Returns vortex index based on UID
		* @ingroup	graph
		* @param	uid UID of vortex
		* @return	Index of vortex with UID uid
		*/
	    unsigned int getVortexIdxUid(unsigned int uid);
		/**
		* @brief	Returns edge index based on UID
		* @ingroup	graph
		* @param	uid UID of edge
		* @return	Index of edge with UID uid
		*/
	    unsigned int getEdgeIdxUid(unsigned int uid);

//##############################################################################

		/**
		* @brief	Returns vortex based on UID
		* @ingroup	graph
		* @param	uid UID of vortex
		* @return	Vortex with UID uid
		* @bug		Fails for vortex UIDs that do not exist. Assumes they do.
		*/
	    std::shared_ptr<Node> getVortexUid(unsigned int uid);
		/**
		* @brief	Returns edge based on UID
		* @ingroup	graph
		* @param	uid UID of edge
		* @return	edge with UID uid
		* @bug		Fails for vortex UIDs that do not exist. Assumes they do.
		*/
	    std::shared_ptr<Edge> getEdgeUid(unsigned int uid);

//##############################################################################

		/**
		* @brief	Calculates distance between vortices
		* @ingroup	graph
		* @param	n1 Vortex (node) 1
		* @param	n2 Vortex (node) 2
		* @return	double of intervortex distance
		*/
	    double getVortexDistance(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);
		/**
		* @brief	Calculates distance between vortices
		* @ingroup	graph
		* @param	n1 Vortex (node) 1
		* @param	n2 Vortex (node) 2
		* @return	double of intervortex distance
		*/
	    double getVortexDistanceD(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);

//##############################################################################

		/**
		* @brief	Set vortex (node) at a specific position
		* @ingroup	graph
		* @param	idx Index to set vortex
		* @param	n Vortex to set
		*/
	    void setVortex(unsigned int idx, std::shared_ptr<Node> n);
		/**
		* @brief	Set edge at a specific position
		* @ingroup	graph
		* @param	idx Index to set vortex
		* @param	n Vortex to set
		*/
	    void setEdge(unsigned int idx, std::shared_ptr<Edge> e);

//##############################################################################

		/**
		* @brief	Add vortex (node) to the lattice (graph).
		* @ingroup	graph
		* @param	n Vortex to add
		*/
	    void addVortex(std::shared_ptr<Node> n);
		/**
		* @brief	Add edge e to the lattice (graph). Assumes edge has nodes already.
		* @ingroup	graph
		* @param	e Edge to add
		*/
	    void addEdge(std::shared_ptr<Edge> e);
		/**
		* @brief	Add edge between vortex (node) n1 and n2 to the lattice (graph).
		* @ingroup	graph
		* @param	n1 Vortex to add edge to as connection 1
		* @param	n2 Vortex to add edge to as connection 2
		*/
	    void addEdge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);
		/**
		* @brief	Add edge e between vortex (node) n1 and n2 to the lattice (graph).
		* @ingroup	graph
		* @param	e Edge to add
		* @param	n1 Vortex to add edge to as connection 1
		* @param	n2 Vortex to add edge to as connection 2
		*/
	    void addEdge(std::shared_ptr<Edge> e, std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);

//##############################################################################

		/**
		* @brief	Remove vortex n from the lattice
		* @ingroup	graph
		* @param	n Vortex (node) to remove
		*/
	    void removeVortex(std::shared_ptr<Node> n);
		/**
		* @brief	Remove vortex from the lattice based on index idx
		* @ingroup	graph
		* @param	idx Vortex at index idx to remove
		*/
	    void removeVortexIdx(unsigned int idx);
		/**
		* @brief	Remove vortex from the lattice based on UID uid
		* @ingroup	graph
		* @param	uid Vortex with UID uid to remove
		*/
	    void removeVortexUid(unsigned int uid);

//##############################################################################

		/**
		* @brief	Remove edge between n1 and n2 from the lattice
		* @ingroup	graph
		* @param	n1 Vortex 1
		* @param	n2 Vortex 2
		*/
	    void removeEdge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);
		/**
		* @brief	Remove edge from the lattice
		* @ingroup	graph
		* @param	e Edge to remove
		*/
	    void removeEdge(std::shared_ptr<Edge> e);
		/**
		* @brief	Remove edge from the lattice based on UID uid
		* @ingroup	graph
		* @param	uid Edge with UID uid to remove
		*/
	    void removeEdgeUid(unsigned int uid);
		/**
		* @brief	Remove edge from the lattice based index idx
		* @ingroup	graph
		* @param	idx Edge at index idx to remove
		*/
	    void removeEdgeIdx(unsigned int idx);
		/**
		* @brief	Remove all edges from vortex
		* @ingroup	graph
		* @param	n1 Vortex (node) to lose all edges
		*/
	    void removeEdges(std::shared_ptr<Node> n1);



		/**
		* @brief	Create a vortex at specified location with given winding
		* @ingroup	graph
		* @param	posx Position x to create vortex
		* @param	posy Position y to create vortex
		* @param	Winding of vortex to create
		*/
	    void createVortex(double posx, double posy, int winding);
		/**
		* @brief	Destroy a vortex with UID uid
		* @ingroup	graph
		* @param	uid Vortex to destroy with UID uid
		*/
	    void destroyVortex(unsigned int uid);



		/**
		* @brief	Generate edges between vortices closer than int radius
		* @ingroup	graph
		* @param	radius Radius cutoff for creating edge connections
		*/
	    void createEdges(unsigned int radius);
		/**
		* @brief	Generate edges between vortices closer than double radius
		* @ingroup	graph
		* @param	radius Radius cutoff for creating edge connections
		*/
	    void createEdges(double radius);

//##############################################################################

		/**
		* @brief	Output adjacency matrix.
		* @ingroup	graph
		* @param	*mat Modifyable UInt array for adjacency matrix
		*/
	    void genAdjMat(unsigned int *mat);  //generate adjacency matrix
		/**
		* @brief	Output adjacency matrix.
		* @ingroup	graph
		* @param	*mat Modifyable double array for adjacency matrix
		*/
	    void genAdjMat(double *mat);  //generate adjacency matrix
		/**
		* @brief	Format adjacency matrix with Mathematica-friendly output
		* @ingroup	graph
		* @param	*mat Adjacency matrix from genAdjMat(*uint)
		*/
	    void adjMatMtca(unsigned int *mat);
		/**
		* @brief	Format adjacency matrix with Mathematica-friendly output
		* @ingroup	graph
		* @param	*mat Adjacency matrix from genAdjMat(*double)
		*/
	    void adjMatMtca(double *mat);

//##############################################################################

		/**
		* @brief	Swap the elements at indices based on UID uid1 and uid2
		* @ingroup	graph
		* @param	uid1 UID of element 1 to swap
		* @param	uid2 UID of element 2 to swap
		*/
	    void swapIdxUid(unsigned int uid1, unsigned int uid2);
		/**
		* @brief	Swap the elements at indices idx1 and idx2
		* @ingroup	graph
		* @param	idx1 Index of element 1
		* @param	idx2 Index of element 2
		*/
	    void swapIdx(unsigned int idx1, unsigned int idx2);

//##############################################################################

		/**
		* @brief	Checks if vortex n1 and n2 are connected
		* @ingroup	graph
		* @param	n1 Vortex (node) 1
		* @param	n2 Vortex (node) 2
		* @return	Weak pointer to the connecting edge
		*/
	    std::weak_ptr<Edge> isConnected(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);
    };
}
#endif //LATTICEGRAPH_LATTICE_H
