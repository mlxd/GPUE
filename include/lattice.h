/*** lattice.h - GPUE: Split Operator based GPU solver for Nonlinear
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
	    //std::vector <Node> vortices; //The nodes
	    std::vector <std::shared_ptr <Node> > vortices; //The nodes
	    std::vector <std::shared_ptr <Edge> > edges; //The nodes

    public:
	    Lattice();
	    ~Lattice();

	    /**
	     * Returns the vectors for vortices or edges
	     */
	    std::vector <std::shared_ptr <Node> > &getVortices(); //get vortices vector
	    std::vector <std::shared_ptr <Edge> > &getEdges(); //get edges vector

	    /**
	     * Returns vortex/edge specified at index idx
	     */
	    std::shared_ptr<Node> getVortexIdx(unsigned int idx);
	    std::shared_ptr<Edge> getEdgeIdx(unsigned int idx);

	    /**
	     * Returns the index of vortex/edge with UID uid
	     */
	    unsigned int getVortexIdxUid(unsigned int uid); //get vortex index by uid
	    unsigned int getEdgeIdxUid(unsigned int uid); //get edge index by uid

	    /**
	     * Returns vortex/edge with UID uid
	     */
	    std::shared_ptr<Node> getVortexUid(unsigned int uid); // gets vortex by uid. Assumes it exists.
	    std::shared_ptr<Edge> getEdgeUid(unsigned int uid); // gets edge by uid. Assumes it exists.

	    /**
	     * Calculates distance between vortices
	     */
	    double getVortexDistance(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2); //compare distances between nodes;
	    double getVortexDistanceD(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2); //compare distances between nodes;

	    /**
	     * Set vortex/edge at a specific position
	     */
	    void setVortex(unsigned int idx, std::shared_ptr<Node> n);
	    void setEdge(unsigned int idx, std::shared_ptr<Edge> e);

	    /**
	     * Add vortex/edge to the lattice.
	     */
	    void addVortex(std::shared_ptr<Node> n);
	    void addEdge(std::shared_ptr<Edge> e);
	    void addEdge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);
	    void addEdge(std::shared_ptr<Edge> e, std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);

	    /**
	     * Remove vortex/edge from the lattice
	     */
	    void removeVortex(std::shared_ptr<Node> n);
	    void removeVortexIdx(unsigned int idx);
	    void removeVortexUid(unsigned int idx);
	    void removeEdge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);
	    void removeEdge(std::shared_ptr<Edge> e);
	    void removeEdgeUid(unsigned int uid);
	    void removeEdgeIdx(unsigned int idx);
	    void removeEdges(std::shared_ptr<Node> n1);

	    /**
	     * Modify wavefunction
	     */
	    void createVortex(double posx, double posy, int winding);
	    void destroyVortex(unsigned int uid);

	    /**
	     * Generate edges between vortices closer than radius
	     */
	    void createEdges(unsigned int radius); //Check all nodes and see if an edge can be created between pairs
	    void createEdges(double radius); //Check all nodes and see if an edge can be created between pairs

	    /**
	     * Generate adjacency matrix. Format with Mathematic-friendly output
	     */
	    void genAdjMat(unsigned int *mat);  //generate adjacency matrix
	    void genAdjMat(double *mat);  //generate adjacency matrix
	    void adjMatMtca(unsigned int *mat);
	    void adjMatMtca(double *mat);

	    void swapIdxUid(unsigned int uid1, unsigned int uid2);
	    void swapIdx(unsigned int idx1, unsigned int idx2);

	    std::weak_ptr<Edge> isConnected(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2);
    };
}
#endif //LATTICEGRAPH_LATTICE_H
