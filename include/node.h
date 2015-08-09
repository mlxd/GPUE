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
	    std::vector<std::weak_ptr <Edge> > edges; //all connected edges

    public:
	    unsigned int uid;

	    Node();
	    ~Node();

	    Node(Tracker::Vortex &data);

	    unsigned int getUid();
	    Tracker::Vortex &getData();
	    std::vector<std::weak_ptr <Edge> > &getEdges(); //returns all connected edges
	    std::weak_ptr<Edge> getEdge(int idx); //returns edge at index idx

	    std::shared_ptr<Node> getConnectedNode(std::shared_ptr<Edge> e); //Return the node on the other side of the edge
	    void getConnectedNodes(unsigned int &nodes); //get all connected nodes to this

	    void setData(Tracker::Vortex &data);

	    //void addEdge(std::shared_ptr<Node> n, int dir, double weight);
	    void addEdge(std::weak_ptr<Edge> e);

	    void removeEdge(std::shared_ptr<Node> n); //remove edge connecting this to Node n
	    void removeEdge(unsigned int uid); //remove edge with UID uid
	    void removeEdge(std::weak_ptr<Edge> e); //remove edge at index idx
	    void removeEdges(); //remove all edges
    };
}
#endif //LATTICEGRAPH_NODE_H
