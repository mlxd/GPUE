/*** edge.h - GPUE: Split Operator based GPU solver for Nonlinear
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

	    std::weak_ptr<Node> n1, n2; //Points to the connected nodes

	    int direction; //1 n1->n2, 0 undirected, -1 n2->n1

	    double weight; // 0.0 if not relevant, otherwise +/- reals

    public:
	    unsigned int uid;

	    Edge();
	    ~Edge();

	    Edge(std::weak_ptr<Node> n1, std::weak_ptr<Node> n2);
	    Edge(std::weak_ptr<Node> n1, std::weak_ptr<Node> n2, int dir, double weight);

	    unsigned int getUid();
	    unsigned int& getSuid();
	    int getDirection(); //Return nodes in order of flow
	    double getWeight();
	    std::weak_ptr<Node> getVortex(int idx);

	    void setDirection(int direction);
	    void setWeight(double weight);

	    void updateVortex(int idx, std::weak_ptr<Node> n_new);

	    bool isMember(std::weak_ptr<Node> n);


	    bool operator < (std::shared_ptr<Edge> e) const{
		    return uid < e->getUid();
	    }
    };

}
#endif //LATTICEGRAPH_EDGE_H
