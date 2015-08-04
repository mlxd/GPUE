/*** edge.cc - GPUE: Split Operator based GPU solver for Nonlinear 
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

#include "../include/edge.h"

using namespace LatticeGraph;

Edge::Edge() : uid(++suid){
}

Edge::~Edge(){
	this->n1 = NULL;
	this->n2 = NULL;
}

Edge::Edge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2) : uid(++suid){
	this->n1 = n1;
	this->n2 = n2;
	this->direction = 0;
	this->weight = 0;
}

Edge::Edge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2, int dir, double weight) : uid(++suid){
	this->n1 = n1;
	this->n2 = n2;
}

/***
 * Returns id of Edge
 */
unsigned int Edge::getUid(){
	return uid;
}

int Edge::getDirection(){
	return this->direction;
}

/***
 * Get Node at index idx.
 */
std::shared_ptr<Node> Edge::getNode(int idx){
	if(idx>0){
		return this->n1;
	}
	else{
		return this->n2;
	}
}

/***
 * Sets Edge direction for a directed graph.
 */
void Edge::setDirection(int direction){
	this->direction = direction;
}

/***
 * Sets Edge weight between nodes.
 */
void Edge::setWeight(double weight){
	this->weight = weight;
}

/***
 * Replaces Node n1 or n2 with new Node n_new.
 */
void Edge::updateNode(int idx, std::shared_ptr<Node> n_new ){
	if(idx>0){
		this->n1 = n_new;
	}
	else{
		this->n2 = n_new;
	}
}

double Edge::getWeight(){
		return this->weight;
}