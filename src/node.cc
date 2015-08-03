/*** node.cc - GPUE: Split Operator based GPU solver for Nonlinear 
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

#include "../include/node.h"

using namespace LatticeGraph;

Node::Node():uid(++suid){
}

Node::~Node(){
	//Need to delete all edges connecting this node + others
	//this->Node::removeEdges(); //May be necessary, but for-each loop fail with it
	//delete this;
}

Node::Node(Tracker::Vortex& data):uid(++suid){
	this->data = data;
}

unsigned int Node::getUid(){
	return uid;
}

Tracker::Vortex& Node::getData(){
	return this->data;
}

std::vector<Edge*>& Node::getEdges(){
	return this->edges;
}

Edge* Node::getEdge(int idx){
	return this->edges.at(idx);
}

void Node::setData(Tracker::Vortex& data){
	this->data = data;
}

void Node::addEdge(Edge &e){
	this->edges.push_back(&e);
}

void Node::removeEdge(unsigned int uid){
	Node *n;
	for (int ii=0; ii < this->Node::edges.size(); ++ii){
		if(this->Node::getEdge(ii)->getUid() == uid){
			n = this->Node::getConnectedNode(this->Node::getEdge(ii));
			for(int jj=0; jj<n->getEdges().size(); ++jj){
				if(n->getEdge(jj)->getUid() == uid) {
					Edge *e = n->getEdge(jj);
					n->getEdges().erase(n->getEdges().begin() + jj);
					delete e;
					break;
				}
			}
			this->Node::getEdges().erase(this->Node::getEdges().begin()+ii);
			break;
		}
	}
}

void Node::removeEdge(Node & n) {
	for(Edge *e1 : this->Node::getEdges()){
		for(Edge *e2 : n.getEdges()){
			if (Node::getConnectedNode(e1)->getUid() == e2->getUid()){
				Node::removeEdge(e2->getUid());
				return;
			}
		}
	}
}

void Node::removeEdge(Edge &edge){
	Node::removeEdge(edge.getUid());
}

void Node::removeEdges(){
	for(int ii=0; ii<this->getEdges().size(); ++ii){
		this->Node::removeEdge(*this->Node::getEdge(ii));
	}
	this->Node::getEdges().clear();
}

#include <iostream>

Node* Node::getConnectedNode(Edge *e){
	std::cout << "e->getNode(0)->getUid()" << e->getNode(0)->getUid() << std::endl;
	std::cout << "e->getNode(1)->getUid()" << e->getNode(1)->getUid() << std::endl;
	std::cout << "   this->Node::getUid()" << this->Node::getUid() << std::endl;
	//exit(1);
	return (e->getNode(0)->getUid() != this->Node::getUid()) ? e->getNode(0) :  e->getNode(1) ;
}
