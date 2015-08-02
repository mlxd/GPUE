/*** lattice.cc - GPUE: Split Operator based GPU solver for Nonlinear 
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

#include "../include/lattice.h"
#include <iostream>

using namespace LatticeGraph;

Lattice::Lattice(){
	
}
Lattice::~Lattice(){
	this->getVortices().clear();
}

Node& Lattice::getElement(int idx){
	return getVortices().at(idx);
}

/***
 * Gets the location of the Node with UID uid.
 */
int Lattice::getElementUidIdx(unsigned int uid){
	for (int ii=0; ii< getVortices().size(); ++ii){
		if(this->Lattice::getElement(ii).getUid()== uid){
			return ii;
		}
	}
	return -1;
}

/***
 * Gets the the Node with UID uid. Assumes Node exists.
 */
Node& Lattice::getElementUid(unsigned int uid){
	for (int ii=0; ii < this->Lattice::getVortices().size(); ++ii){
		if(this->Lattice::getElement(ii).getUid()== uid){
			return this->Lattice::getElement(ii);
		}
	}
}

void Lattice::setElement(int idx, Node &n){
	this->Lattice::getVortices().at(idx) = n;
}

void Lattice::addNode(Node &n){
	this->Lattice::getVortices().push_back(n);
}

void Lattice::removeNode(Node &n){
	this->Lattice::removeNode(this->Lattice::getElementUid(n.getUid()));
}

void Lattice::removeNode(int idx){
	Node *n = &this->Lattice::getElement(idx);
	this->Lattice::getVortices().erase(this->Lattice::getVortices().begin() + idx);
	delete n;
}

void Lattice::addEdge(Edge& e, Node &n1, Node &n2){
	n1.addEdge(e);
	n2.addEdge(e);
}

void Lattice::createEdges(){
	for(int ii=0; ii< this->Lattice::getVortices().size(); ++ii){
		for(int jj=ii; ii< this->Lattice::getVortices().size(); ++jj){
			if(Lattice::getNodeDistance(this->getElement(ii),this->getElement(jj))) {
				Edge *e = new Edge ( this->getElement(ii), this->getElement(jj) );
				this->Lattice::addEdge(*e,this->getElement(ii),this->getElement(jj));
				e = NULL;
			}
		}
	}
}

std::vector< Node >& Lattice::getVortices(){
	return this->vortices;
}

void Lattice::removeEdge(Node &n1, Node &n2){
	n1.removeEdge(n2);
}

void Lattice::removeEdges(Node &n1){
	n1.removeEdges();
}

/**
 *
 * Problem with nodes not returning the correct connection.
 */
void Lattice::genAdjMat(unsigned int *mat){
	int idx1, idx2, idx;
	idx1 = 0; idx2 = 0; idx=0;
	for(Node n : this->Lattice::getVortices()){
		idx1=this->getElementUidIdx(n.getUid());
		for(Edge *e : n.getEdges()){
			idx2 = this->getElementUidIdx(n.getConnectedNode(e)->getUid());
			idx = idx1*this->Lattice::getVortices().size() + idx2;
			std::cout << "IDX1=" << idx1 << "  IDX2=" << n.getConnectedNode(e)->getUid() << std::endl;
			if(idx1 != idx2){
				mat[idx] = 1;
			}
		}
	}
}

/**
 *
 * Outputs adjacency matrix in format for copy/paste into Mathematica.
 */
void Lattice::adjMatMtca(unsigned int *mat){
	unsigned int size = this->Lattice::getVortices().size();
	std::cout << "{";
	for(int i = 0; i < size; ++i){
		std::cout << "{";
		for(int j = 0; j < size; ++j){
			std::cout << mat[i*size + j];
			if(j<size-1)
				std::cout <<",";
			else
				std::cout << "}";
		}
		if(i<size-1)
			std::cout <<",";
		std::cout << std::endl;
	}
	std::cout << "}" << std::endl;
}

double Lattice::getNodeDistance(Node &n1,Node &n2){
	return sqrt(n1.getData().coords.x - n2.getData().coords.x) + (n1.getData().coords.y - n2.getData().coords.y);
}