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

//######################################################################################################################

#include "../include/lattice.h"
#include <iostream>

using namespace LatticeGraph;

//######################################################################################################################
//####################################    Ceiling Cat & Basement Cat     ###############################################
//######################################################################################################################

Lattice::Lattice(){
}

Lattice::~Lattice(){
	this->getVortices().clear();
	this->getEdges().clear();
}

//######################################################################################################################
//####################################            Get stuff              ###############################################
//######################################################################################################################

std::vector< std::shared_ptr<Node> >& Lattice::getVortices(){
	return this->vortices;
}

std::shared_ptr<Node> Lattice::getVortexIdx(unsigned int idx){
	return getVortices().at(idx);
}

/***
 * Gets the location of the Node with UID uid.
 */
unsigned int Lattice::getVortexIdxUid(unsigned int uid){
	for (int ii=0; ii< getVortices().size(); ++ii){
		if(this->Lattice::getVortexIdx(ii)->getUid()== uid){
			return ii;
		}
	}
	return -1;
}

/***
 * Gets the the Node with UID uid. Assumes Node exists.
 */
std::shared_ptr<Node> Lattice::getVortexUid(unsigned int uid){
	for (std::shared_ptr<Node> n : this->Lattice::getVortices()){
		if(n->getUid()== uid){
			return n;
		}
	}
	return std::shared_ptr<Node>();
}

double Lattice::getVortexDistance(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
	return sqrt(pow(n1->getData().coords.x - n2->getData().coords.x,2)
	            +  pow(n1->getData().coords.y - n2->getData().coords.y,2));
}

double Lattice::getVortexDistanceD(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
	return sqrt(pow(n1->getData().coordsD.x - n2->getData().coordsD.x,2)
	            +  pow(n1->getData().coordsD.y - n2->getData().coordsD.y,2));
}

std::shared_ptr<Edge> Lattice::getEdgeIdx(unsigned int idx){
	return getEdges().at(idx);
}

/***
 * Gets the location of the Edge with UID uid.
 */
unsigned int Lattice::getEdgeIdxUid(unsigned int uid){
	for (int ii=0; ii< getEdges().size(); ++ii){
		if(this->Lattice::getEdgeIdx(ii)->getUid()== uid){
			return ii;
		}
	}
	return -1;
}

/***
 * Gets the the Edge with UID uid. Assumes Node exists.
 */
std::shared_ptr<Edge> Lattice::getEdgeUid(unsigned int uid){
	for (std::shared_ptr<Edge> e : this->Lattice::getEdges()){
		if(e->getUid()== uid){
			return e;
		}
	}
	return NULL;
}

std::vector< std::shared_ptr<Edge> >& Lattice::getEdges(){
	return this->edges;
}

//######################################################################################################################
//####################################             Set stuff             ###############################################
//######################################################################################################################

void Lattice::setVortex(unsigned int idx, std::shared_ptr<Node> n){
	this->Lattice::getVortices().at(idx)=(n);
}

void Lattice::setEdge(unsigned int idx, std::shared_ptr<Edge> e){
	this->Lattice::getEdges().at(idx)=(e);
}

//######################################################################################################################
//####################################              + stuff              ###############################################
//######################################################################################################################


void Lattice::createEdges(unsigned int radius){
	std::shared_ptr<Edge> e;
	double dist = 0.0;
	for(int ii=0; ii< this->Lattice::getVortices().size(); ++ii){
		//std::cout << "Got here ii " << ii << std::endl;
		for(int jj=ii+1; jj < this->Lattice::getVortices().size(); ++jj){
			dist = Lattice::getVortexDistance(this->getVortexIdx(ii),this->getVortexIdx(jj));
			if(dist < radius ) {
				//std::cout << "Got here jj " << jj << std::endl;
				e.reset(new Edge ( this->getVortexIdx(ii), this->getVortexIdx(jj) ));
				e->setWeight(dist);
				this->Lattice::addEdge(e,this->getVortexIdx(ii),this->getVortexIdx(jj));
			}
		}
	}
}
void Lattice::createEdges(double radius){
	std::shared_ptr<Edge> e;
	double dist = 0.0;
	for(int ii=0; ii< this->Lattice::getVortices().size(); ++ii){
		//std::cout << "Got here ii " << ii << std::endl;
		for(int jj=ii+1; jj < this->Lattice::getVortices().size(); ++jj){
			dist = Lattice::getVortexDistance(this->getVortexIdx(ii),this->getVortexIdx(jj));
			if( dist < radius ) {
				//std::cout << "Got here jj " << jj << std::endl;
				e.reset(new Edge ( this->getVortexIdx(ii), this->getVortexIdx(jj) ));
				e->setWeight(dist);
				this->Lattice::addEdge(e,this->getVortexIdx(ii),this->getVortexIdx(jj));
			}
		}
	}
}


void Lattice::addVortex(std::shared_ptr<Node> n){
	this->Lattice::getVortices().push_back((n));
}

void Lattice::addEdge(std::shared_ptr<Edge> e){
	this->addEdge(e, e->getVortex(0).lock(), e->getVortex(1).lock());
}

void Lattice::addEdge(std::shared_ptr<Edge> e, std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
	this->Lattice::getEdges().push_back(e);
	std::weak_ptr<Edge> e1 = e;
	std::weak_ptr<Edge> e2 = e;
	n1->addEdge(e1);
	n2->addEdge(e2);
}

//######################################################################################################################
//####################################              - stuff              ###############################################
//######################################################################################################################

void Lattice::removeVortex(std::shared_ptr<Node> n){
	for(std::weak_ptr<Edge> e : n->getEdges()){
		if(e.lock()){
			std::cout << "UID: Removing Vortex{" << n->getUid() <<"}." << std::endl;
			this->removeEdge(e.lock());
			this->Lattice::getVortices().erase(this->Lattice::getVortices().begin() + this->getVortexIdxUid(n->getUid()));
		}
		else{
			std::cout << "Cannot remove UID:Edge{"<< e.lock()->getUid() << "}, does not exist." << std::endl;
		}
	}
}

void Lattice::removeVortexUid(unsigned int uid){
	auto vtx = this->getVortexUid(uid);
	if(vtx){
		this->Lattice::removeVortex(vtx);
	}
	else{
		std::cout << "Cannot remove UID:Vortex{"<< uid << "}, does not exist." << std::endl;
	}
}

void Lattice::removeVortexIdx(unsigned int idx){
	auto vtx = this->getVortexIdx(idx);
	if(vtx){
		this->Lattice::removeVortex(vtx);
	}
	else{
		std::cout << "Cannot remove IDX:Vortex["<< idx << "], does not exist." << std::endl;
	}
}


void Lattice::removeEdge(std::shared_ptr<Edge> e){
	std::cout << "Removing Edge{" << e->getUid() <<"} connecting Node{" << e->getVortex(0).lock()->getUid() << "} and Node{" << e->getVortex(1).lock()->getUid() << "}." << std::endl;
	e->getVortex(0).lock()->removeEdgeUid(e->getUid());
	e->getVortex(1).lock()->removeEdgeUid(e->getUid());
	this->Lattice::getEdges().erase(this->Lattice::getEdges().begin() + this->Lattice::getEdgeIdxUid(e->getUid()));
}

void Lattice::removeEdgeIdx(unsigned int idx){
	std::weak_ptr<Edge> e = this->getEdgeIdx(idx);
	if (auto el = e.lock()) {
		this->Lattice::removeEdge(el);
	}
	else{
		std::cout << "Cannot remove IDX:Edge[" << idx << "], does not exist." << std::endl;
	}
}

void Lattice::removeEdgeUid(unsigned int uid) {
	std::weak_ptr<Edge> e = this->getEdgeUid(uid);
	if (auto el = e.lock()) {
		this->Lattice::removeEdge(el);
	}
	else{
		std::cout << "Cannot remove UID:Edge{" << uid << "}, does not exist." << std::endl;
	}
}

void Lattice::removeEdge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
	std::weak_ptr<Edge> e = this->Lattice::isConnected(n1,n2);
	if(e.lock()){
		this->Lattice::removeEdge(e.lock());
	}
	else{
		std::cout << "Node{" << n1->getUid() << "} and Node{" << n2->getUid() << "} were unconnected." << std::endl;
	}

}

void Lattice::removeEdges(std::shared_ptr<Node> n1){
	//n1->removeEdges();
}


void Lattice::createVortex(double posx, double posy, int winding){

}

void Lattice::destroyVortex(unsigned int uid){
	this->Lattice::getVortexUid(uid);
}



//######################################################################################################################
//####################################         Generate stuff            ###############################################
//######################################################################################################################

/**
 * Create adjacency matrix
 */
void Lattice::genAdjMat(unsigned int *mat){
	int idx1, idx2, idx;
	idx1 = 0; idx2 = 0; idx=0;
	for(std::shared_ptr<Node> n : this->Lattice::getVortices()){
		idx1=this->getVortexIdxUid(n->getUid());
		for(std::weak_ptr<Edge> e : n->getEdges()){
			idx2 = this->getVortexIdxUid(n->getConnectedNode(e.lock())->getUid());
			//std::cout << "this=" << n->getUid() << "   connected=" << n->getConnectedNode(e.lock())->getUid() << std::endl;
			idx = idx1*this->Lattice::getVortices().size() + idx2;
			//std::cout << "idx1=" << idx1 << "   idx2=" << idx2 << " idx=" << idx << "\n" << std::endl;
			mat[idx] = 1;
		}
	}
}

void Lattice::genAdjMat(double *mat){
	int idx1, idx2, idx;
	idx1 = 0; idx2 = 0; idx=0;
	for(std::shared_ptr<Node> n : this->Lattice::getVortices()){
		idx1=this->getVortexIdxUid(n->getUid());
		for(std::weak_ptr<Edge> e : n->getEdges()){
			idx2 = this->getVortexIdxUid(n->getConnectedNode(e.lock())->getUid());
			//std::cout << "this=" << n->getUid() << "   connected=" << n->getConnectedNode(e.lock())->getUid() << std::endl;
			idx = idx1*this->Lattice::getVortices().size() + idx2;
			//std::cout << "idx1=" << idx1 << "   idx2=" << idx2 << " idx=" << idx << "\n" << std::endl;
			mat[idx] = this->Lattice::getVortexDistance(n, this->getVortexIdx(idx2));
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
	for(int ii = 0; ii < size; ++ii){
		std::cout << "{";
		for(int jj = 0; jj < size; ++jj){
			std::cout << mat[ii*size + jj];
			if(jj<size-1)
				std::cout <<",";
			else
				std::cout << "}";
		}
		if(ii<size-1)
			std::cout <<",";
		std::cout << std::endl;
	}
	std::cout << "}" << std::endl;
}
void Lattice::adjMatMtca(double *mat){
	unsigned int size = this->Lattice::getVortices().size();
	std::cout << "{";
	for(int ii = 0; ii < size; ++ii){
		std::cout << "{";
		for(int jj = 0; jj < size; ++jj){
			std::cout << mat[ii*size + jj];
			if(jj<size-1)
				std::cout <<",";
			else
				std::cout << "}";
		}
		if(ii<size-1)
			std::cout <<",";
		std::cout << std::endl;
	}
	std::cout << "}" << std::endl;
}

//######################################################################################################################
//####################################           Check stuff             ###############################################
//######################################################################################################################

std::weak_ptr<Edge> Lattice::isConnected(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){

	if(n1->getUid() != n2->getUid()){
		for(std::weak_ptr<Edge> e1 : n1->getEdges()){
			if(e1.lock()->isMember(n2)){
				return e1;
			}
		}
	}
	return std::weak_ptr<Edge> ();
}

//######################################################################################################################
//####################################          Modify stuff             ###############################################
//######################################################################################################################

void Lattice::swapIdxUid(unsigned int uid1, unsigned int uid2) {
	Lattice::swapIdx(this->getVortexIdxUid(uid1),this->getVortexIdxUid(uid2));
}
void Lattice::swapIdx(unsigned int idx1, unsigned int idx2) {
	std::swap(this->getVortices().at(idx1),this->getVortices().at(idx2));
}
//void Lattice::swapVort(std::shared_ptr<Node> v1, std::shared_ptr<Node> v2) {

//}

//######################################################################################################################
