
//######################################################################################################################

#include "../include/edge.h"
#include <memory>
#include <iostream>

using namespace LatticeGraph;

//######################################################################################################################
//####################################    Ceiling Cat & Basement Cat     ###############################################
//######################################################################################################################

Edge::~Edge(){
}

Edge::Edge() : uid(++suid){
}

Edge::Edge(std::weak_ptr<Node> n1, std::weak_ptr<Node> n2) : uid(++suid){
	this->n1 = n1;
	this->n2 = n2;
	this->direction = 0;
	this->weight = 0;
}

Edge::Edge(std::weak_ptr<Node> n1, std::weak_ptr<Node> n2, int dir, double weight) : uid(++suid){
	this->n1 = n1;
	this->n2 = n2;
	this->direction = dir;
	this->weight = weight;
}

//######################################################################################################################
//####################################            Get stuff              ###############################################
//######################################################################################################################

/***
 * Returns UID of Edge.
 */
unsigned int Edge::getUid(){
	return uid;
}

/***
 * Returns direction of Edge.
 */
int Edge::getDirection(){
	return this->direction;
}
/***
 * Get Node at index idx.
 */
std::weak_ptr<Node> Edge::getVortex(int idx){
	try{
		return (idx==0) ? n1 : n2;
	}
	catch(std::exception e){
		return std::weak_ptr<Node>();
	}
}

/***
 * Get weight of edge.
 */
double Edge::getWeight(){
	return this->weight;
}

//######################################################################################################################
//####################################             Set stuff             ###############################################
//######################################################################################################################

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

//######################################################################################################################
//####################################           Update stuff            ###############################################
//######################################################################################################################

/***
 * Replaces Node n1 or n2 with new Node n_new.
 */
void Edge::updateVortex(int idx, std::weak_ptr<Node> n_new ){
	if(idx>0){
		this->n1 = n_new;
	}
	else{
		this->n2 = n_new;
	}
}

//######################################################################################################################
//####################################           Check stuff             ###############################################
//######################################################################################################################

bool Edge::isMember(std::weak_ptr<Node> n){
	return ( this->n1.lock()->getUid() == n.lock()->getUid() || this->n2.lock()->getUid() == n.lock()->getUid() ) ? true : false;
}

//######################################################################################################################
