
//##############################################################################

#include "../include/node.h"
#include <memory>
#include <iostream>

using namespace LatticeGraph;

//##############################################################################
//#####################    Ceiling Cat & Basement Cat     ######################
//##############################################################################

Node::Node():uid(++suid){
}

Node::~Node(){
    this->removeEdges();
}

Node::Node(Vtx::Vortex& data):uid(++suid){
    this->data = data;
}

//##############################################################################
//#####################            Get stuff              ######################
//##############################################################################

unsigned int Node::getUid(){
    return uid;
}

Vtx::Vortex& Node::getData(){
    return this->data;
}

std::vector<std::weak_ptr <Edge> >& Node::getEdges(){
    return this->edges;
}

std::weak_ptr<Edge> Node::getEdge(int idx) {
    return this->edges.at(idx);
}

std::shared_ptr<Node> Node::getConnectedNode(std::shared_ptr<Edge> e){
    return (e->getVortex(0).lock()->getUid() != this->Node::getUid()) ? e->getVortex(0).lock() :  e->getVortex(1).lock() ;
}

//##############################################################################
//#####################             Set stuff             ######################
//##############################################################################

void Node::setData(Vtx::Vortex& data){
    this->data = data;
}

//##############################################################################
//#####################             +/- stuff             ######################
//##############################################################################

void Node::addEdge(std::weak_ptr<Edge> e){
    this->edges.push_back(e);
}

void Node::removeEdgeUid(unsigned int uid){
    for (size_t ii=0; ii < this->Node::edges.size(); ++ii){
        if(this->Node::getEdge(ii).lock()->getUid() == uid){
            this->Node::getEdges().erase(this->Node::getEdges().begin()+ii);
            break;
        }
    }
}

void Node::removeEdgeIdx(unsigned int idx){
    this->Node::getEdges().erase(this->Node::getEdges().begin()+idx);
}

void Node::removeEdge(std::shared_ptr<Node> n) {
    for(std::weak_ptr<Edge> e1 : this->Node::getEdges()){
        for(std::weak_ptr<Edge> e2 : n->getEdges()){
            if (Node::getConnectedNode(e1.lock())->getUid() == e2.lock()->getUid()){
                this->Node::removeEdgeUid(e2.lock()->getUid());
                return;
            }
        }
    }
}

void Node::removeEdge(std::weak_ptr<Edge> edge){
    this->Node::removeEdgeUid(edge.lock()->getUid());
}

void Node::removeEdges(){
    for(size_t ii=0; ii<this->getEdges().size(); ++ii){
        this->Node::removeEdge(this->Node::getEdge(ii));
    }
    this->Node::getEdges().clear();
}

//##############################################################################
