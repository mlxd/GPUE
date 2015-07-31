//
// Created by Lee James O'Riordan on 23/07/15.
//

#define NUM_VORT 6

#include "../include/lattice.h"
#include "../include/node.h"
#include "../include/edge.h"
#include <iostream>

using namespace LatticeGraph;
unsigned int Edge::suid = 0;
unsigned int Node::suid = 0;

int main(){
	Lattice *l = new Lattice();

	Node *n;//
	Edge *e;

	Tracker::Vortex *v = (Tracker::Vortex*) malloc(sizeof(Tracker::Vortex)*NUM_VORT);
	for(int i=0; i < NUM_VORT; ++i){
		v[i].coords.x = (i+1)*NUM_VORT;
		v[i].coords.y = (i+1)*NUM_VORT;
		n = new Node(v[i]);
		l->Lattice::addNode(*n);
		//std::cout << "UID=" << l->getElement(i).getUid() << std::endl;
	}
/*
	for(int i=0; i < 99; i=i+1){
		e = new Edge(l->getElement(i),l->getElement(i+1),0,0 );
		l->addEdge(*e,l->getElement(i),l->getElement(i+1));
	}
*/	for(int i=1; i < NUM_VORT; i++){
		std::cout << i << std::endl;
		e = new Edge(l->getElement(0),l->getElement(i),0,0 );
		l->addEdge(*e,l->getElement(0),l->getElement(i));
		e = NULL;
	}
/*
	Edge *e = new Edge(l->getElement(0),l->getElement(1),0,0 );
	l->addEdge(*e,l->getElement(0),l->getElement(1));

	e = new Edge(l->getElement(2),l->getElement(3),0,0 );
	l->addEdge(*e,l->getElement(2),l->getElement(3));

	e = new Edge(l->getElement(4),l->getElement(5),0,0 );
	l->addEdge(*e,l->getElement(4),l->getElement(5));

	e = new Edge(l->getElement(6),l->getElement(7),0,0 );
	l->addEdge(*e,l->getElement(6),l->getElement(7));

	e = new Edge(l->getElement(8),l->getElement(9),0,0 );
	l->addEdge(*e,l->getElement(8),l->getElement(9));
 */

	for(int i=0; i<NUM_VORT; ++i){
		//std::cout << "UID_Node=" << l->getElementUid(i).getUid() << std::endl;
		for(Edge *a : l->getElement(i).getEdges()){
			std::cout << "UID_Element_Edge=" << l->getElement(i).getUid() << "_" << a->getUid() << std::endl;
		}
	}

	unsigned int* mat = (unsigned int*) calloc(l->getVortices().size()*l->getVortices().size(),sizeof(int));
	l->genAdjMat(mat);

	std::cout << "{";
	for(int i = 0; i < NUM_VORT; ++i){
		std::cout << "{";
		for(int j = 0; j < NUM_VORT; ++j){
			std::cout << mat[i*NUM_VORT + j];
			if(j<NUM_VORT-1)
				std::cout <<",";
			else
				std::cout << "}";
		}
		if(i<NUM_VORT-1)
			std::cout <<",";
		std::cout << std::endl;
	}
	std::cout << "}" << std::endl;
}
