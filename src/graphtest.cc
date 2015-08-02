//
// Created by Lee James O'Riordan on 23/07/15.
//

#define NUM_VORT 3

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
	Edge *e1,*e2,*e3;

	Tracker::Vortex *v = (Tracker::Vortex*) malloc(sizeof(Tracker::Vortex)*NUM_VORT);
	for(int i=0; i < NUM_VORT; ++i) {
		v[i].coords.x = (i + 1) * NUM_VORT;
		v[i].coords.y = (i + 1) * NUM_VORT;
		n = new Node(v[i]);
		l->Lattice::addNode(*n);
		//std::cout << "UID=" << l->getElement(i).getUid() << std::endl;
	}
	n=NULL;
/*
	e = new Edge(l->getElement(0),l->getElement(1),0,0 );
	l->addEdge(*e,l->getElement(0),l->getElement(1));

	e = new Edge(l->getElement(0),l->getElement(2),0,0 );
	l->addEdge(*e,l->getElement(0),l->getElement(2));

	//Edge e(l->getElement(0),l->getElement(3),0,0 );
	e = new Edge(l->getElement(0),l->getElement(3),0,0 );
	l->addEdge(*e,l->getElement(0),l->getElement(3));

	e = new Edge(l->getElement(0),l->getElement(4),0,0 );
	l->addEdge(*e,l->getElement(0),l->getElement(4));

	e = new Edge(l->getElement(0),l->getElement(5),0,0 );
	l->addEdge(*e,l->getElement(0),l->getElement(5));
*/

	e1 = new Edge(l->getElement(0),l->getElement(1),0,0 );
	l->addEdge(*e1,l->getElement(0),l->getElement(1));

	e2 = new Edge(l->getElement(0),l->getElement(2),0,0 );
	l->addEdge(*e2,l->getElement(0),l->getElement(2));

	e3 = new Edge(l->getElement(1),l->getElement(2),0,0 );
	l->addEdge(*e3,l->getElement(1),l->getElement(2));



	unsigned int* mat = (unsigned int*) calloc(l->getVortices().size()*l->getVortices().size(),sizeof(int));
	l->genAdjMat(mat);
	l->adjMatMtca(mat);

	free(mat); free(n); free(e); free(v);
}
