//
// Created by Lee James O'Riordan on 11/08/15.
//

#ifndef GPUE_1_MANIP_H
#define GPUE_1_MANIP_H
#include<cuda.h>
#include <math.h>
#include "constants.h"

namespace WFC{
	void phaseWinding(double *phi, int winding, double *x, double *y, double dx, double dy, double posx, double posy, int dim);
	void applyPhase(double *phi, double2 *wfc, int dim);

}

#endif //GPUE_1_MANIP_H
