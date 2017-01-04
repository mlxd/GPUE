///@endcond
//##############################################################################
/**
 *  @file    manip.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    11/08/2015
 *  @version 0.1
 *
 *  @brief Routines for manipulating the wavefunction.
 *
 *  @section DESCRIPTION
 *  The functions herein can be used to modify the wavefunction, in particular
 *  by applying phase windings. This can be used in combination with the
 *  vortex lattice graph to eliminate or add vortices at arbitrary positions.
 */
 //##############################################################################

#ifndef MANIP_H
#define MANIP_H
#include<cuda.h>
#include <math.h>
#include "constants.h"

namespace WFC {

	/**
	* @brief	Generate phase winding for condensate imprint.
	* @ingroup	wfc
	* @param	phi Array to pass the phase profile by reference
	* @param	winding How many units of circulation. $2\pi winding$
	* @param	x X-space grid
	* @param	y Y-space grid
	* @param	dx Increment along X
	* @param	dy Increment along Y
	* @param	posx Location of phase singularity (centre), x position
	* @param	posy Location of phase singularity (centre), y position
	* @param	dim Length of X and Y dimension. Assumes square grid.
	*/
    void phaseWinding(double *phi, int winding, double *x, double *y, double dx, double dy, double posx, double posy, int dim);

	/**
	* @brief	Generate phase winding for condensate imprint, multiple locations
	* @ingroup	wfc
	* @param	phi Array to pass the phase profile by reference
	* @param	winding How many units of circulation. $2\pi winding$
	* @param	x X-space grid
	* @param	y Y-space grid
	* @param	dx Increment along X
	* @param	dy Increment along Y
	* @param	posx Location of phase singularities (centre), X positions
	* @param	posy Location of phase singularities (centre), Y positions
	* @param	dim Length of X and Y dimension. Assumes square grid.
	* @todo		Not yet implemented!!!
	*/
    void phaseWinding(double *phi, int winding, double *x, double *y, double dx, double dy, double *posx, double *posy, int sites, int dim);

	/**
	* @brief	Applies the generated phase to the wavefunction
	* @ingroup	wfc
	* @param	phi Array to pass the phase profile by reference
	* @param	wfc Wavefunction to receive the phase imprint
	* @param	dim Length of X and Y dimension. Assumes square grid.
	*/
    void applyPhase(double *phi, double2 *wfc, int dim);
}

#endif //MANIP_H
