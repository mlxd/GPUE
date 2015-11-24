///@cond LICENSE
/*** manip.h - GPUE: Split Operator based GPU solver for Nonlinear
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
