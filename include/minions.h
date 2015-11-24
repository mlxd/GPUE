///@cond LICENSE
/*** minions.h - GPUE: Split Operator based GPU solver for Nonlinear
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
 *  @file    minions.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Helper functions for evaluating and manipulating data.
 *
 *  @section DESCRIPTION
 *  Some useful functions for carrying out trivial tasks as needed.
 */
//##############################################################################

#ifndef MINIONS_H
#define MINIONS_H

#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include "tracker.h"

namespace Minions{
	/**
	* @brief	Calculates $|z|^2$, where $z \element \mathbb{C}$
	* @ingroup	data
	* @param	in Complex number
	* @return	$|z|^2$
	*/
	double psi2(double2 in);

	/**
	* @brief	Returns the minimum value in the array
	* @ingroup	data
	* @param	grid Array of values
	* @param	len Length of grid
	* @return	$\min{\mathrm{grid}}$
	*/
	double minValue(double* grid,int len);
	/**
	* @brief	Returns the maximum value in the array
	* @ingroup	data
	* @param	grid Array of values
	* @param	len Length of grid
	* @return	$\max{\mathrm{grid}}$
	*/
	double maxValue(double* grid,int len);

	/**
	* @brief	Computes average of the array
	* @ingroup	data
	* @param	in Array of values
	* @param	len Length of array
	* @return	$\frac{1}{\mathrm{len}}\displaystyle\sum\limits_{i=1}^{\mathrm{len}}$
	*/
	double sumAvg(double* in, int len);

	/** id magic hackery */
	/**
	* @brief	id magic hackery. Double precision fast inverse square-root. Useless, but necessary to have.
	* @ingroup	data
	* @param	in Value to calculate fast inverse square root
	* @return	$\sqrt{\mathrm{in}}^{-1}$, but really fast
	*/
	double fInvSqRt(double in);
	//float fInvSqRt(float);

	/**
	* @brief	Swap the position of vortices
	* @ingroup	data
	* @param	vCoords Pointer to vortex array
	* @param	src Source vortex to swap
	* @param	dest Destination for vortex swap
	*/
	void coordSwap(struct Vtx::Vortex *vCoords, int src, int dest);

	/** More complex helper functions **/
	/**
	* @brief	Calculates $|z|$, where $z \element \mathbb{C}$
	* @ingroup	data
	* @param	in Complex number
	* @return	$|z|^2$
	*/
	double complexMag(double2 in);
	/**
	* @brief	Calculates $|z|^2$, where $z \element \mathbb{C}$
	* @ingroup	data
	* @param	in Complex number
	* @return	$|z|^2$
	*/
	double complexMag2(double2 in);
	/**
	* @brief	Calculates complex multiplication of input parameters
	* @ingroup	data
	* @param	in1 Complex number 1
	* @param	in1 Complex number 2
	* @return	in1 * in2
	*/
	double2 complexMult(double2 in1, double2 in2);
	/**
	* @brief	Calculates real * complex
	* @ingroup	data
	* @param	comp Complex number multiplicand
	* @param	scale Mulitplier
	* @return	scale * comp
	*/
	double2 complexScale(double2 comp, double scale);
	/**
	* @brief	Calculates complex conjugate
	* @ingroup	data
	* @param	c Complex number
	* @return	conj(c)
	*/
	double2 conj(double2 c);
	/**
	* @brief	Calculates complex division
	* @ingroup	data
	* @param	num Complex numerator
	* @param	num Complex denominator
	* @return	num/den
	*/
	double2 complexDiv(double2 num, double2 den);
}

#endif
