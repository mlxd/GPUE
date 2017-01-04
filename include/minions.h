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
