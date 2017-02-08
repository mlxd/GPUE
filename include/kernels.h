///@endcond
//##############################################################################
/**
 *  @file    kernels.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief GPU kernel definitions
 *
 *  @section DESCRIPTION
 *  Kernel definitions for all CUDA-enabled routines for solving GPE.
 */
//##############################################################################

#ifndef KERNELS_H
#define KERNELS_H
#include<stdio.h>

/**
* @brief	Indexing of threads on grid
* @ingroup	gpu
*/
unsigned int getGid3d3d();;

/**
* @brief	Indexing of blocks on device
* @ingroup	gpu
*/
__device__ unsigned int getBid3d3d();
/**
* @brief	Indexing of threads in a block on device
* @ingroup	gpu
*/
__device__ unsigned int getTid3d3d();

//##############################################################################
/**
* Helper functions for complex numbers
*/
//##############################################################################

/**
* @brief	Calculates magnitude of complex number. $|a + ib|$
* @ingroup	gpu
* @param	in Complex number
* @return	Magnitude of complex number
*/
__device__ double complexMagnitude(double2 in);
/**
* @brief	Return the squared magnitude of a complex number. $|(a+\textrm{i}b)*(a-\textrm{i}b)|$
* @ingroup	gpu
* @param	in Complex number
* @return	Absolute-squared complex number
*/
__device__ double complexMagnitudeSquared(double2 in);
/**
* @brief	Returns conjugate of the a complex number
* @ingroup	gpu
* @param	in Number to be conjugated
* @return	Conjugated complex number
*/
__device__ double2 conjugate(double2 in);
/**
* @brief	Multiply real scalar by a complex number
* @ingroup	gpu
* @param	scalar Scalar multiplier
* @param	comp Complex multiplicand
* @return	Result of scalar * comp
*/
__device__ double2 realCompMult(double scalar, double2 comp);

//##############################################################################
/**
 * Multiplication for linear, non-linear and phase-imprinting of the condensate.
 */
//##############################################################################

/**
* @brief	Kernel for complex multiplication
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplcation result
*/
__global__ void cMult(cufftDoubleComplex* in1, cufftDoubleComplex* in2, cufftDoubleComplex* out);

/**
* @brief	Kernel for multiplcation with real array and complex array
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplcation result
*/
__global__ void cMultPhi(cufftDoubleComplex* in1, double* in2, cufftDoubleComplex* out);

/**
* @brief	Kernel for complex multiplication with nonlinear density term
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplcation result
* @param	dt Timestep for evolution
* @param	mass Atomic species mass
* @param	omegaZ Trapping frequency along z-dimension
* @param	gState If performing real (1) or imaginary (0) time evolution
* @param	N Number of atoms in condensate
*/
__global__ void cMultDensity(double2* in1, double2* in2, double2* out, double dt, double mass, int gstate, double gDenConst);

//##############################################################################

/**
* @brief	Hold vortex at specified position. Not implemented. cMultPhi should implement required functionality.
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplcation result
*/
__global__ void pinVortex(cufftDoubleComplex* in1, cufftDoubleComplex* in2, cufftDoubleComplex* out);

//##############################################################################

/**
* @brief        Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in Complex field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecMult(double2 *in, double *factor, double2 *out);

/**
* @brief	Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup	gpu
* @param	in Complex field to be scaled (divided, not multiplied)
* @param	factor Scaling factor to be used
* @param	out Pass by reference output for result
*/
__global__ void scalarDiv(double2* in, double factor, double2* out);

/**
* @brief        Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in Complex field to be scaled (multiplied, not divided)
* @param        factor Scaling factor to be used
* @param        out Pass by reference output for result
*/
__global__ void scalarMult(double2* in, double factor, double2* out);

/**
* @brief        Complex field raised to a power
* @ingroup      gpu
* @param        in Complex field to be scaled (multiplied, not divided)
* @param        power parameter
* @param        out Pass by reference output for result
*/
__global__ void scalarPow(double2* in, double param, double2* out);

/**
* @brief        Conjugate of double2*.
* @ingroup      gpu
* @param        in Complex field to be conjugated
* @param        out Pass by reference output for result
*/
__global__ void vecConjugate(double2 *in, double2 *out);

/**
* @brief	Complex field scaling and renormalisation. Not implemented. Use scalarDiv
* @ingroup	gpu
*/
__global__ void scalarDiv1D(double2*, double2*);
/**
* @brief	Complex field scaling and renormalisation. Not implemented. Use scalarDiv
* @ingroup	gpu
*/
__global__ void scalarDiv2D(double2*, double2*);
/**
* @brief	Used as part of multipass to renormalise the wavefucntion
* @ingroup	gpu
* @param	in Complex field to be renormalised
* @param	dr Smallest area element of grid (dx*dy)
* @param	pSum GPU array used to store intermediate results during parallel summation
*/
__global__ void scalarDiv_wfcNorm(double2* in, double dr, double2* pSum, double2* out);

//##############################################################################

/**
* @brief	Not implemented
* @ingroup	gpu
* @param	in Input field
* @param	out Output values
*/
__global__ void reduce(double2* in, double* out);
/**
* @brief        Performs wavefunction renormalisation using parallel summation and applying scalarDiv_wfcNorm
* @ingroup      gpu
* @param        input Wavefunction to be renormalised
* @param        output Pass by reference return of renormalised wavefunction
* @param        pass Number of passes performed by routine
*/
__global__ void thread_test(double* input, double* output);
/**
* @brief	Performs wavefunction renormalisation using parallel summation and applying scalarDiv_wfcNorm
* @ingroup	gpu
* @param	input Wavefunction to be renormalised
* @param	output Pass by reference return of renormalised wavefunction
* @param	pass Number of passes performed by routine
*/
__global__ void multipass(double2* input, double2* output, int pass);

//##############################################################################

/**
* @brief	Calculates angular momentum. Not fully implemented. Handled in post-processing instead.
* @ingroup	gpu
* @param	omega Harmonic trap rotation frequency
* @param	dt Time-step for evolution
* @param	wfc Wavefunction
* @param	xpyypx L_z operator
* @param	out Output of calculation
*/
__global__ void angularOp(double omega, double dt, double2* wfc, double* xpyypx, double2* out);

// Kernel to perform 2d transposition
__global__ void transpose2d(double *indata, double *outdata);

__global__ void naivetranspose2d(int xDim, int yDim, 
                            const double *indata, double *outdata);

// Kernel to perform 2d transposition
__global__ void transpose2d2(const double2 *indata, double2 *outdata);

__global__ void naivetranspose2d2(int xDim, int yDim, 
                            const double2 *indata, double2 *outdata);

//##############################################################################
/**
 * Non-implemented functions.
 */
 //##############################################################################

/**
* @brief	Calculates energy of the current state during evolution. Not implemented.
* @ingroup	gpu
* @param	wfc Wavefunction
* @param	op Operator to calculate energy for.
* @param	dt Time-step for evolution
* @param	energy Energy result output
* @param	gnd_state Wavefunction
* @param	op_space Check if position space with non-linear term or not.
* @param	sqrt_omegaz_mass sqrt(omegaZ/mass), part of the nonlin interaction term.
*/
__global__ void energyCalc(double2 *wfc, double2 *op, double dt, double2 *energy, int gnd_state, int op_space, double sqrt_omegaz_mass, double gDenConst);
/**
* @brief	Performs bra-ket state multiplication. Not fully implemented.
* @ingroup	gpu
* @param	in1 Bra
* @param	in2 Ket
* @return	<Bra|Ket>
*/
inline __device__ double2 braKetMult(double2 in1, double2 in2);

//template<typename T> __global__ void pSumT(T* in1, T* output, int pass);

/**
* @brief	Performs parallel sum. Not verified. I use multipass instead.
* @ingroup	gpu
* @param	in1 That which must be summed
* @param	output That which has been summed
* @param	pass Number of passes
*/
__global__ void pSum(double* in1, double* output, int pass);

//template<double> __global__ void pSumT(double* in1, double* output, int pass);

#endif
