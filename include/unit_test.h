///@endcond
//##############################################################################
/**
 *  @file    evolution.h
 *  @author  James Ryan Schloss (leios)
 *  @date    5/31/2016
 *  @version 0.1
 *
 *  @brief function for evolution.
 *
 *  @section DESCRIPTION
 *  These functions and variables are necessary for carrying out the GPUE
 *	simulations. This file will be re-written in an improved form in some
 *	future release.
 */
//##############################################################################

#ifndef UNIT_H
#define UNIT_H

#include "ds.h"
/* Function declarations */
/*
 * arg1 = Function result code from CUDA CUFFT calls.
 * arg2 = String data for name of function called. Prints value to stdout.
 */

// UPDATE LIST LATER
 /**
 * @brief	performs all necessary unit tests to ensure proper function
 * @ingroup	data
 */
void test_all();

 /**
 * @brief        Performs a simple integral with trapezoidal sum
 * @ingroup      data
 * @param        array to be integrated
 * @param        dimension of the array to be integrated
 * @param        dx of the array (assuming we are integrating along x)
 */
double trapz(double *array, int dimension, double dx);

#endif
