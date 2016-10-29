///@cond LICENSE
/*** operators.h - GPUE: Split Operator based GPU solver for Nonlinear
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
 *  @file    operators.h
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

#ifndef OPERATORS_H
#define OPERATORS_H

#include "../include/ds.h"
#include "../include/constants.h"
#include <sys/stat.h>
#include <unordered_map>
#include <boost/math/special_functions.hpp>

// function to find Bz, the curl of the gauge field
 /**
 * @brief       Finds Bz, the curl of the gauge field
 * @ingroup     data
 * @param       Grid simulation data
 * @param       gauge field Ax
 * @param       gauge field Ay
 * @return      Bz, the curl of A
 */
double *curl2d(Grid &par, double *Ax, double *Ay);

// UPDATE LIST LATER
 /**
 * @brief	determines K for the standard rotational case
 * @ingroup	data
 * @param	Grid simulation data
 * @param       location in x, y, z
 * @return	K at that location
 */
double rotation_K(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines K for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      K at that location
 */
double rotation_K_dimensionless(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines K for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      K at that location
 */
double rotation_gauge_K(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines V for the standard harmonic case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double harmonic_V(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines V for the standard harmonic case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double harmonic_V_dimensionless(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines V for the standard harmonic case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double harmonic_gauge_V(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAx for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double rotation_pAx(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double rotation_pAy(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAx for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double rotation_Ax(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double rotation_Ay(Grid &par, Op &opr, int i, int j, int k);


 /**
 * @brief       determines pAx for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double test_Ax(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double test_Ay(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard fiber2d case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double fiber2d_Ay(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard fiber2d case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double fiber2d_Ax(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines electric field around the fiber
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double LP01_E_squared(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
std::unordered_map<std::string, double> read_matlab_data(int index);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double dynamic_Ax(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double dynamic_Ay(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double dynamic_Az(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double* file_A(Grid &par, std::string filename);

// Function to check whether a file exists
std::string filecheck(std::string filename);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
void parse_equation(Grid par, std::string &equation, double &val, 
                    int i, int j, int k);

#endif
