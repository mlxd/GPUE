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
double rotation_K3d(Grid &par, Op &opr, int i, int j, int k);

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
double harmonic_V3d(Grid &par, Op &opr, int i, int j, int k);

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
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double rotation_pAz(Grid &par, Op &opr, int i, int j, int k);

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
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double rotation_Az(Grid &par, Op &opr, int i, int j, int k);

 /**
 * @brief       determines pAy for the standard rotational case
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
double constant_A(Grid &par, Op &opr, int i, int j, int k);

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
void file_A(std::string filename, double *A);

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

/*----------------------------------------------------------------------------//
* WFC
*-----------------------------------------------------------------------------*/

 /**
 * @brief       creates the initial wavefunction for 2d
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
cufftDoubleComplex standard_wfc_2d(Grid &par, double Phi,
                                   int i, int j, int k);

 /**
 * @brief       creates the initial wavefunction for 2d
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
cufftDoubleComplex standard_wfc_3d(Grid &par, double Phi,
                                   int i, int j, int k);

 /**
 * @brief       creates the initial wavefunction for 2d
 * @ingroup     data
 * @param       Grid simulation data
 * @param       location in x, y, z
 * @return      V at that location
 */
cufftDoubleComplex torus_wfc(Grid &par, double Phi,
                             int i, int j, int k);

#endif
