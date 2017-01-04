///@endcond
//##############################################################################
/**
 *  @file    tracker.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Vortex lattice position, orientation, and tracking.
 *
 *  @section DESCRIPTION
 *  These functions determine vortex positions, calculate intervortex separation
 *	at the max density, determine vortex lattice angle, and all useful routines
 *	to know where they are hiding.
 */
 //##############################################################################

#ifndef TRACKER_H
#define TRACKER_H
#ifdef __linux
	#include<omp.h>
#elif __APPLE__
#endif
#include<math.h>
#include<stdio.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include "vort.h"

//##############################################################################

/** See the source file for more info on functions.**/
namespace Tracker {

    /** Vortex is used to track specific individual vortices.
	 * coords tracks x,y positions.
	 * sign indicates direction of vortex rotation.
	 * wind indicates the unit charge of the vortex.
	 */


    /**
	 * Values from solving Ax=b for vortex least squares core finder.
	 */
    static const double lsq[3][4] = {{-0.5, 0.5,  -0.5, 0.5},
                                     {-0.5, -0.5, 0.5,  0.5},
                                     {0.75, 0.25, 0.25, -0.25}};


	/**
	* @brief	Find vortex locations in the condensate.
	* @ingroup	data
	* @param	marker Matrix for vortex locations to nearest grid point
	* @param	wfc Wavefunction
	* @param	radius Vortex search radius in condensate
	* @param	xDim Length of X dimension
	* @param	x X grid
	* @param	timestep Timestep in simulation
	* @return	Number of found vortices of either rotation direction
	*/
    int findVortex(int *marker, double2* wfc, double radius, int xDim, double *x, int timestep);

	/**
	* @brief	Accepts matrix of vortex locations as argument, returns array of x,y coordinates of locations and winding
	* @ingroup	data
	* @param	marker Matrix containing vortex grid locations
	* @param	vLocation Array to store vortex data
	* @param	xDim Length of X dimension
	* @param	wfc Wavefunction
	*/
    void vortPos(int *marker, struct Vtx::Vortex *vLocation, int xDim, double2 *wfc);

	/**
	* @brief	Accepts matrix of vortex locations as argument, returns array of x,y coordinates of locations and winding
	* @ingroup	data
	* @param	marker Matrix containing vortex grid locations
	* @param	vLocation Array to store vortex data
	* @param	xDim Length of X dimension
	* @param	wfc Wavefunction
	*/
    void olPos(int *marker, int2 *vLocation, int xDim);

	/**
	* @brief	Changes in vortex positions. Not implemented. See vort.py for current tracking.
	* @ingroup	data
	* @param	cMarker Current vortex location matrix
	* @param	pMarker Previous vortex location matrix
	* @param	x X grid
	* @param	tolerance Maximum change acceptable for a vortex to have moved
	* @param	numVortices Number of vortices
	* @param	xDim Length of X dimension
	* @return	Vortex struct pointer.
	*/
    struct Vtx::Vortex *vortPosDelta(int *cMarker, int2 *pMarker, double *x, double tolerance, int numVortices, int xDim);

	/**
	* @brief	Determines the most central vortex in the condensate
	* @ingroup	data
	* @param	cArray Array of vortices
	* @param	length Number of vortices
	* @param	xDim Length of X dimension
	* @return	Central vortex struct
	*/
    struct Vtx::Vortex vortCentre(struct Vtx::Vortex *cArray, int length, int xDim);

	/**
	* @brief	Determines the rotation angle of the vortex lattice
	* @ingroup	data
	* @param	vortCoords Array of vortices
	* @param	central Central vortex in lattice
	* @param	numVort Number of vortices counted
	* @return	$0 \leq \theta \leq 2\pi$, though between $0 \leq \theta \leq \pi\3$ is sufficient.
	*/
    double vortAngle(struct Vtx::Vortex *vortCoords, struct Vtx::Vortex central, int numVort);

	/**
	* @brief	Determines average inter-vortex separation about the condensate centre
	* @ingroup	data
	* @param	vArray Vortices
	* @param	centre Central vortex in lattice
	* @param	length Number of counted vortices
	* @return	Separation distance
	*/
    double vortSepAvg(struct Vtx::Vortex *vArray, struct Vtx::Vortex centre, int length);


    double sigVOL(int2 *vArr, int2 *opLatt, double *x, int numVort);

	/**
	* @brief	Finds optical lattice maxima locations. Deprecated.
	* @ingroup	data
	* @param	marker Matrix of lattice maxima locations
	* @param	V Optical lattice potential
	* @param	radius Search radius for maxima
	* @param	x X grid
	*/
    int findOLMaxima(int *marker, double *V, double radius, int xDim, double *x);

	/**
	* @brief	Ensures the vortices are tracked and arranged in the right order based on minimum distance between previous and current positions
	* @ingroup	data
	* @param	vCoordsC Current vortex locations
	* @param	vCoordsP Previous vortex locations
	* @param	length Number of vortices tracked
	*/
    void vortArrange(struct Vtx::Vortex *vCoordsC, struct Vtx::Vortex *vCoordsP, int length);

	/**
	* @brief	Checks given coordinate for phase singularity of +ve winding
	* @ingroup	data
	* @param	vLoc Vortex location
	* @param	wfc Wavefunction
	* @param	xDim Length of X dimension
	* @return	1 for singularity. 0 otherwise
	*/
    int phaseTest(int2 vLoc, double2 *wfc, int xDim);

	/**
	* @brief	Least-squares vortex code estimation. Loosely based on c42f's vortutils code
	* @ingroup	data
	* @param	vortCoords Array of vortices. Result returned in struct double coordinates
	* @param	wfc Wavefunction
	* @param	numVort Number of vortices
	* @param	xDim Length of X dimension
	*/
    void lsFit(struct Vtx::Vortex *vortCoords, double2 *wfc, int numVort, int xDim);
}

#endif
