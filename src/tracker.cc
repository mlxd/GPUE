/*** tracker.cc - GPUE: Split Operator based GPU solver for Nonlinear 
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

#include "../include/tracker.h"
#include "../include/fileIO.h"
#include "../include/minions.h"
#include "../include/constants.h"

/**
 *  Contains all the glorious info you need to track vortices and see what they are up to.
 **/
namespace Tracker {
	char bufferT[1024];

	/**
	 * Determines the vortex separation at the centre of the lattice.
	 */
	double vortSepAvg(struct Vortex *vArray, struct Tracker::Vortex centre, int length){
		double result=0.0;// = sqrt( pow(centre.x - v_array[0].x,2) + pow(centre.y - v_array[0].y,2));
		double min = 0.0;
		double min_tmp = 0.0;
		int index=0;
		min = sqrt( pow(centre.coords.x - vArray[0].coords.x,2) + pow(centre.coords.y - vArray[0].coords.y,2));
		for (int j=1; j<length; ++j){
			min_tmp	= sqrt( pow(centre.coords.x - vArray[j].coords.x,2) + pow(centre.coords.y - vArray[j].coords.y,2));
			if(min > min_tmp && min_tmp > 1e-7){
				min = min_tmp;
				index = j;
			}
		}
		return min;
	}
	
	/**
	 * Finds the maxima of the optical lattice. Deprecated.
	 */
	int findOLMaxima(int *marker, double *Vopt, double radius, int xDim, double* x){
		double gridValues[9];
		int2 mIndex[1024];
		int2 index;
		int i,j,found;
		found=0;
		for (i=1; i<xDim-1; ++i ){
			for(j=1; j<xDim-1;++j){
				if(sqrt(x[i]*x[i] + x[j]*x[j]) < radius){			
					gridValues[0] = Vopt[(i-1)*xDim + (j-1)];
					gridValues[1] = Vopt[(i-1)*xDim + j];
					gridValues[2] = Vopt[(i-1)*xDim + (j+1)];
					gridValues[3] = Vopt[i*xDim + (j-1)];
					gridValues[4] = Vopt[i*xDim + j];
					gridValues[5] = Vopt[i*xDim + (j+1)];
					gridValues[6] = Vopt[(i+1)*xDim + (j-1)];
					gridValues[7] = Vopt[(i+1)*xDim + j];
					gridValues[8] = Vopt[(i+1)*xDim + (j+1)];
					if(fabs((gridValues[4]-Minions::maxValue(gridValues,9))/gridValues[4]) <= 1e-7){
						//printf ("%d,%d\n",i,j);
						(marker)[i*xDim + j] = 1;
						index.x=i;
						index.y=j;
						mIndex[found] = index;
						++found;
					}
				}
			}
		}
		return found;	
	}

	#ifdef VORT_MIN
	int findVortex(int *marker, double2* wfc, double radius, int xDim, double* x, int timestep){
		double gridValues[9];
		int2 vIndex[1024];
		int2 index;
		int i,j,found;
		found=0;
	//	#pragma omp parallel for private(j)
		for (i=1; i<xDim-1; ++i ){
			for(j=1; j<xDim-1;++j){
				if(sqrt(x[i]*x[i] + x[j]*x[j]) < radius){			
					gridValues[0] = Minions::psi2(wfc[(i-1)*xDim + (j-1)]);
					gridValues[1] = Minions::psi2(wfc[(i-1)*xDim + j]);
					gridValues[2] = Minions::psi2(wfc[(i-1)*xDim + (j+1)]);
					gridValues[3] = Minions::psi2(wfc[(i)*xDim + (j-1)]);
					gridValues[4] = Minions::psi2(wfc[(i)*xDim + j]);
					gridValues[5] = Minions::psi2(wfc[(i)*xDim + (j+1)]);
					gridValues[6] = Minions::psi2(wfc[(i+1)*xDim + (j-1)]);
					gridValues[7] = Minions::psi2(wfc[(i+1)*xDim + j]);
					gridValues[8] = Minions::psi2(wfc[(i+1)*xDim + (j+1)]);
					if(fabs((gridValues[4]-Minions::minValue(gridValues,9))/gridValues[4]) < 1e-7){
						//printf ("%d,%d\n",i,j);
						(marker)[i*xDim + j] = 1;
						index.x=i;
						index.y=j;
						vIndex[found] = index;
						found++;
					}
				}
			}
		}
		return found;	
	}
	#else 
	/**
	 * Phase winding method to determine vortex positions. Calculates the phase around a loop and checks if ~ +/-2Pi. 
	 */
	int findVortex(int *marker, double2* wfc, double radius, int xDim, double *x, int timestep){
			double2 *g = (double2*) malloc(sizeof(double2)*4);
			double *phiDelta = (double*) malloc(sizeof(double)*4);
		int i,j,found;
		int cond_x, cond_y;
		cond_x = 0; cond_y = 0;
		found = 0;
		long rnd_value = 0;
		double sum = 0.0;
			for ( i=0; i < xDim-1; ++i ){
					for( j=0; j < xDim-1; ++j ){
							if(sqrt(x[i]*x[i] + x[j]*x[j]) < radius){
									g[0] = Minions::complexScale( Minions::complexDiv( wfc[i*xDim + j],     	wfc[(i+1)*xDim + j] ), 		( Minions::complexMag( wfc[(i+1)*xDim + j])		/ Minions::complexMag( wfc[i*xDim + j] )));
									g[1] = Minions::complexScale( Minions::complexDiv( wfc[(i+1)*xDim + j], 	wfc[(i+1)*xDim + (j+1)] ), 	( Minions::complexMag( wfc[(i+1)*xDim + (j+1)])	/ Minions::complexMag( wfc[(i+1)*xDim + j] )));	
									g[2] = Minions::complexScale( Minions::complexDiv( wfc[(i+1)*xDim + (j+1)],	wfc[i*xDim + (j+1)] ), 		( Minions::complexMag( wfc[i*xDim + (j+1)])		/ Minions::complexMag( wfc[(i+1)*xDim + (j+1)] )));		
									g[3] = Minions::complexScale( Minions::complexDiv( wfc[i*xDim + (j+1)], 	wfc[i*xDim + j] ), 			( Minions::complexMag( wfc[i*xDim + j])			/ Minions::complexMag( wfc[i*xDim + (j+1)] )));

					for (int k=0; k<4; ++k){
						phiDelta[k] = atan2( g[k].y, g[k].x );
						if(phiDelta[k] <= -PI){
							phiDelta[k] += 2*PI;
						}
					}
					sum = phiDelta[0] + phiDelta[1] + phiDelta[2] + phiDelta[3];
					rnd_value = lround(sum/(2*PI));
					if( sum >= 1.9*PI && cond_x <= 0 && cond_y <= 0){
						marker[i*xDim + j] = rnd_value;
						++found;
						sum = 0.0;
						cond_x = 2; cond_y = 2;
									}
					else if( sum <= -1.9*PI && cond_x <= 0 && cond_y <= 0 )  {
						marker[i*xDim + j] = -rnd_value;
						++found;
						sum = 0.0;
						cond_x = 2; cond_y = 2;

					}
					--cond_x;
					--cond_y;
							}
					}
			}
		return found;
	}
	#endif

	/** 
	 * Accepts matrix of vortex locations as argument, returns array of x,y coordinates of locations and first encountered vortex angle 
	 */
	void olPos(int *marker, int2 *olLocation, int xDim){
		int i,j;
		unsigned int counter=0;
		for(i=0; i<xDim; ++i){
			for(j=0; j<xDim; ++j){
				if((marker)[i*xDim + j] == 1){
					(olLocation)[ counter ].x=i;
					(olLocation)[ counter ].y=j;
					++counter;
				}
			}
		}
	}

	/**  
	 * Tests the phase winding of the wavefunction, looking for vortices 
	 */
	int phaseTest(int2 vLoc, double2* wfc, int xDim){
		int result = 0;
		double2 gridValues[4];
		double phiDelta[4];
		double sum=0.0;
		int i=vLoc.x, j=vLoc.y;
		gridValues[0] = Minions::complexScale( Minions::complexDiv(wfc[i*xDim + j],wfc[(i+1)*xDim + j]), 			(Minions::complexMag(wfc[(i+1)*xDim + j])	 / Minions::complexMag(wfc[i*xDim + j])));
		gridValues[1] = Minions::complexScale( Minions::complexDiv(wfc[(i+1)*xDim + j],wfc[(i+1)*xDim + (j+1)]), 	(Minions::complexMag(wfc[(i+1)*xDim + (j+1)])/ Minions::complexMag(wfc[(i+1)*xDim + j])));
		gridValues[2] = Minions::complexScale( Minions::complexDiv(wfc[(i+1)*xDim + (j+1)],wfc[i*xDim + (j+1)]), 	(Minions::complexMag(wfc[i*xDim + (j+1)])	 / Minions::complexMag(wfc[(i+1)*xDim + (j+1)])));
		gridValues[3] = Minions::complexScale( Minions::complexDiv(wfc[i*xDim + (j+1)],wfc[i*xDim + j]), 			(Minions::complexMag(wfc[i*xDim + j])		 / Minions::complexMag(wfc[i*xDim + (j+1)])));

		for (int k=0; k<4; ++k){
			phiDelta[k] = atan2(gridValues[k].y,gridValues[k].x);
					if(phiDelta[k] <= -PI){
						phiDelta[k] += 2*PI;
			}
		}
		sum = phiDelta[0] + phiDelta[1] + phiDelta[2] + phiDelta[3];
		if(sum >=1.8*PI){
			result = 1;
		}
		return result;
	}

	/** 
	 * Accepts matrix of vortex locations as argument, returns array of x,y coordinates of locations and first encountered vortex angle 
	 */
	void vortPos(int *marker, struct Tracker::Vortex *vLocation, int xDim, double2 *wfc){
		int i,j;
		unsigned int counter=0;
		for(i=0; i<xDim; ++i){
			for(j=0; j<xDim; ++j){
				if( abs((marker)[i*xDim + j]) >= 1){
					(vLocation)[ counter ].coords.x=i;
					(vLocation)[ counter ].coords.y=j;
					(vLocation)[ counter ].sign = ( signbit(abs(marker[i*xDim + j])) == 0 ) ? 1 : -1;
					(vLocation)[ counter ].wind = abs(marker[i*xDim + j]);
					++counter;
				}
			}
		}
	}

	/**
	 * Ensures the vortices are tracked and arranged in the right order based on minimum distance between previous and current positions
	 */
	void vortArrange(struct Tracker::Vortex *vCoordsC, struct Vortex *vCoordsP, int length){
		int dist, dist_t;
		int i, j, index;
		for ( i = 0; i < length; ++i ){
			dist = 0x7FFFFFFF; //arbitrary big value
			index = i;
			for ( j = i; j < length ; ++j){
				dist_t = ( (vCoordsP[i].coords.x - vCoordsC[j].coords.x)*(vCoordsP[i].coords.x - vCoordsC[j].coords.x) + (vCoordsP[i].coords.y - vCoordsC[j].coords.y)*(vCoordsP[i].coords.y - vCoordsC[j].coords.y) );
				if(dist > dist_t ){
					dist = dist_t;
					index = j;
				}
			}
			Minions::coordSwap(vCoordsC,index,i);
		}
	}

	/** 
	 * Determines the coords of the vortex closest to the central position. Useful for centering the optical lattice over v. lattice*
	*/
	struct Vortex vortCentre(struct Tracker::Vortex *cArray, int length, int xDim){
		int i, j, counter=0;
		int valX, valY;
		double valueTest, value = 0.0;
		valX = (cArray)[0].coords.x - ((xDim/2)-1);
		valY = (cArray)[0].coords.y - ((xDim/2)-1);
		value = sqrt( valX*valX + valY*valY );//Calcs the sqrt(x^2+y^2) from central position. try to minimise this value
		for ( i=1; i<length; ++i ){
			valX = (cArray)[i].coords.x - ((xDim/2)-1);
			valY = (cArray)[i].coords.y - ((xDim/2)-1);
			valueTest = sqrt(valX*valX + valY*valY);
			if(value > valueTest){ 
				value = valueTest;
				counter = i;
			}
		}
		return (cArray)[counter];
	}

	/** 
	 * Determines the angle of the vortex lattice relative to the x-axis
	 */
	double vortAngle(struct Tracker::Vortex *vortCoords, struct Vortex central, int numVort){
		int location;
		double sign=1.0;
		double minVal=1e300;//(pow(central.x - vortCoords[0].x,2) + pow(central.y - vortCoords[0].y,2));
		for (int i=0; i < numVort; ++i){
			if (minVal > (pow(central.coords.x - vortCoords[i].coords.x,2) + pow(central.coords.y - vortCoords[i].coords.y,2)) && abs(central.coords.x - vortCoords[i].coords.x) > 2e-6 && abs(central.coords.y - vortCoords[i].coords.y) > 2e-6){
				minVal = (pow(central.coords.x - vortCoords[i].coords.x,2) + pow(central.coords.y - vortCoords[i].coords.y,2));
				location = i;
			}
		}
		double ang=(fmod(atan2( (vortCoords[location].coords.y - central.coords.y), (vortCoords[location].coords.x - central.coords.x) ),PI/3));
		printf("Angle=%e\n",ang);
		return PI/3 - ang;
		
		//return PI/2 + fmod(atan2(vortCoords[location].y-central.y, vortCoords[location].x - central.x), PI/3);
		//return PI/2 - sign*acos( ( (central.x - vortCoords[location].x)*(central.x - vortCoords[location].x) ) / ( minVal*(central.x - vortCoords[location].x) ) );
	}

	/**
	 * Sigma of vortex lattice and optical lattice
	 */
	double sigVOL(struct Tracker::Vortex *vArr, int2 *opLatt, double *x, int numVort){
		double sigma = 0.0;
		double dx = abs(x[1]-x[0]);
		for (int i=0; i<numVort; ++i){
			sigma += pow( abs( sqrt( (vArr[i].coords.x - opLatt[i].x)*(vArr[i].coords.x - opLatt[i].x) + (vArr[i].coords.y - opLatt[i].y)*(vArr[i].coords.y - opLatt[i].y) )*dx),2);
		}
		sigma /= numVort;
		return sigma;
	}

}
