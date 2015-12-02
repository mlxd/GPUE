/*
* split_op.cu - GPUE: Split Operator based GPU solver for Nonlinear 
Schrodinger Equation, Copyright (C) 2011-2015, Lee J. O'Riordan 
<loriordan@gmail.com>, Tadhg Morgan, Neil Crowley. All rights reserved.

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

#include "../include/split_op.h"
#include "../include/kernels.h"
#include "../include/constants.h"
#include "../include/fileIO.h"
#include "../include/tracker.h"
#include "../include/minions.h"
#include "../include/ds.h"

#include "../include/lattice.h"
#include "../include/node.h"
#include "../include/edge.h"
#include "../include/manip.h"
#include "../include/vort.h"
#include <iostream>

unsigned int LatticeGraph::Edge::suid = 0;
unsigned int LatticeGraph::Node::suid = 0;

char buffer[100];
int verbose; //Print more info. Not curently implemented.
int device; //GPU ID choice.
int kick_it; //Kicking mode: 0 = off, 1 = multiple, 2 = single
int graph=0; //Generate graph from vortex lattice.
double gammaY; //Aspect ratio of trapping geometry.
double omega; //Rotation rate of condensate
double timeTotal;
double angle_sweep; //Rotation angle of condensate relative to x-axis
Params *paramS;
Array params;
double x0_shift, y0_shift; //Optical lattice shift parameters.
double Rxy; //Condensate scaling factor.
double a0x, a0y; //Harmonic oscillator length in x and y directions
double sepMinEpsilon=0.0; //Minimum separation for epsilon.

/*
 * Checks CUDA routines have exitted correctly.
 */
int isError(int result, char* c){
	if(result!=0){printf("Error has occurred for method %s with return type %d\n",c,result);
		exit(result);
	}
	return result;
}
int initialise(double omegaX, double omegaY, int N){
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	unsigned int xD=1,yD=1,zD=1;
	threads = 128;
	unsigned int b = xDim*yDim/threads;  //number of blocks in simulation
	unsigned long long maxElements = 65536*65536ULL; //largest number of elements

	if( b < (1<<16) ){
		xD = b;
	}
	else if( (b >= (1<<16) ) && (b <= (maxElements)) ){
		int t1 = log(b)/log(2);
		float t2 = (float) t1/2;
		t1 = (int) t2;
		if(t2 > (float) t1){
			xD <<= t1;
			yD <<= (t1 + 1);
		}
		else if(t2 == (float) t1){
			xD <<= t1;
			yD <<= t1;
		}
	}
	else{
		printf("Outside range of supported indexing");
		exit(-1);
	}
	printf("Compute grid dimensions chosen as X=%d	Y=%d\n",xD,yD);
	
	grid.x=xD; 
	grid.y=yD; 
	grid.z=zD; 
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	unsigned int i,j; //Used in for-loops for indexing
	
	unsigned int gSize = xDim*yDim;
	double xOffset, yOffset;
	xOffset=0.0;//5.0e-6;
	yOffset=0.0;//5.0e-6;
	
	mass = 1.0; //Rb 87 mass, kg
	appendData(&params,"Mass",mass);
	a_s = 0.0;//4.67e-9;
	appendData(&params,"a_s",a_s);

	double sum = 0.0;

	a0x = sqrt(1.0/2.0);
	a0y = sqrt(1.0/2.0);
	appendData(&params,"a0x",a0x);
	appendData(&params,"a0y",a0y);
	
	xMax = 10*a0x;//10*bec_length;//6*Rxy*a0x;
	yMax = 10*a0y;//10*bec_length;//
	appendData(&params,"xMax",xMax);
	appendData(&params,"yMax",yMax);

	double pxMax, pyMax;
	pxMax = (PI/xMax)*(xDim>>1);
	pyMax = (PI/yMax)*(yDim>>1);
	appendData(&params,"pyMax",pyMax);
	appendData(&params,"pxMax",pxMax);
	
	dx = xMax/(xDim>>1);
	dy = yMax/(yDim>>1);
	appendData(&params,"dx",dx);
	appendData(&params,"dy",dy);
	
	double dpx, dpy;
	dpx = PI/(xMax);
	dpy = PI/(yMax);
	appendData(&params,"dpx",dpx);
	appendData(&params,"dpy",dpy);

	//printf("a0x=%e  a0y=%e \n dx=%e   dx=%e\n R_xy=%e\n",a0x,a0y,dx,dy,Rxy);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	//double *x,*y,*xp,*yp;
	x = (double *) malloc(sizeof(double) * xDim);
	y = (double *) malloc(sizeof(double) * yDim);
	xp = (double *) malloc(sizeof(double) * xDim);
	yp = (double *) malloc(sizeof(double) * yDim);

	/*
	 * R-space and K-space grids
	 */
	for(i=0; i<xDim/2; ++i){
		x[i] = -xMax + (i)*dx;		
		x[i + (xDim/2)] = (i)*dx;
		
		y[i] = -yMax + (i)*dy;		
		y[i + (yDim/2)] = (i)*dy;
		
		xp[i] = (i)*dpx;
		xp[i + (xDim/2)] = -pxMax + (i)*dpx;
		
		yp[i] = (i)*dpy;
		yp[i + (yDim/2)] = -pyMax + (i)*dpy;
	}

/*
	 % maximum values
	  pxmax=pi*Ngx/(2*xmax);
	   % spacing in position and momentum space
	    dx=2*xmax/Ngx;   dpx=2*pxmax/Ngx;
	     % grid vectors, position and momentum space
	      x=(1:Ngx)'*dx-xmax;   pxn=((1:Ngx)*dpx-pxmax)';   
	       % reordination needed for the fourier transform
	        px(Ngx/2+2:Ngx,1)=pxn(1:Ngx/2-1,1);   px(1:Ngx/2+1,1)=pxn(Ngx/2:Ngx,1);
*/	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	/* Initialise wavefunction, momentum, position, angular momentum, imaginary and real-time evolution operators . */
	Energy = (double*) malloc(sizeof(double) * gSize);
	r = (double *) malloc(sizeof(double) * gSize);
	Phi = (double *) malloc(sizeof(double) * gSize);
	wfc = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	wfc_backup = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * (gSize/threads));
	K = (double *) malloc(sizeof(double) * gSize);
	V = (double *) malloc(sizeof(double) * gSize);
	V_opt = (double *) malloc(sizeof(double) * gSize);
	GK = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	GV = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EK = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EV = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EV_opt = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	xPy = (double *) malloc(sizeof(double) * gSize);
	yPx = (double *) malloc(sizeof(double) * gSize);
	ExPy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EyPx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EappliedField = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	
	/* Initialise wfc, EKp, and EVr buffers on GPU */
	cudaMalloc((void**) &Energy_gpu, sizeof(double) * gSize);
	cudaMalloc((void**) &wfc_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &Phi_gpu, sizeof(double) * gSize);
	cudaMalloc((void**) &K_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &V_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &xPy_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &yPx_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &par_sum, sizeof(cufftDoubleComplex) * (gSize/threads));
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

	#ifdef __linux
	int cores = omp_get_num_procs();
	appendData(&params,"Cores_Total",cores);
	appendData(&params,"Cores_Max",cores/2); //Assuming dev system specifics (Xeon with HT -> cores detected / 2)
	omp_set_num_threads(cores/2);
	#pragma omp parallel for private(j)
	#endif
	for( i=0; i < xDim; i++ ){
		for( j=0; j < yDim; j++ ){
			Phi[(i*yDim + j)] = fmod(l*atan2(y[j], x[i]),2*PI);
			
			wfc[(i*yDim + j)].x =  (1/sqrt(2))*pow(1/PI,0.5) * exp( -0.5*( x[i]*x[i] + y[j]*y[j] ) )*(1+2*x[i]/sqrt(2));
			wfc[(i*yDim + j)].y = 0.0;
				
			V[(i*yDim + j)] = 0.5*( x[i]*x[i] + y[j]*y[j] );
			K[(i*yDim + j)] = 0.5*( xp[i]*xp[i] + yp[j]*yp[j] );

			GV[(i*yDim + j)].x = exp( -V[(i*xDim + j)]*(gdt/(2)));
			GK[(i*yDim + j)].x = exp( -K[(i*xDim + j)]*(gdt/1));
			GV[(i*yDim + j)].y = 0.0;
			GK[(i*yDim + j)].y = 0.0;
			
			xPy[(i*yDim + j)] = x[i]*yp[j];
			yPx[(i*yDim + j)] = -y[j]*xp[i];
			
			EV[(i*yDim + j)].x=cos( -V[(i*xDim + j)]*(dt/(2)));
			EV[(i*yDim + j)].y=sin( -V[(i*xDim + j)]*(dt/(2)));
			EK[(i*yDim + j)].x=cos( -K[(i*xDim + j)]*(dt/1));
			EK[(i*yDim + j)].y=sin( -K[(i*xDim + j)]*(dt/1));
			
			ExPy[(i*yDim + j)].x=cos(-omega*omegaX*xPy[(i*xDim + j)]*dt);
			ExPy[(i*yDim + j)].y=sin(-omega*omegaX*xPy[(i*xDim + j)]*dt);
			EyPx[(i*yDim + j)].x=cos(-omega*omegaX*yPx[(i*xDim + j)]*dt);
			EyPx[(i*yDim + j)].y=sin(-omega*omegaX*yPx[(i*xDim + j)]*dt);
	
			sum += (wfc[(i*xDim + j)].x*wfc[(i*xDim + j)].x + wfc[(i*xDim + j)].y*wfc[(i*xDim + j)].y);
		}
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	//hdfWriteDouble(xDim, V, 0, "V_0"); //HDF not required for current projects. Removed.
	//hdfWriteComplex(xDim, wfc, 0, "wfc_0");
	FileIO::writeOutDouble(buffer,"V",V,xDim*yDim,0);
	//FileIO::writeOutDouble(buffer,"V_opt",V_opt,xDim*yDim,0);
	FileIO::writeOutDouble(buffer,"K",K,xDim*yDim,0);
	FileIO::writeOutDouble(buffer,"xPy",xPy,xDim*yDim,0);
	FileIO::writeOutDouble(buffer,"yPx",yPx,xDim*yDim,0);
	FileIO::writeOut(buffer,"WFC",wfc,xDim*yDim,0);
	FileIO::writeOut(buffer,"ExPy",ExPy,xDim*yDim,0);
	FileIO::writeOut(buffer,"EyPx",EyPx,xDim*yDim,0);
	FileIO::writeOutDouble(buffer,"Phi",Phi,xDim*yDim,0);
	FileIO::writeOutDouble(buffer,"r",r,xDim*yDim,0);
	FileIO::writeOutDouble(buffer,"x",x,xDim,0);
	FileIO::writeOutDouble(buffer,"y",y,yDim,0);
	FileIO::writeOutDouble(buffer,"px",xp,xDim,0);
	FileIO::writeOutDouble(buffer,"py",yp,yDim,0);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

	//free(V); 
	free(K); free(r); //free(Phi);

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	if (gsteps>0.0){
		sum=sqrt(sum*dx*dy);
		//#pragma omp parallel for reduction(+:sum) private(j)
		for (i = 0; i < xDim; i++){
			for (j = 0; j < yDim; j++){
				wfc[(i*yDim + j)].x = (wfc[(i*yDim + j)].x)/(sum);
				wfc[(i*yDim + j)].y = (wfc[(i*yDim + j)].y)/(sum);
			}
		}
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	result = cufftPlan2d(&plan_2d, xDim, yDim, CUFFT_Z2Z);
	if(result != CUFFT_SUCCESS){
		printf("Result:=%d\n",result);
		printf("Error: Could not execute cufftPlan2d(%s ,%d, %d).\n", "plan_2d", (unsigned int)xDim, (unsigned int)yDim);
		return -1;
	}

	result = cufftPlan1d(&plan_1d, xDim, CUFFT_Z2Z, yDim);
	if(result != CUFFT_SUCCESS){
		printf("Result:=%d\n",result);
		printf("Error: Could not execute cufftPlan3d(%s ,%d ,%d ).\n", "plan_1d", (unsigned int)xDim, (unsigned int)yDim);
		return -1;
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	return 0;
}

int evolve( cufftDoubleComplex *gpuWfc, 
			cufftDoubleComplex *gpuMomentumOp,
			cufftDoubleComplex *gpuPositionOp,
			void *gpu1dyPx,
			void *gpu1dxPy,
			cufftDoubleComplex *gpuParSum,			 
			int gridSize, int numSteps, int threads, 
			unsigned int gstate, int lz, int nonlin, int printSteps, int N, unsigned int ramp){

	//Because no two operations are created equally. Multiplimultiplication is faster than divisions.
	double renorm_factor_2d=1.0/pow(gridSize,0.5);
	double renorm_factor_1d=1.0/pow(xDim,0.5);

	clock_t begin, end;
	double time_spent;
	double Dt;
	if(gstate==0){
		Dt = gdt;
		printf("Timestep for groundstate solver set as: %E\n",Dt);
	}
	else{
		Dt = dt;
		printf("Timestep for evolution set as: %E\n",Dt);
	}
	begin = clock();
	double omega_0=omega*omegaX;

	#if 0 
	/** Determines the initial average density at the condensate central position and calculates a value for the healing length from this. Used thereafter as the lower limit for distances between vortices. **/
	int gridSum = 1<<6;
	double *densitySubset = (double*) malloc(sizeof(double)*gridSum);
	#pragma omp parallel for private(k)
	for (int j=0; j<gridSum; ++j){
		for (int k=0; k<gridSum; ++k){
			densitySubset[j*gridSum + k] = Minions::psi2(wfc[ ( (yDim/2) - (gridSum/2) + j )*yDim  + ( (xDim/2)  - (gridSum/2) + k )]);
		}
	}
	xi = 1/sqrt(8*PI*a_s*Minions::sumAvg(densitySubset,gridSum)/(dx*dy));//defined central condensate density
	printf("Avg healing length at centre=%E\n",xi);
	#endif

	/** ** ####################################################################################################### ** **/
	/** **					        HERE BE DRAGONS OF THE MOST DANGEROUS KIND!			                           ** **/
	/** ** ####################################################################################################### ** **/

	//Double buffering and will attempt to thread free and calloc operations to hide time penalty. Or may not bother.
	int num_vortices[2] = {0,0};
	int num_latt_max = 0;
	int* vortexLocation; //binary matrix of size xDim*yDim, 1 for vortex at specified index, 0 otherwise
	int* olMaxLocation = (int*) calloc(xDim*yDim,sizeof(int));

	struct Vtx::Vortex central_vortex; //vortex closest to the central position
	double vort_angle; //Angle of vortex lattice. Add to optical lattice for alignment.
	struct Vtx::Vortex *vortCoords = NULL; //array of vortex coordinates from vortexLocation 1's

	struct Vtx::Vortex *vortCoordsP = NULL; //Previous array of vortex coordinates from vortexLocation 1's

	int2 *olCoords = NULL; //array of vortex coordinates from vortexLocation 1's
	int2 *vortDelta = NULL;

	LatticeGraph::Lattice lattice; //Vortex lattice graph.
	double* adjMat;
	
	double vortOLSigma=0.0;
	double sepAvg = 0.0;
	
	int num_kick = 0;
	double t_kick = (2*PI/omega_0)/(6*Dt);
	
	for(int i=0; i < numSteps; ++i){
		if ( ramp == 1 ){
			omega_0=omegaX*((omega-0.39)*((double)i/(double)(numSteps)) + 0.39); //Adjusts omega for the appropriate trap frequency.
		}
		if(i % printSteps == 0) { //Print-out at pre-determined rate. Vortex & wfc analysis performed here also.
			printf("Step: %d	Omega: %lf\n", i, omega_0 / omegaX);
			cudaMemcpy(wfc, gpuWfc, sizeof(cufftDoubleComplex) * xDim * yDim, cudaMemcpyDeviceToHost);
			end = clock();
			time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
			printf("Time spent: %lf\n", time_spent);
			char *fileName = "";
			printf("ramp=%d		gstate=%d	rg=%d		\n", ramp, gstate, ramp | (gstate << 1));
			switch (ramp | (gstate << 1)) {
				case 0: //Groundstate solver, constant Omega value.
					fileName = "wfc_0_const";
			        break;
				case 1: //Groundstate solver, ramped Omega value.
					fileName = "wfc_0_ramp";
			        break;
				case 2: //Real-time evolution, constant Omega value.
					fileName = "wfc_ev";
					if(gpe==0){
						break;
					}
			        vortexLocation = (int *) calloc(xDim * yDim, sizeof(int));
			        num_vortices[0] = Tracker::findVortex(vortexLocation, wfc, 1e-4, xDim, x, i);

			        if (i == 0) { //If initial step, locate vortices, least-squares to find exact centre, calculate lattice angle, generate optical lattice.
				        vortCoords = (struct Vtx::Vortex *) malloc(
						        sizeof(struct Vtx::Vortex) * (2 * num_vortices[0]));
				        vortCoordsP = (struct Vtx::Vortex *) malloc(
						        sizeof(struct Vtx::Vortex) * (2 * num_vortices[0]));
				        Tracker::vortPos(vortexLocation, vortCoords, xDim, wfc);
				        Tracker::lsFit(vortCoords, wfc, num_vortices[0], xDim);
				        central_vortex = Tracker::vortCentre(vortCoords, num_vortices[0], xDim);
				        vort_angle = Tracker::vortAngle(vortCoords, central_vortex, num_vortices[0]);
				        appendData(&params, "Vort_angle", vort_angle);
				        optLatSetup(central_vortex, V, vortCoords, num_vortices[0],
				                    vort_angle + PI * angle_sweep / 180.0, laser_power * HBAR * sqrt(omegaX * omegaY),
				                    V_opt, x, y);
				        sepAvg = Tracker::vortSepAvg(vortCoords, central_vortex, num_vortices[0]);
				        if (kick_it == 2) {
					        printf("Kicked it 1\n");
					        cudaMemcpy(V_gpu, EV_opt, sizeof(cufftDoubleComplex) * xDim * yDim, cudaMemcpyHostToDevice);
				        }
				        FileIO::writeOutDouble(buffer, "V_opt_1", V_opt, xDim * yDim, 0);
				        FileIO::writeOut(buffer, "EV_opt_1", EV_opt, xDim * yDim, 0);
				        appendData(&params, "Central_vort_x", (double) central_vortex.coords.x);
				        appendData(&params, "Central_vort_y", (double) central_vortex.coords.y);
				        appendData(&params, "Central_vort_winding", (double) central_vortex.wind);
				        appendData(&params, "Num_vort", (double) num_vortices[0]);
				        FileIO::writeOutParam(buffer, params, "Params.dat");
			        }
			        else if (num_vortices[0] > num_vortices[1]) {
				        printf("Number of vortices increased from %d to %d\n", num_vortices[1], num_vortices[0]);
				        Tracker::vortPos(vortexLocation, vortCoords, xDim, wfc);
				        Tracker::lsFit(vortCoords, wfc, num_vortices[0], xDim);
			        }
			        else {
				        Tracker::vortPos(vortexLocation, vortCoords, xDim, wfc);
				        Tracker::lsFit(vortCoords, wfc, num_vortices[0], xDim);
				        Tracker::vortArrange(vortCoords, vortCoordsP, num_vortices[0]);
			        }

			        if (graph == 1) {

				        for (unsigned int ii = 0; ii < num_vortices[0]; ++ii) {
					        std::shared_ptr<LatticeGraph::Node> n(new LatticeGraph::Node(vortCoords[ii]));
					        lattice.addVortex(std::move(n));
				        }
				        unsigned int *uids = (unsigned int *) malloc(
						        sizeof(unsigned int) * lattice.getVortices().size());
				        for (int a = 0; a < lattice.getVortices().size(); ++a) {
					        uids[a] = lattice.getVortexIdx(a)->getUid();
				        }
				        if(i==0) {
					        //Lambda for vortex annihilation/creation.
					        auto killIt=[&](int idx) {
					            WFC::phaseWinding(Phi, 1, x, y, dx, dy, lattice.getVortexUid(idx)->getData().coordsD.x,
					                          lattice.getVortexUid(idx)->getData().coordsD.y, xDim);
					            cudaMemcpy(Phi_gpu, Phi, sizeof(double) * xDim * yDim, cudaMemcpyHostToDevice);
					            cMultPhi <<<grid, threads>>> (gpuWfc, Phi_gpu, gpuWfc);
				            };
					        //killIt(44); //Kills vortex with UID 44


				        }
				        lattice.createEdges(1.5 * 2e-5 / dx);
				        adjMat = (double *) calloc(lattice.getVortices().size() * lattice.getVortices().size(),
				                                   sizeof(double));
				        lattice.genAdjMat(adjMat);
				        FileIO::writeOutAdjMat(buffer, "graph", adjMat, uids, lattice.getVortices().size(), i);
				        free(adjMat);
				        free(uids);
				        lattice.getVortices().clear();
				        lattice.getEdges().clear();
				        //exit(0);
			        }

			        FileIO::writeOutVortex(buffer, "vort_arr", vortCoords, num_vortices[0], i);
			        printf("Located %d vortices\n", num_vortices[0]);
			        printf("Sigma=%e\n", vortOLSigma);
			        free(vortexLocation);
			        num_vortices[1] = num_vortices[0];
			        memcpy(vortCoordsP, vortCoords, sizeof(int2) * num_vortices[0]);
			        //exit(1);
			        break;
				case 3:
					fileName = "wfc_ev_ramp";
			        break;
				default:
					break;
			}
			if (write_it) {
				FileIO::writeOut(buffer, fileName, wfc, xDim * yDim, i);
			}
			//printf("Energy[t@%d]=%E\n",i,energy_angmom(gpuPositionOp, gpuMomentumOp, dx, dy, gpuWfc,gstate));
/*			cudaMemcpy(V_gpu, V, sizeof(double)*xDim*yDim, cudaMemcpyHostToDevice);
			cudaMemcpy(K_gpu, K, sizeof(double)*xDim*yDim, cudaMemcpyHostToDevice);
			cudaMemcpy(V_gpu, , sizeof(double)*xDim*yDim, cudaMemcpyHostToDevice);
			cudaMemcpy(K_gpu, K, sizeof(double)*xDim*yDim, cudaMemcpyHostToDevice);
*/		}
	
	/** ** ####################################################################################################### ** **/
	/** ** ####################################################################################################### ** **/
	/** ** 							More F'n' Dragons!				       ** **/
	/** ** ####################################################################################################### ** **/
		if(i % ((int)t_kick+1) == 0 && num_kick<=6 && gstate==1 && kick_it == 1 ){
			cudaMemcpy(V_gpu, EV_opt, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
			++num_kick;
		}
	/** ** ####################################################################################################### ** **/

		/*
		 * U_r(dt/2)*wfc
		 */ 
		if(nonlin == 1){
			cMultDensity<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc,0.5*Dt,mass,omegaZ,gstate,N*interaction);
		}
		else {
			cMult<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc);
		}
				
		/*
		 * U_p(dt)*fft2(wfc)
		 */		
		result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD);
		scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc); //Normalise
		cMult<<<grid,threads>>>(gpuMomentumOp,gpuWfc,gpuWfc);
		result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE);
		scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc); //Normalise
		
		/*
		 * U_r(dt/2)*wfc
		 */	
		if(nonlin == 1){
			cMultDensity<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc,Dt*0.5,mass,omegaZ,gstate,N*interaction);
		}
		else {
			cMult<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc);
		}
		if( kick_it !=0 && ( (i % (int)(t_kick+1) == 0 && num_kick<=6 && gstate==1) || (kick_it >= 1 && i==0)) ){
			cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
			printf("Got here: Cuda memcpy EV into GPU\n");
		}
		/**************************************************************/
		/* Angular momentum xPy-yPx   */
		if(lz == 1){
			switch(i%2 | (gstate<<1)){
				case 0: //Groundstate solver, even step
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_xPy
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				angularOp<<<grid,threads>>>(omega_0, Dt, gpuWfc, (double*) gpu1dxPy, gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
			
				result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD); //2D forward
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE); //1D inverse to wfc_yPx
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				angularOp<<<grid,threads>>>(omega_0, Dt, gpuWfc, (double*) gpu1dyPx, gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_PxPy
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE); //2D Inverse
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
				break;
				
				case 1:	//Groundstate solver, odd step
				result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD); //2D forward
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE); //1D inverse to wfc_yPx
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				angularOp<<<grid,threads>>>(omega_0, Dt, gpuWfc, (double*) gpu1dyPx, gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_PxPy
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE); //2D Inverse
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
				
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_xPy
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				angularOp<<<grid,threads>>>(omega_0, Dt, gpuWfc, (double*) gpu1dxPy, gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				break;
				
				case 2: //Real time evolution, even step
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_xPy
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				cMult<<<grid,threads>>>(gpuWfc, (cufftDoubleComplex*) gpu1dxPy, gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
			
				result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD); //2D forward
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE); //1D inverse to wfc_yPx
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				cMult<<<grid,threads>>>(gpuWfc, (cufftDoubleComplex*) gpu1dyPx, gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_PxPy
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE); //2D Inverse
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
				break;
				
				case 3:	//Real time evolution, odd step
				result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD); //2D forward
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE); //1D inverse to wfc_yPx
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				cMult<<<grid,threads>>>(gpuWfc, (cufftDoubleComplex*) gpu1dyPx, gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_PxPy
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE); //2D Inverse
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
				
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_xPy
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				cMult<<<grid,threads>>>(gpuWfc, (cufftDoubleComplex*) gpu1dxPy, gpuWfc);
				result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
				scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
				break;
			
			}
		}
		/**************************************************************/
	
		if(gstate==0){
			parSum(gpuWfc, gpuParSum, xDim, yDim, threads);
		}
	}
	return 0;
}

/*
 * Used to perform parallel summation on WFC for normalisation.
 */
void parSum(double2* gpuWfc, double2* gpuParSum, int xDim, int yDim, int threads){
		int grid_tmp = xDim*yDim;
		int block = grid_tmp/threads;
		int thread_tmp = threads;
		int pass = 0;
		while((double)grid_tmp/threads > 1.0){
			if(grid_tmp == xDim*yDim){
				multipass<<<block,threads,threads*sizeof(double2)>>>(&gpuWfc[0],&gpuParSum[0],pass); 
			}
			else{
				multipass<<<block,thread_tmp,thread_tmp*sizeof(double2)>>>(&gpuParSum[0],&gpuParSum[0],pass);
			}
			grid_tmp /= threads;
			block = (int) ceil((double)grid_tmp/threads);
			pass++;
		}
		thread_tmp = grid_tmp;
		multipass<<<1,thread_tmp,thread_tmp*sizeof(double2)>>>(&gpuParSum[0],&gpuParSum[0], pass);
		scalarDiv_wfcNorm<<<grid,threads>>>(gpuWfc, dx*dy, gpuParSum, gpuWfc);
}

/**
** Matches the optical lattice to the vortex lattice. Moire super-lattice project.
**/
void optLatSetup(struct Vtx::Vortex centre, double* V, struct Vtx::Vortex *vArray, int num_vortices, double theta_opt, double intensity, double* v_opt, double *x, double *y){
	int i,j;
	double sepMin = Tracker::vortSepAvg(vArray,centre,num_vortices);
	sepMin = sepMin*(1 + sepMinEpsilon);
	appendData(&params,"Vort_sep",(double)sepMin);
	/*
	* Defining the necessary k vectors for the optical lattice
	*/
	double k_mag = ((2*PI/(sepMin*dx))/2)*(2/sqrt(3)); // Additional /2 as a result of lambda/2 period
	double2* k = (double2*) malloc(sizeof(double2)*3);
	appendData(&params,"kmag",(double)k_mag);
	k[0].x = k_mag * cos(0*PI/3 + theta_opt);
	k[0].y = k_mag * sin(0*PI/3 + theta_opt);
	k[1].x = k_mag * cos(2*PI/3 + theta_opt);
	k[1].y = k_mag * sin(2*PI/3 + theta_opt);
	k[2].x = k_mag * cos(4*PI/3 + theta_opt);
	k[2].y = k_mag * sin(4*PI/3 + theta_opt);
	
	double2 *r_opt = (double2*) malloc(sizeof(double2)*xDim);

/*	for (int ii = 0; ii < xDim; ++ii){
		r_opt[ii].x = 0.0 + (xDim/sepMin)*PI*(ii-centre.coords.x)/(xDim-1);
		r_opt[ii].y = 0.0 + (xDim/sepMin)*PI*(ii-centre.coords.y)/(yDim-1);
	}
*/
	FileIO::writeOut(buffer,"r_opt",r_opt,xDim,0);
	appendData(&params,"k[0].x",(double)k[0].x);
	appendData(&params,"k[0].y",(double)k[0].y);
	appendData(&params,"k[1].x",(double)k[1].x);
	appendData(&params,"k[1].y",(double)k[1].y);
	appendData(&params,"k[2].x",(double)k[2].x);
	appendData(&params,"k[2].y",(double)k[2].y);

	double x_shift = dx*(9+(0.5*xDim-1) - centre.coords.x);//sin(theta_opt)*(sepMin);
	double y_shift = dy*(0+(0.5*yDim-1) - centre.coords.y);//cos(theta_opt)*(sepMin);

	printf("Xs=%e\nYs=%e\n",x_shift,y_shift);

	//#pragma omp parallel for private(j)
	for ( j=0; j<yDim; ++j ){
		for ( i=0; i<xDim; ++i ){
			v_opt[j*xDim + i] = intensity*(
					    	  pow( ( cos( k[0].x*( x[i] + x_shift ) + k[0].y*( y[j] + y_shift ) ) ), 2)
					  	+ pow( ( cos( k[1].x*( x[i] + x_shift ) + k[1].y*( y[j] + y_shift ) ) ), 2)
					  	+ pow( ( cos( k[2].x*( x[i] + x_shift ) + k[2].y*( y[j] + y_shift ) ) ), 2)
			/*		    	  pow( abs( cos( k[0].x*( r_opt[i].x + x_shift ) + k[0].y*( r_opt[j].y + y_shift ) ) ), 2)
					  	+ pow( abs( cos( k[1].x*( r_opt[i].x + x_shift ) + k[1].y*( r_opt[j].y + y_shift ) ) ), 2)
					  	+ pow( abs( cos( k[2].x*( r_opt[i].x + x_shift ) + k[2].y*( r_opt[j].y + y_shift ) ) ), 2)
			*/		    );
			EV_opt[(j*xDim + i)].x=cos( -(V[(j*xDim + i)] + v_opt[j*xDim + i])*(dt/(2*HBAR)));
			EV_opt[(j*xDim + i)].y=sin( -(V[(j*xDim + i)] + v_opt[j*xDim + i])*(dt/(2*HBAR)));
		}
	}
}

/**
** Calculates energy and angular momentum of current state. Implementation not fully finished.
**/
double energy_angmom(double *Energy, double* Energy_gpu, double2 *V_op, double2 *K_op, double dx, double dy, double2 *gpuWfc, int gState){
	double renorm_factor_2d=1.0/pow(xDim*yDim,0.5);
	double result=0;

	for (int i=0; i < xDim*yDim; ++i){
		Energy[i] = 0.0; 
	}
	
	
/*	cudaMalloc((void**) &energy_gpu, sizeof(double2) * xDim*yDim);

	energyCalc<<<grid,threads>>>( gpuWfc, V_op, 0.5*dt, energy_gpu, gState,1,i 0.5*sqrt(omegaZ/mass));
	result = cufftExecZ2Z( plan_2d, gpuWfc, gpuWfc, CUFFT_FORWARD );
	scalarDiv<<<grid,threads>>>( gpuWfc, renorm_factor_2d, gpuWfc ); //Normalise

	energyCalc<<<grid,threads>>>( gpuWfc, K_op, dt, energy_gpu, gState,0, 0.5*sqrt(omegaZ/mass));
	result = cufftExecZ2Z( plan_2d, gpuWfc, gpuWfc, CUFFT_INVERSE );
	scalarDiv<<<grid,threads>>>( gpuWfc, renorm_factor_2d, gpuWfc ); //Normalise
	
	err=cudaMemcpy(energy, energy_gpu, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyDeviceToHost);
	
	for(int i=0; i<xDim*yDim; i++){
		result += energy[i].x;
		//printf("En=%E\n",result*dx*dy);
	}
	return result*dx*dy;
*/
	
}


//###################################################################################################################
//###################################################################################################################

/*
 * Used to perform parallel summation using templates from c++
 */
template<typename T> void parSum(T *gpuToSumArr, T *gpuParSum, int xDim, int yDim, int threads){
                int grid_tmp = xDim*yDim;
                int block = grid_tmp/threads;
                int thread_tmp = threads;
                int pass = 0;
                while((double)grid_tmp/threads > 1.0){
                        if(grid_tmp == xDim*yDim){
                                multipass<<<block,threads,threads*sizeof(T)>>>(&gpuToSumArr[0],&gpuParSum[0],pass);
                        }
                        else{
                                multipass<<<block,thread_tmp,thread_tmp*sizeof(T)>>>(&gpuParSum[0],&gpuParSum[0],pass);
                        }
                        grid_tmp /= threads;
                        block = (int) ceil((double)grid_tmp/threads);
                        pass++;
                }
                thread_tmp = grid_tmp;
                multipass<<<1,thread_tmp,thread_tmp*sizeof(double2)>>>(&gpuParSum[0],&gpuParSum[0], pass);
                scalarDiv_wfcNorm<<<grid,threads>>>(gpuToSumArr, dx*dy, gpuParSum, gpuToSumArr);
}
//###################################################################################################################
//###################################################################################################################
int parseArgs(int argc, char** argv){
	int opt;
	while ((opt = getopt (argc, argv, "d:x:y:w:G:g:e:T:t:n:p:r:o:L:l:s:i:P:X:Y:O:k:W:U:V:S:a:")) != -1) {
		switch (opt)
		{
			case 'x':
				xDim = atoi(optarg);
				printf("Argument for x is given as %d\n",xDim);
				appendData(&params,"xDim",(double)xDim);
				break;
			case 'y':
				yDim = atoi(optarg);
				printf("Argument for y is given as %d\n",yDim);
				appendData(&params,"yDim",(double)yDim);
				break;
			case 'w':
				omega = atof(optarg);
				printf("Argument for OmegaRotate is given as %E\n",omega);
				appendData(&params,"omega",omega);
				break;
			case 'G':
				gammaY = atof(optarg);
				printf("Argument for gamma is given as %E\n",gammaY);
				appendData(&params,"gammaY",gammaY);
				break;
			case 'g':
				gsteps = atof(optarg);
				printf("Argument for Groundsteps is given as %ld\n",gsteps);
				appendData(&params,"gsteps",gsteps);
				break;
			case 'e':
				esteps = atof(optarg);
				printf("Argument for EvSteps is given as %ld\n",esteps);
				appendData(&params,"esteps",esteps);
				break;
			case 'T':
				gdt = atof(optarg);
				printf("Argument for groundstate Timestep is given as %E\n",gdt);
				appendData(&params,"gdt",gdt);
				break;
			case 't':
				dt = atof(optarg);
				printf("Argument for Timestep is given as %E\n",dt);
				appendData(&params,"dt",dt);
				break;
			case 'd':
				device = atoi(optarg);
				printf("Argument for device is given as %d\n",device);
				appendData(&params,"device",device);
				break;
			case 'n':
				atoms = atof(optarg);
				printf("Argument for atoms is given as %ld\n",atoms);
				appendData(&params,"atoms",atoms);
				break;
			case 'r':
				read_wfc  = atoi(optarg);
				printf("Argument for ReadIn is given as %d\n",read_wfc);
				appendData(&params,"read_wfc",(double)read_wfc);
				break;
			case 'p':
				print = atoi(optarg);
				printf("Argument for Printout is given as %d\n",print);
				appendData(&params,"print_out",(double)print);
				break;
			case 'L':
				l = atof(optarg);
				printf("Vortex winding is given as : %E\n",l);
				appendData(&params,"winding",l);
				break;
			case 'l':
				ang_mom = atoi(optarg);
				printf("Angular Momentum mode engaged: %d\n",ang_mom);
				appendData(&params,"corotating",(double)ang_mom);
				break;
			case 's':
				gpe = atoi(optarg);
				printf("Non-linear mode engaged: %d\n",gpe);
				appendData(&params,"gpe",gpe);
				break;
			case 'o':
				omegaZ = atof(optarg);
				printf("Argument for OmegaZ is given as %E\n",omegaZ);
				appendData(&params,"omegaZ",omegaZ);
				break;
			case 'i':
				interaction = atof(optarg);
				printf("Argument for interaction scaling is %E\n",interaction);
				appendData(&params,"int_scaling",interaction);
				break;
			case 'P':
				laser_power = atof(optarg);
				printf("Argument for laser power is %E\n",laser_power);
				appendData(&params,"laser_power",laser_power);
				break;
			case 'X':
				omegaX = atof(optarg);
				printf("Argument for omegaX is %E\n",omegaX);
				appendData(&params,"omegaX",omegaX);
				break;
			case 'Y':
				omegaY = atof(optarg);
				printf("Argument for omegaY is %E\n",omegaY);
				appendData(&params,"omegaY",omegaY);
				break;
			case 'O':
				angle_sweep = atof(optarg);
				printf("Argument for angle_sweep is %E\n",angle_sweep);
				appendData(&params,"angle_sweep",angle_sweep);
				break;
			case 'k':
				kick_it = atoi(optarg);
				printf("Argument for kick_it is %i\n",kick_it);
				appendData(&params,"kick_it",kick_it);
				break;
			case 'W':
				write_it = atoi(optarg);
				printf("Argument for write_it is %i\n",write_it);
				appendData(&params,"write_it",write_it);
				break;
			case 'U':
				x0_shift = atof(optarg);
				printf("Argument for x0_shift is %lf\n",x0_shift);
				appendData(&params,"x0_shift",x0_shift);
				break;
			case 'V':
				y0_shift = atof(optarg);
				printf("Argument for y0_shift is %lf\n",y0_shift);
				appendData(&params,"y0_shift",y0_shift);
				break;
			case 'S':
				sepMinEpsilon = atof(optarg);
				printf("Argument for sepMinEpsilon is %lf\n",sepMinEpsilon);
				appendData(&params,"sepMinEpsilon",sepMinEpsilon);
				break;
			case 'a':
				graph = atoi(optarg);
				printf("Argument for graph is %d\n",graph);
				appendData(&params,"graph",graph);
				break;
			case '?':
				if (optopt == 'c') {
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				} else if (isprint (optopt)) {
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				} else {
					fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
				}
				return -1;
			default:
				abort ();
		}
	}
	return 0;
}

void delta_define(double *x, double *y, double x0, double y0, double *delta){
	for (unsigned int i=0; i<xDim; ++i){
		for (unsigned int j=0; j<yDim; ++j){
			delta[j*xDim + i] = 1e6*HBAR*exp( -( pow( x[i] - x0, 2)  +  pow( y[j] - y0, 2) )/(5*dx*dx) );
			EV_opt[(j*xDim + i)].x=cos( -(V[(j*xDim + i)] + delta[j*xDim + i])*(dt/(2*HBAR)));
			EV_opt[(j*xDim + i)].y=sin( -(V[(j*xDim + i)] + delta[j*xDim + i])*(dt/(2*HBAR)));
		}
	}
}


int main(int argc, char **argv){
	
	time_t start,fin;
	time(&start);
	printf("Start: %s\n", ctime(&start));
	initArr(&params,32);
	//appendData(&params,ctime(&start),0.0);
	parseArgs(argc,argv);
	cudaSetDevice(device);
	//************************************************************//
	/*
	* Initialise the Params data structure to track params and variables
	*/
	//************************************************************//
	//paramS = (Params *) malloc(sizeof(Params));
	//strcpy(paramS->data,"INIT");
	//paramS->next=NULL;

	initialise(omegaX,omegaY,atoms);
	timeTotal = 0.0;
	//************************************************************//
	/*
	* Groundstate finder section
	*/
	//************************************************************//
	FileIO::writeOutParam(buffer, params, "Params.dat");
	if(read_wfc == 1){
		printf("Loading wavefunction...");
		wfc=FileIO::readIn("wfc_load","wfci_load",xDim, yDim);
		printf("Wavefunction loaded.\n");
	}
	
	double2 ph;
	double x_0,y_0;
	x_0 = 0;//(0.5*xDim)*dx;
	y_0 = 0;//(0.5*yDim)*dy;
/*	for(int i=0; i < xDim; i++ ){
		for(int j=0; j < yDim; j++ ){
			ph.x = cos( fmod( 0*atan2( y[j] - y_0, x[i] - x_0 ), 2*PI) );
			ph.y = -sin( fmod( 0*atan2( y[j] - y_0, x[i] - x_0 ), 2*PI) );
			wfc[(i*yDim + j)] = Minions::complexMult( wfc[(i*yDim + j)], ph );
		}
	}
	printf("l=%e\n",l);
*/	if(gsteps > 0){
		err=cudaMemcpy(K_gpu, GK, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(V_gpu, GV, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(xPy_gpu, xPy, sizeof(double)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(yPx_gpu, yPx, sizeof(double)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		
		evolve(wfc_gpu, K_gpu, V_gpu, yPx_gpu, xPy_gpu, par_sum, xDim*yDim, gsteps, threads, 0, ang_mom, gpe, print, atoms, 0);
		cudaMemcpy(wfc, wfc_gpu, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyDeviceToHost);
	}

	free(GV); free(GK); free(xPy); free(yPx);

	//************************************************************//
	/*
	* Evolution
	*/
	//************************************************************//
	if(esteps > 0){
		err=cudaMemcpy(xPy_gpu, ExPy, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(yPx_gpu, EyPx, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(xPy_gpu, ExPy, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(yPx_gpu, EyPx, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(K_gpu, EK, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
		err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		if(err!=cudaSuccess)
			exit(1);
			
		//delta_define(x, y, (523.6667 - 512 + x0_shift)*dx, (512.6667 - 512  + y0_shift)*dy, V_opt);
		FileIO::writeOutDouble(buffer,"V_opt",V_opt,xDim*yDim,0);
		evolve(wfc_gpu, K_gpu, V_gpu, yPx_gpu, xPy_gpu, par_sum, xDim*yDim, esteps, threads, 1, ang_mom, gpe, print, atoms, 0);
	
	}
	free(EV); free(EK); free(ExPy); free(EyPx);
	free(x);free(y);
	cudaFree(wfc_gpu); cudaFree(K_gpu); cudaFree(V_gpu); cudaFree(yPx_gpu); cudaFree(xPy_gpu); cudaFree(par_sum);

	time(&fin);
	//appendData(&params,ctime(&fin),0.0);
	printf("Finish: %s\n", ctime(&fin));
	printf("Total time: %ld seconds\n ",(long)fin-start);
	//appendData(&params,"t_duration",fin-start);
	return 0;
}
