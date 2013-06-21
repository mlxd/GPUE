/*
* split_op.cu - GPUE: Split Operator based GPU solver for Nonlinear 
* Schrodinger Equation, Copyright (C) 2012, Lee J. O'Riordan, Tadhg 
* Morgan, Neil Crowley. 

* This library is free software; you can redistribute it and/or modify 
* it under the terms of the GNU Lesser General Public License as 
* published by the Free Software Foundation; either version 2.1 of the 
* License, or (at your option) any later version. This library is 
* distributed in the hope that it will be useful, but WITHOUT ANY 
* WARRANTY; without even the implied warranty of MERCHANTABILITY or 
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
* License for more details. You should have received a copy of the GNU 
* Lesser General Public License along with this library; if not, write 
* to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
* Boston, MA 02111-1307 USA 
*/

#include "../include/split_op.h"
#include "../include/kernels.h"
#include "../include/constants.h"
#include "../include/fileIO.h"
#include "../include/tracker.h"

__const__ double FFT_RENORM_2D = 0.0;
__const__ double FFT_RENORM_1D = 0.0;

int verbose;
int device;
double gammaY;
double omega;

int isError(int result, char* c){
	if(result!=0){printf("Error has occurred for method %s with return type %d\n",c,result);
		exit(result);
	}
	return result;
}

int initialise(double omegaX, double omegaY, int N){
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	int xD=1,yD=1;//,zD=1;
	threads = 256;
	int b = xDim*yDim/threads;  //number of blocks in simulation
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
	else
		printf("Outside range of supported indexing");
	printf("Compute grid dimensions chosen as X=%d	Y=%d\n",xD,yD);
	
	grid.x=xD; 
	grid.y=yD; 
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	int i,j; //Used in for-loops for indexing
	
	int gSize = xDim*yDim;
	
	double xOffset, yOffset;
	xOffset=0.0;//5.0e-6;
	yOffset=0.0;//5.0e-6;
	
	mass = 1.4431607e-25; //Rb 87 mass, kg
	double sum = 0.0;
	
	double a0x = sqrt(HBAR/(2*mass*omegaX));
	double a0y = sqrt(HBAR/(2*mass*omegaY));
	
	double Rxy;
	Rxy = pow(15,0.2)*pow(N*4.67e-9*sqrt(mass*pow(omegaX*omegaY,0.5)/HBAR),0.2);
	
	double xMax = 4*Rxy*a0x;
	double yMax = 4*Rxy*a0y;

	omega = omega*omegaX;
	
	double pxMax, pyMax;
	pxMax = (PI/xMax)*(xDim>>1);
	pyMax = (PI/yMax)*(yDim>>1);
	
	//double dx, dy, dz;
	dx = xMax/(xDim>>1);
	dy = yMax/(yDim>>1);
	
	double dpx, dpy;
	dpx = PI/(xMax);
	dpy = PI/(yMax);

	printf("a0x=%e  a0y=%e \n dx=%e   dx=%e\n R_xy=%e\n",a0x,a0y,dx,dy,Rxy);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	double *x,*y;
	double *xpn,*ypn,*xp,*yp;
	x = (double *) malloc(sizeof(double) * xDim);
	y = (double *) malloc(sizeof(double) * yDim);
	
	xp = (double *) malloc(sizeof(double) * xDim);
	yp = (double *) malloc(sizeof(double) * yDim);
	xpn = (double *) malloc(sizeof(double) * xDim);
	ypn = (double *) malloc(sizeof(double) * yDim);

	/*
	 * Pos and Mom grids
	 */
	for(int i=0; i < (xDim/2); i++){
		x[i] = -xMax + i*dx;		
		x[i + (xDim/2)] = i*dx;
		
		y[i] = -yMax + i*dy;		
		y[i + (yDim/2)] = i*dy;
		
		xp[i] = i*dpx;
		xp[i + (xDim/2)] = -pxMax +i*dpx;
		
		yp[i] = i*dpy;
		yp[i + (yDim/2)] = -pyMax +i*dpy;
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	
	/* Initialise wavefunction, momentum and position operators on host */
	r = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	Phi = (double *) malloc(sizeof(double) * gSize);
	wfc = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	wfc_backup = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * (gSize/threads));
	K = (double *) malloc(sizeof(double) * gSize);
	V = (double *) malloc(sizeof(double) * gSize);
	GK = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	GV = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EK = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EV = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	xPy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	yPx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	GxPy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	GyPx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	ExPy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EyPx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	EappliedField = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
	
	/* Initialise wfc, EKp, and EVr buffers on GPU */
	cudaMalloc((void**) &wfc_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &K_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &V_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &xPy_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &yPx_gpu, sizeof(cufftDoubleComplex) * gSize);
	cudaMalloc((void**) &par_sum, sizeof(cufftDoubleComplex) * (gSize/threads));
	cudaMalloc((void**) &r_gpu, sizeof(double2) * (gSize));
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	double l;
	int app_field = 2; //1 Gaussian, else cosine
	for( i=0; i < xDim; i++ ){
		for( j=0; j < yDim; j++ ){
			r[(i*yDim + j)].x = x[i] - xOffset; 
			r[(i*yDim + j)].y = y[j] - yOffset;
		
			//Remember, you are going from -PI to +PI, not 0 to 2PI
			if(x[i]>=0){
				Phi[(i*yDim + j)] = atan((y[j] + dx/10)/(x[i]) ) - PI/2.0;
			}
			else
				Phi[(i*yDim + j)] = atan((y[j] + dx/10)/(x[i]) ) + PI/2.0;
			l = 2;

			Phi[(i*yDim + j)] = fmod(l*Phi[(i*yDim + j)],2*PI);
			wfc[(i*yDim + j)].x = exp(-( pow((x[j])/(Rxy*a0x),2) + pow((y[i])/(Rxy*a0y),2) ) )*cos(Phi[(i*yDim + j)]);
			wfc[(i*yDim + j)].y = -exp(-( pow((x[j])/(Rxy*a0x),2) + pow((y[i])/(Rxy*a0y),2) ) )*sin(Phi[(i*yDim + j)]);
				
			V[(i*yDim + j)] = 0.5*mass*( pow(omegaX*(x[j]+xOffset),2) + pow(omegaY*(y[i]+yOffset),2) );
			K[(i*yDim + j)] = (HBAR*HBAR/(2*mass))*(xp[j]*xp[j] + yp[i]*yp[i]);

			GV[(i*yDim + j)].x = exp( -V[(i*yDim + j)]*(dt/HBAR));
			GK[(i*yDim + j)].x = exp( -K[(i*yDim + j)]*(dt/(2*HBAR)));
			GV[(i*yDim + j)].y = 0.0;
			GK[(i*yDim + j)].y = 0.0;
			
			xPy[(i*yDim + j)].x = x[i]*yp[j];
			xPy[(i*yDim + j)].y = 0.0;
			yPx[(i*yDim + j)].x = y[j]*xp[i];
			yPx[(i*yDim + j)].y = 0.0;
			
			GxPy[(i*yDim + j)].x = exp( -omega*xPy[(i*yDim + j)].x*dt);
			GxPy[(i*yDim + j)].y = 0.0;
			GyPx[(i*yDim + j)].x = exp( omega*yPx[(i*yDim + j)].x*dt);
			GyPx[(i*yDim + j)].y = 0.0;
			
			EV[(i*yDim + j)].x=cos( -V[(i*yDim + j)]*(dt/HBAR));
			EV[(i*yDim + j)].y=sin( -V[(i*yDim + j)]*(dt/HBAR));
			EK[(i*yDim + j)].x=cos( -K[(i*yDim + j)]*(dt/HBAR));
			EK[(i*yDim + j)].y=sin( -K[(i*yDim + j)]*(dt/HBAR));
			
			ExPy[(i*yDim + j)].x=cos(-omega*xPy[(i*yDim + j)].x*dt);
			ExPy[(i*yDim + j)].y=sin(-omega*xPy[(i*yDim + j)].x*dt);
			EyPx[(i*yDim + j)].x=cos(omega*yPx[(i*yDim + j)].x*dt);
			EyPx[(i*yDim + j)].y=sin(omega*yPx[(i*yDim + j)].x*dt);
	
			if(app_field==1){
				EappliedField[(i*yDim + j)].x = cos(fmod(140*(-1)*cos(x[i]/(dx*150))*cos(x[i]/(dx*150)),2*PI));
				EappliedField[(i*yDim + j)].y = sin(fmod(140*(-1)*cos(x[i]/(dx*150))*cos(x[i]/(dx*150)),2*PI));
			}
			else if(app_field==2){
				EappliedField[(i*yDim + j)].x = cos(fmod(30.0*exp( -(pow((x[i]-1.75e-5)/(dx*10),2) + pow((y[j]+0e-5)/(dy*10),2)) ),2*PI));
				EappliedField[(i*yDim + j)].y = sin(fmod(30.0*exp( -(pow((x[i]-1.75e-5)/(dx*10),2) + pow((y[j]+0e-5)/(dy*10),2)) ),2*PI));
			}
			sum+=sqrt(wfc[(i*yDim + j)].x*wfc[(i*yDim + j)].x + wfc[(i*yDim + j)].y*wfc[(i*yDim + j)].y);
		}
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	writeOutDouble(buffer,"V",V,xDim*yDim,0);
	writeOutDouble(buffer,"K",K,xDim*yDim,0);
	writeOut(buffer,"xPy",xPy,xDim*yDim,0);
	writeOut(buffer,"yPx",yPx,xDim*yDim,0);
	writeOut(buffer,"ExPy",ExPy,xDim*yDim,0);
	writeOut(buffer,"EyPx",EyPx,xDim*yDim,0);
	writeOutDouble(buffer,"Phi",Phi,xDim*yDim,0);
	writeOut(buffer,"r",r,xDim*yDim,0);
	writeOutDouble(buffer,"x",x,xDim,0);
	writeOutDouble(buffer,"y",y,yDim,0);
	writeOutDouble(buffer,"px",xp,xDim,0);
	writeOutDouble(buffer,"py",yp,yDim,0);
	writeOut(buffer,"ypx",xPy,xDim*yDim,0);
	writeOut(buffer,"xpy",yPx,yDim*xDim,0);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

	sum=sqrt(sum*dx*dy);
	for (i = 0; i < xDim; i++){
		for (j = 0; j < yDim; j++){
			wfc[(i*yDim + j)].x = (wfc[(i*yDim + j)].x)/(sum);
			wfc[(i*yDim + j)].y = (wfc[(i*yDim + j)].y)/(sum);
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
	
	free(x);free(y);free(xp);free(yp);free(xpn);free(ypn);
	return 0;
}

int evolve( cufftDoubleComplex *gpuWfc, 
			cufftDoubleComplex *gpuMomentumOp,
			cufftDoubleComplex *gpuPositionOp,
			cufftDoubleComplex *gpu1dyPx,
			cufftDoubleComplex *gpu1dxPy,
			cufftDoubleComplex *gpuParSum,			 
			int gridSize, int numSteps, int threads, 
			int gstate, int lz, int nonlin, int printSteps, int N){

	//double sum;
	double renorm_factor_2d=1.0/pow(gridSize,0.5);
	double renorm_factor_1d=1.0/pow(xDim,0.5);

	clock_t begin, end;
	double time_spent;

	begin = clock();
	for(int i=0; i < numSteps; i++){
		printf("Step: %d\n",i);
		if(i % printSteps == 0){
			cudaMemcpy(wfc, wfc_gpu, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyDeviceToHost);
			if(gstate == 0){
				end = clock();
				time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
				printf("Time spent: %lf\n",time_spent);
				writeOut(buffer, "wfc_0",wfc, xDim*yDim, i);
			}
			else
				writeOut(buffer, "wfc",wfc, xDim*yDim, i);
		}
		/*
		 * FFT 1 - Forward: r -> p
		 */ 
		result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD);
		isError(result,"Z2Z 1");
		scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc); //Normalise
		cMult<<<grid,threads>>>(gpuMomentumOp,gpuWfc,gpuWfc); //EKp complex Mult
		
		/*
		 * FFT 2 - Inverse: p -> r
		 */ 
		result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE);
		isError(result,"Z2Z 2");
		scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc); //Normalise
		if(nonlin == 1){
			cMultDensity<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc,dt,mass,omegaZ,gstate,N*interaction);
		}
		else {
			cMult<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc);
		}
				
		/*
		 * FFT 3 - Forward: r -> p
		 */		
		result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD);
		scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc); //Normalise
		cMult<<<grid,threads>>>(gpuMomentumOp,gpuWfc,gpuWfc);
		
		/*
		 * FFT 4 - Forward: p -> r
		 */	
		result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE);
		scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc); //Normalise

		/**************************************************************/
		if(lz == 1){
			result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_xPy
			scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
			cMult<<<grid,threads>>>(gpu1dxPy,gpuWfc,gpuWfc); //wfc_xPy complex Mult
			result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
			scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
			
			result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD); //2D forward
			scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
			result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE); //1D inverse to wfc_yPx
			scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
			cMult<<<grid,threads>>>(gpu1dyPx,gpuWfc,gpuWfc); //wfc_xPy complex Mult
			result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); // wfc_PxPy
			scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
			result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE); //2D Inverse
			scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
		}
		/**************************************************************/
		if(0){
			//pinVortex<<<grid,threads>>>(gpuWfc,gpuWfc,r_gpu);
		}
		parSum(gpuWfc, wfc_backup, gpuParSum, xDim, yDim, threads);
	}
	return 0;
}

/*
 * Used to perform parallel summation on WFC and normalise
 */
void parSum(double2* gpuWfc, double2* wfc_backup, double2* gpuParSum, int xDim, int yDim, int threads){
		int grid_tmp = xDim*yDim;
		int block = grid_tmp/threads;
		int thread_tmp = threads;
		int pass = 0;
		while((double)grid_tmp/threads > 1.0){
			//printf("G_TMP=%d	B=%d	T=%d\n",grid_tmp, block, thread_tmp);
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

double energyCalc(double2* wfcGpu, double2* potGpu, double2* kinGpu ){
	/*dx
	dy
	dr = dx*dy;
	double E_K, E_V;
	*/
	
	return 0.0;
}

int parseArgs(int argc, char** argv){

	int opt;
	while ((opt = getopt (argc, argv, "d:x:y:w:g:e:t:n:p:r:o:l:s:i:")) != -1) {
		switch (opt)
		{
			case 'v':
				verbose = 1;
				break;
			case 'x':
				xDim = atoi(optarg);
				printf("Argument for x is given as %s\n",optarg);
				break;
			case 'y':
				yDim = atoi(optarg);
				printf("Argument for y is given as %s\n",optarg);
				break;
			case 'w':
				omega = atof(optarg);
				printf("Argument for OmegaRotate is given as %s\n",optarg);
				break;
			case 'g':
				gsteps = atof(optarg);
				printf("Argument for Groundsteps is given as %ld\n",gsteps);
				break;
			case 'e':
				esteps = atof(optarg);
				printf("Argument for EvSteps is given as %ld\n",esteps);
				break;
			case 't':
				dt = atof(optarg);
				printf("Argument for Timestep is given as %E\n",dt);
				break;
			case 'd':
				device = atoi(optarg);
				printf("Argument for device is given as %d\n",device);
				break;
			case 'n':
				atoms = atof(optarg);
				printf("Argument for atoms is given as %ld\n",atoms);
				break;
			case 'r':
				read_wfc  = atoi(optarg);
				printf("Argument for ReadIn is given as %d\n",read_wfc);
				break;
			case 'p':
				print = atof(optarg);
				printf("Argument for Printout is given as %d\n",print);
				break;
			case 'l':
				ang_mom = atoi(optarg);
				printf("Angular Momentum mode engaged: %d\n",ang_mom);
				break;
			case 's':
				gpe = atoi(optarg);
				printf("Non-linear mode engaged: %d\n",gpe);
				break;
			case 'o':
				omegaZ = atof(optarg);
				printf("Argument for OmegaZ is given as %E\n",omegaZ);
				break;
			case 'i':
				interaction = atof(optarg);
				printf("Argument for interaction scaling is %E",interaction);
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

int main(int argc, char **argv)
{

	double value2d = 1.0/pow(xDim*yDim,0.5);
	double value1d = 1.0/pow(xDim,0.5);
	cudaMemcpyToSymbol(FFT_RENORM_2D, &value2d, sizeof(double));
	cudaMemcpyToSymbol(FFT_RENORM_1D, &value1d, sizeof(double));
	time_t start,fin;
	gsteps = 200000;
	esteps = 0;
	time(&start);
	printf("Start: %s\n", ctime(&start));
	parseArgs(argc,argv);
	printf("Device =%d\n",device);
	cudaSetDevice(device);

	initialise(2*PI,2*PI,atoms);
	
	//************************************************************//
	/*
	* Groundstate finder section
	*/
	//************************************************************//
	
	if(read_wfc == 1){
		readIn("wfc_0","wfci_0",xDim, yDim, wfc);
	}
	
	if(gsteps > 0){
		cudaMemcpy(r_gpu, r, sizeof(double2)*xDim*yDim, cudaMemcpyHostToDevice);
		cudaMemcpy(K_gpu, GK, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		cudaMemcpy(V_gpu, GV, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		cudaMemcpy(xPy_gpu, GxPy, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		cudaMemcpy(yPx_gpu, GyPx, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		
		evolve(wfc_gpu, K_gpu, V_gpu, yPx_gpu, xPy_gpu, par_sum, xDim*yDim, gsteps, 256, 0, ang_mom, gpe, print, atoms);
	}
	
	cudaMemcpy(V_gpu, EappliedField, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
	cMult<<<grid,threads>>>(V_gpu,wfc_gpu,wfc_gpu);

	//************************************************************//
	/*
	* Evolution
	*/
	//************************************************************//
	if(esteps > 0){
		cudaMemcpy(xPy_gpu, ExPy, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		cudaMemcpy(yPx_gpu, EyPx, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		cudaMemcpy(K_gpu, EK, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
		cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*xDim*yDim, cudaMemcpyHostToDevice);
			
		evolve(wfc_gpu, K_gpu, V_gpu, yPx_gpu, xPy_gpu, par_sum, xDim*yDim, esteps, 256, 1, ang_mom, gpe, print, atoms);
	}
	time(&fin);
	printf("Finish: %s\n", ctime(&fin));
	printf("Total time: %ld seconds\n ",(long)fin-start);
	return 0;
}
