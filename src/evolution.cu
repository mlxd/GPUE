/*
* evolution.cu - GPUE: Split Operator based GPU solver for Nonlinear 
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

#include "../include/evolution.h"
#include "../include/split_op.h"
#include "../include/kernels.h"
#include "../include/constants.h"
#include "../include/fileIO.h"
#include "../include/lattice.h"
#include "../include/manip.h"
#include <string>

void evolve( cufftDoubleComplex *gpuWfc, cufftDoubleComplex *gpuMomentumOp,
            cufftDoubleComplex *gpuPositionOp, void *gpu1dyPx, void *gpu1dxPy,
            cufftDoubleComplex *gpuParSum, int numSteps, Cuda cupar,
            unsigned int gstate, unsigned int ramp, Grid &par, char *buffer){

    // Re-establishing variables from parsed Grid class
    double omega = par.dval("omega");
    double angle_sweep = par.dval("angle_sweep");
    double gdt = par.dval("gdt");
    double dt = par.dval("dt");
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double omegaZ = par.dval("omegaZ");
    double mass = par.dval("mass");
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double interaction = par.dval("interaction");
    double laser_power = par.dval("laser_power");
    int kick_it = par.ival("kick_it");
    int write_it = par.ival("write_it");
    int graph = par.ival("graph");
    int N = par.ival("atoms");
    int printSteps = par.ival("print");
    int nonlin = par.ival("gpe");
    int lz = par.ival("ang_mom");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int gridSize = xDim * yDim;

    // getting data from Cuda class
    cufftResult result = cupar.cufftResultval("result");
    cufftHandle plan_1d = cupar.cufftHandleval("plan_1d");
    cufftHandle plan_2d = cupar.cufftHandleval("plan_2d");
    int threads = par.ival("threads");
    dim3 grid = cupar.dim3val("grid");
/*
    cufftResult result = cupar.result;
    cufftHandle plan_1d = cupar.plan_1d;
    cufftHandle plan_2d = cupar.plan_2d;
    int threads = cupar.threads;
*/

    // Because no two operations are created equally. 
    // Multiplication is faster than divisions.
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
    /** Determines the initial average density at the condensate central 
    * position and calculates a value for the healing length from this. Used 
    * thereafter as the lower limit for distances between vortices. **/
    int gridSum = 1<<6;
    double *densitySubset = (double*) malloc(sizeof(double)*gridSum);
    #pragma omp parallel for private(k)
    for (int j=0; j<gridSum; ++j){
        for (int k=0; k<gridSum; ++k){
            densitySubset[j*gridSum + k] = Minions::psi2(wfc[ ( (yDim/2) - 
                (gridSum/2) + j )*yDim  + ( (xDim/2)  - (gridSum/2) + k )]);
        }
    }
    // defined central condensate density
    xi = 1/sqrt(8*PI*a_s*Minions::sumAvg(densitySubset,gridSum)/(dx*dy));
    printf("Avg healing length at centre=%E\n",xi);
    #endif

    /** ** ############################################################## ** **/
    /** **         HERE BE DRAGONS OF THE MOST DANGEROUS KIND!            ** **/
    /** ** ############################################################## ** **/

    // Double buffering and will attempt to thread free and calloc operations to
    // hide time penalty. Or may not bother.
    int num_vortices[2] = {0,0};

    // binary matrix of size xDim*yDim, 
    // 1 for vortex at specified index, 0 otherwise
    int* vortexLocation;
    int* olMaxLocation = (int*) calloc(xDim*yDim,sizeof(int));

    struct Vtx::Vortex central_vortex; //vortex closest to the central position

    // Angle of vortex lattice. Add to optical lattice for alignment.
    double vort_angle;

    // array of vortex coordinates from vortexLocation 1's
    struct Vtx::Vortex *vortCoords = NULL;

    //Previous array of vortex coordinates from vortexLocation 1's
    struct Vtx::Vortex *vortCoordsP = NULL;

    LatticeGraph::Lattice lattice; //Vortex lattice graph.
    double* adjMat;
    
    double vortOLSigma=0.0;
    double sepAvg = 0.0;
    
    int num_kick = 0;
    double t_kick = (2*PI/omega_0)/(6*Dt);
    
    for(int i=0; i < numSteps; ++i){
        if ( ramp == 1 ){
            //Adjusts omega for the appropriate trap frequency.
            omega_0=omegaX*((omega-0.39)*((double)i/(double)(numSteps)) + 0.39);
        }

        // Print-out at pre-determined rate.
        // Vortex & wfc analysis performed here also.
        if(i % printSteps == 0) { 
            printf("Step: %d    Omega: %lf\n", i, omega_0 / omegaX);
            cudaMemcpy(wfc, gpuWfc, sizeof(cufftDoubleComplex) * xDim * yDim, 
                       cudaMemcpyDeviceToHost);
            end = clock();
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
            printf("Time spent: %lf\n", time_spent);
            std::string fileName = "";
            printf("ramp=%d        gstate=%d    rg=%d        \n", 
                   ramp, gstate, ramp | (gstate << 1));
            switch (ramp | (gstate << 1)) {
                case 0: //Groundstate solver, constant Omega value.
                    fileName = "wfc_0_const";
                    break;
                case 1: //Groundstate solver, ramped Omega value.
                    fileName = "wfc_0_ramp";
                    break;
                case 2: //Real-time evolution, constant Omega value.
                    fileName = "wfc_ev";
                    vortexLocation = (int *) calloc(xDim * yDim, sizeof(int));
                    num_vortices[0] = Tracker::findVortex(vortexLocation, wfc,
                                                         1e-4, xDim, par.x, i);

                    // If initial step, locate vortices, least-squares to find
                    // exact centre, calculate lattice angle, generate optical 
                    // lattice.
                    if (i == 0) {
                        vortCoords = (struct Vtx::Vortex *) malloc(
                                sizeof(struct Vtx::Vortex) * 
                                (2 * num_vortices[0]));
                        vortCoordsP = (struct Vtx::Vortex *) malloc(
                                sizeof(struct Vtx::Vortex) * 
                                (2 * num_vortices[0]));
                        Tracker::vortPos(vortexLocation, vortCoords, xDim, wfc);
                        Tracker::lsFit(vortCoords, wfc, num_vortices[0], xDim);
                        central_vortex = Tracker::vortCentre(vortCoords, 
                                                             num_vortices[0], 
                                                             xDim);
                        vort_angle = Tracker::vortAngle(vortCoords, 
                                                        central_vortex, 
                                                        num_vortices[0]);
                        par.store("Vort_angle", vort_angle);
                        optLatSetup(central_vortex, V, vortCoords, 
                                    num_vortices[0], 
                                    vort_angle + PI * angle_sweep / 180.0,
                                    laser_power * HBAR * sqrt(omegaX * omegaY),
                                    V_opt, par.x, par.y, par);
                        sepAvg = Tracker::vortSepAvg(vortCoords, central_vortex,
                                                     num_vortices[0]);
                        if (kick_it == 2) {
                            printf("Kicked it 1\n");
                            cudaMemcpy(V_gpu, EV_opt, 
                                       sizeof(cufftDoubleComplex) * xDim * yDim,
                                       cudaMemcpyHostToDevice);
                        }
                        FileIO::writeOutDouble(buffer, "V_opt_1", V_opt, 
                                               xDim * yDim, 0);
                        FileIO::writeOut(buffer, "EV_opt_1", EV_opt, 
                                         xDim * yDim, 0);
                        par.store("Central_vort_x", 
                                  (double) central_vortex.coords.x);
                        par.store("Central_vort_y", 
                                  (double) central_vortex.coords.y);
                        par.store("Central_vort_winding", 
                                  (double) central_vortex.wind);
                        par.store("Num_vort", (double) num_vortices[0]);
                        FileIO::writeOutParam(buffer, par, "Params.dat");
                    }
                    else if (num_vortices[0] > num_vortices[1]) {
                        printf("Number of vortices increased from %d to %d\n", 
                               num_vortices[1], num_vortices[0]);
                        Tracker::vortPos(vortexLocation, vortCoords, xDim, wfc);
                        Tracker::lsFit(vortCoords, wfc, num_vortices[0], xDim);
                    }
                    else {
                        Tracker::vortPos(vortexLocation, vortCoords, xDim, wfc);
                        Tracker::lsFit(vortCoords, wfc, num_vortices[0], xDim);
                        Tracker::vortArrange(vortCoords, vortCoordsP, 
                                             num_vortices[0]);
                    }

                    if (graph == 1) {

                        for (int ii = 0; ii < num_vortices[0]; ++ii) {
                            std::shared_ptr<LatticeGraph::Node> 
                                n(new LatticeGraph::Node(vortCoords[ii]));
                            lattice.addVortex(std::move(n));
                        }
                        unsigned int *uids = (unsigned int *) malloc(
                                sizeof(unsigned int) *
                                lattice.getVortices().size());
                        for (size_t a=0; a < lattice.getVortices().size(); ++a){
                            uids[a] = lattice.getVortexIdx(a)->getUid();
                        }
                        if(i==0) {
                            //Lambda for vortex annihilation/creation.
                            auto killIt=[&](int idx) {
                                WFC::phaseWinding(Phi, 1, par.x, par.y, dx, dy,
                                                  lattice.getVortexUid(idx)->
                                                      getData().coordsD.x,
                                                  lattice.getVortexUid(idx)->
                                                      getData().coordsD.y,xDim);
                                cudaMemcpy(Phi_gpu, Phi, 
                                           sizeof(double) * xDim * yDim, 
                                           cudaMemcpyHostToDevice);
                                cMultPhi <<<grid, threads>>> (gpuWfc, Phi_gpu, 
                                                              gpuWfc);
                            };
                            //killIt(44); //Kills vortex with UID 44


                        }
                        lattice.createEdges(1.5 * 2e-5 / dx);
                        adjMat = (double *)calloc(lattice.getVortices().size() *
                                                  lattice.getVortices().size(),
                                                   sizeof(double));
                        lattice.genAdjMat(adjMat);
                        FileIO::writeOutAdjMat(buffer, "graph", adjMat, uids, 
                                               lattice.getVortices().size(), i);
                        free(adjMat);
                        free(uids);
                        lattice.getVortices().clear();
                        lattice.getEdges().clear();
                        //exit(0);
                    }

                    FileIO::writeOutVortex(buffer, "vort_arr", vortCoords, 
                                           num_vortices[0], i);
                    printf("Located %d vortices\n", num_vortices[0]);
                    printf("Sigma=%e\n", vortOLSigma);
                    free(vortexLocation);
                    num_vortices[1] = num_vortices[0];
                    memcpy(vortCoordsP, vortCoords, 
                           sizeof(int2) * num_vortices[0]);
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
            // printf("Energy[t@%d]=%E\n",i,energy_angmom(gpuPositionOp, 
            //        gpuMomentumOp, dx, dy, gpuWfc,gstate));
/*
            cudaMemcpy(V_gpu, V, sizeof(double)*xDim*yDim, 
                         cudaMemcpyHostToDevice);
            cudaMemcpy(K_gpu, K, sizeof(double)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
            cudaMemcpy(V_gpu, , sizeof(double)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
            cudaMemcpy(K_gpu, K, sizeof(double)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
*/        
        }
    
    /** ** ############################################################## ** **/
    /** **                       More F'n' Dragons!                       ** **/
    /** ** ############################################################## ** **/
        if(i%((int)t_kick+1) == 0 && num_kick<=6 && gstate==1 && kick_it == 1 ){
            cudaMemcpy(V_gpu, EV_opt, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
            ++num_kick;
        }
    /** ** ############################################################## ** **/

        /*
         * U_r(dt/2)*wfc
         */ 
        if(nonlin == 1){
            cMultDensity<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc,0.5*Dt,
                                           mass,omegaZ,gstate,N*interaction);
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
            cMultDensity<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc,Dt*0.5,
                                           mass,omegaZ,gstate,N*interaction);
        }
        else {
            cMult<<<grid,threads>>>(gpuPositionOp,gpuWfc,gpuWfc);
        }
        if( (i % (int)(t_kick+1) == 0 && num_kick<=6 && gstate==1) || 
            (kick_it >= 1 && i==0) ){
            cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
            printf("Got here: Cuda memcpy EV into GPU\n");
        }
        /**************************************************************/
        /* Angular momentum xPy-yPx   */
        if(lz == 1){
            switch(i%2 | (gstate<<1)){
                case 0: //Groundstate solver, even step

                    // wfc_xPy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    angularOp<<<grid,threads>>>(omega_0, Dt, gpuWfc, 
                                                (double*) gpu1dxPy, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
    
                    // 2D forward
                    result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
    
                    // 1D inverse to wfc_yPx
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    angularOp<<<grid,threads>>>(omega_0, Dt, gpuWfc, 
                                                (double*) gpu1dyPx, gpuWfc);
    
                    // wfc_PxPy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
    
                    // 2D Inverse
                    result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
                    break;
                
                case 1:    //Groundstate solver, odd step

                    // 2D forward
                    result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);

                    // 1D inverse to wfc_yPx
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    angularOp<<<grid,threads>>>(omega_0, Dt, gpuWfc, 
                                                (double*) gpu1dyPx, gpuWfc);

                    // wfc_PxPy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);

                    // 2D inverse
                    result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
                    
                    // wfc_xPy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    angularOp<<<grid,threads>>>(omega_0, Dt, gpuWfc, 
                                                (double*) gpu1dxPy, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    break;
                
                case 2: //Real time evolution, even step

                    // wfc_xPy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dxPy, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                
                    //2D forward
                    result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);

                    // 1D inverses to wfc_yPx
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dyPx, gpuWfc);

                    // wfc_PxPy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);

                    // 2D Inverse
                    result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
                    break;
                
                case 3:    //Real time evolution, odd step

                    // 2D forward
                    result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);

                    // 1D inverse to wfc_yPx
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dyPx, gpuWfc);

                    // wfc_PxPy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);

                    // 2D inverse
                    result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
                    
                    // wfc_xPy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dxPy, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarDiv<<<grid,threads>>>(gpuWfc,renorm_factor_1d,gpuWfc);
                    break;
            
            }
        }
        /**************************************************************/
    
        if(gstate==0){
            parSum(gpuWfc, gpuParSum, threads, par, cupar);
        }
    }
}

