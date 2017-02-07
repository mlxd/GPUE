
#include "../include/evolution.h"

void evolve_2d(Wave &wave, Op &opr,
               cufftDoubleComplex *gpuParSum, int numSteps, Cuda &cupar,
               unsigned int gstate, Grid &par, 
               std::string buffer){

    // Re-establishing variables from parsed Grid class
    std::string data_dir = par.sval("data_dir");
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
    double gDenConst = par.dval("gDenConst");
    double DX = par.dval("DX");
    double mask_2d = par.dval("mask_2d");
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double *V = opr.dsval("V");
    double *V_opt = opr.dsval("V_opt");
    double *Phi = wave.dsval("Phi");
    double *gpu1dpAx = opr.dsval("pAx_gpu");
    double *gpu1dpAy = opr.dsval("pAy_gpu");
    double *Phi_gpu = wave.dsval("Phi_gpu");
    int kick_it = par.ival("kick_it");
    bool write_it = par.bval("write_it");
    bool graph = par.bval("graph");
    int N = par.ival("atoms");
    int printSteps = par.ival("printSteps");
    bool nonlin = par.bval("gpe");
    bool lz = par.bval("corotating");
    bool ramp = par.bval("ramp");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int gridSize = xDim * yDim;
    int kill_idx = par.ival("kill_idx");
    cufftDoubleComplex *EV = opr.cufftDoubleComplexval("EV");
    cufftDoubleComplex *wfc = wave.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *EV_opt = opr.cufftDoubleComplexval("EV_opt");
    cufftDoubleComplex *gpuWfc = wave.cufftDoubleComplexval("wfc_gpu");
    cufftDoubleComplex *K_gpu =
        opr.cufftDoubleComplexval("K_gpu");
    cufftDoubleComplex *V_gpu =
        opr.cufftDoubleComplexval("V_gpu");

    std::cout << x[0] << '\t' << EV[0].x << '\t' << wfc[0].x << '\t'
              << EV_opt[0].x << '\t' << '\n';

    // getting data from Cuda class
    cufftResult result = cupar.cufftResultval("result");
    cufftHandle plan_1d = cupar.cufftHandleval("plan_1d");
    cufftHandle plan_2d = cupar.cufftHandleval("plan_2d");
    cufftHandle plan_other2d = cupar.cufftHandleval("plan_other2d");

    dim3 threads = cupar.dim3val("threads");
    dim3 grid = cupar.dim3val("grid");

    // Because no two operations are created equally. 
    // Multiplication is faster than divisions.
    double renorm_factor_2d=1.0/pow(gridSize,0.5);
    double renorm_factor_1d=1.0/pow(xDim,0.5);

    // outputting a bunch of variables just to check thigs out...
    std::cout << omega << '\t' << angle_sweep << '\t' << gdt << '\t'
              << dt << '\t' << omegaX << '\t' << omegaY << '\t' 
              << mass << '\t' << dx << '\t' << dy << '\t' << interaction << '\t'
              << laser_power << '\t' << N << '\t' << xDim << '\t' 
              << yDim << '\n';


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

    // ** ############################################################## ** //
    // **         HERE BE DRAGONS OF THE MOST DANGEROUS KIND!            ** //
    // ** ############################################################## ** //

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

    //std::cout << "numSteps is: " << numSteps << '\n';
    // Iterating through all of the steps in either g or esteps.
    for(int i=0; i < numSteps; ++i){
        if ( ramp ){
            //Adjusts omega for the appropriate trap frequency.
            omega_0=omegaX*((omega-0.39)*((double)i/(double)(numSteps)) + 0.39);
        }

        // Print-out at pre-determined rate.
        // Vortex & wfc analysis performed here also.
        if(i % printSteps == 0) { 
            // If the unit_test flag is on, we need a special case
            printf("Step: %d    Omega: %lf\n", i, omega_0 / omegaX);
            cudaMemcpy(wfc, gpuWfc, sizeof(cufftDoubleComplex) * xDim * yDim, 
                       cudaMemcpyDeviceToHost);

            // Printing out time of iteration
            end = clock();
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
            printf("Time spent: %lf\n", time_spent);
            std::string fileName = "";
            printf("ramp=%d        gstate=%d    rg=%d        \n", 
                   ramp, gstate, ramp | (gstate << 1));
            switch (ramp | (gstate << 1)) {
                case 0: //Groundstate solver, constant Omega value.
                    std::cout << "we are in case 0" << '\n';
                    fileName = "wfc_0_const";
                    break;
                case 1: //Groundstate solver, ramped Omega value.
                    std::cout << "we are in state 1" << '\n';
                    fileName = "wfc_0_ramp";
                    break;
                case 2: //Real-time evolution, constant Omega value.
                    std::cout << "we are in case 2" << '\n';
                    fileName = "wfc_ev";
                    vortexLocation = (int *) calloc(xDim * yDim, sizeof(int));
                    num_vortices[0] = Tracker::findVortex(vortexLocation, wfc,
                                                          mask_2d, xDim, x, i);

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
                                    V_opt, x, y, par, opr);
                        //V = opr.dsval("V");
                        //V_opt = opr.dsval("V_opt");
                        //EV_opt = opr.cufftDoubleComplexval("EV_opt");
                        sepAvg = Tracker::vortSepAvg(vortCoords, central_vortex,
                                                     num_vortices[0]);
                        if (kick_it == 2) {
                            printf("Kicked it 1\n");
                            cudaMemcpy(V_gpu, EV_opt, 
                                       sizeof(cufftDoubleComplex) * xDim * yDim,
                                       cudaMemcpyHostToDevice);
                        }
                        FileIO::writeOutDouble(buffer, data_dir + "V_opt_1",
                                               V_opt, xDim * yDim, 0);
                        FileIO::writeOut(buffer, data_dir + "EV_opt_1", EV_opt, 
                                         xDim * yDim, 0);
                        par.store("Central_vort_x", 
                                  (double) central_vortex.coords.x);
                        par.store("Central_vort_y", 
                                  (double) central_vortex.coords.y);
                        par.store("Central_vort_winding", 
                                  (double) central_vortex.wind);
                        par.store("Num_vort", (double) num_vortices[0]);
                        //std::cout << "writing to file in conditional" << '\n';
                        FileIO::writeOutParam(buffer, par, 
                                              data_dir + "Params.dat");
                    }
                    else if (num_vortices[0] > num_vortices[1]) {
                        printf("Number of vortices increased from %d to %d\n", 
                               num_vortices[1], num_vortices[0]);
                        Tracker::vortPos(vortexLocation, vortCoords, xDim, wfc);
                        Tracker::lsFit(vortCoords, wfc, num_vortices[0], xDim);
                    }
                    // if num_vortices[1] < num_vortices[0] ... Fewer vortices
                    else {
                        Tracker::vortPos(vortexLocation, vortCoords, xDim, wfc);
                        Tracker::lsFit(vortCoords, wfc, num_vortices[0], xDim);
                        Tracker::vortArrange(vortCoords, vortCoordsP, 
                                             num_vortices[0]);
                    }

                    // The following will be modified and moved into a new 
                    // library that works closely with GPUE
                    if (graph) {

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
                            auto killIt=[&](int idx, int winding, 
                                            double delta_x) {
                                WFC::phaseWinding(Phi, 1, x, y, dx, dy,
                                    lattice.getVortexUid(idx)->
                                    getData().coordsD.x 
                                    +cos(angle_sweep + vort_angle)*delta_x,
                                    lattice.getVortexUid(idx)->
                                    getData().coordsD.y
                                    +sin(angle_sweep + vort_angle)*delta_x,
                                    xDim);
                                cudaMemcpy(Phi_gpu, Phi, 
                                           sizeof(double) * xDim * yDim, 
                                           cudaMemcpyHostToDevice);
                                cMultPhi <<<grid, threads>>> (gpuWfc, Phi_gpu, 
                                                              gpuWfc);
                            };
                            if (kill_idx > 0){
                                killIt(kill_idx, 1, DX);
                            }

                        }
                        lattice.createEdges(1.5 * 2e-5 / dx);
                        adjMat = (double *)calloc(lattice.getVortices().size() *
                                                  lattice.getVortices().size(),
                                                   sizeof(double));
                        lattice.genAdjMat(adjMat);
                        FileIO::writeOutAdjMat(buffer, data_dir + "graph", 
                                               adjMat, uids, 
                                               lattice.getVortices().size(), i);
                        free(adjMat);
                        free(uids);
                        lattice.getVortices().clear();
                        lattice.getEdges().clear();
                        //exit(0);
                    }

                    FileIO::writeOutVortex(buffer, data_dir + "vort_arr",
                                           vortCoords, num_vortices[0], i);
                    printf("Located %d vortices\n", num_vortices[0]);
                    printf("Sigma=%e\n", vortOLSigma);
                    free(vortexLocation);
                    num_vortices[1] = num_vortices[0];
                    memcpy(vortCoordsP, vortCoords, 
                           sizeof(int2) * num_vortices[0]);
                    //exit(1);
                    //std::cout << "finished case 2" << '\n';
                    break;

                case 3:
                    fileName = "wfc_ev_ramp";
                    break;
                default:
                    break;
            }

            //std::cout << "writing" << '\n';
            if (write_it) {
                FileIO::writeOut(buffer, data_dir + fileName, 
                                 wfc, xDim * yDim, i);
            }
            //std::cout << "written" << '\n';
            //printf("Energy[t@%d]=%E\n",i,energy_angmom(V_gpu, 
            //       K_gpu, dx, dy, gpuWfc,gstate));
        }

        // No longer writing out

        // ** ########################################################## ** //
        // **                     More F'n' Dragons!                     ** //
        // ** ########################################################## ** //

        // If not already kicked at this time step more than 6 times... kick it!
        if(i%((int)t_kick+1) == 0 && num_kick<=6 && gstate==1 && kick_it == 1 ){
            cudaMemcpy(V_gpu, EV_opt, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
            ++num_kick;
        }
        // ** ########################################################## ** //

        // U_r(dt/2)*wfc
        if(nonlin == 1){
            //std::cout << Dt << '\t' << mass << '\t' << omegaZ << '\t' 
            //          << gstate << '\t' << N*interaction << '\n';
            cMultDensity<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc,0.5*Dt,
                                           mass,gstate,N*interaction,gDenConst);
        }
        else {
            cMult<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc);
        }
                
        // U_p(dt)*fft2(wfc)
        result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_FORWARD);

        // Normalise
        scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
        cMult<<<grid,threads>>>(K_gpu,gpuWfc,gpuWfc);
        result = cufftExecZ2Z(plan_2d,gpuWfc,gpuWfc,CUFFT_INVERSE);

        // Normalise
        scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_2d,gpuWfc);
        
        // U_r(dt/2)*wfc
        if(nonlin == 1){
            cMultDensity<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc,Dt*0.5,
                                           mass,gstate,N*interaction,gDenConst);
        }
        else {
            cMult<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc);
        }

        // If first timestep and kick_it >= 1, kick.
        // Also kick if not kicked enough
        if( (i % (int)(t_kick+1) == 0 && num_kick<=6 && gstate==1) || 
            (kick_it >= 1 && i==0) ){
            cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
            printf("Got here: Cuda memcpy EV into GPU\n");
        }
        // Angular momentum pAy-pAx (if engaged)  //
        if(lz == 1){
            // Multiplying by ramping factor if necessary
            // Note: using scalarPow to do the scaling inside of the exp
            if (ramp ){
                scalarPow<<<grid,threads>>>((cufftDoubleComplex*) gpu1dpAy, 
                                            omega_0/(omega * omegaY),
                                            (cufftDoubleComplex*) gpu1dpAy);
                scalarPow<<<grid,threads>>>((cufftDoubleComplex*) gpu1dpAx, 
                                            omega_0/(omega * omegaX),
                                            (cufftDoubleComplex*) gpu1dpAx);
            }
            switch(i%2 | (gstate<<1)){
                case 0: //Groundstate solver, even step

                    // 1d forward / mult by Ay
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAy, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d, gpuWfc);


                    // 1D FFT to wfc_pAx
                    result = cufftExecZ2Z(plan_other2d,gpuWfc,gpuWfc,
                                          CUFFT_FORWARD);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAx, gpuWfc);
    
                    result = cufftExecZ2Z(plan_other2d,gpuWfc,gpuWfc,
                                          CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d, gpuWfc);
                    break; 
                
                case 1:    //Groundstate solver, odd step
                    // 1D FFT to wfc_pAx
                    result = cufftExecZ2Z(plan_other2d,gpuWfc,gpuWfc,
                                          CUFFT_FORWARD);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAx, gpuWfc);
    
                    result = cufftExecZ2Z(plan_other2d,gpuWfc,gpuWfc,
                                          CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d, gpuWfc);

                    // wfc_pAy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAy, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d, gpuWfc);
                    break; 
                
                case 2: //Real time evolution, even step
                    //std::cout << "RT solver even." << '\n';

                    // wfc_pAy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAy, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                
                    // 1D to wfc_pAx
                    result = cufftExecZ2Z(plan_other2d,gpuWfc,gpuWfc,
                                          CUFFT_FORWARD);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAx, gpuWfc);

                    // wfc_pAy
                    result = cufftExecZ2Z(plan_other2d,gpuWfc,gpuWfc,
                                          CUFFT_INVERSE); 
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);

                    break;
                
                case 3:    //Real time evolution, odd step
                    //std::cout << "RT solver odd." << '\n';

                    // 1D inverse to wfc_pAx
                    result = cufftExecZ2Z(plan_other2d,gpuWfc,gpuWfc,
                                          CUFFT_FORWARD); 
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAx, gpuWfc);

                    // wfc_pAy
                    result = cufftExecZ2Z(plan_other2d,gpuWfc,gpuWfc,
                                          CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);

                    // wfc_pAy
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAy, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    break;
            
            }
        }
    
        if(gstate==0){
            parSum(gpuWfc, gpuParSum, par, cupar);
        }
    }

    // std::cout << "finished evolution" << '\n';
    // Storing wavefunctions for later
    //std::cout << gpuWfc[0].x << '\t' << gpuWfc[0].y << '\n';
    wave.store("wfc", wfc);
    wave.store("wfc_gpu", gpuWfc);

/*
    par.store("omega", omega);
    par.store("angle_sweep", angle_sweep);
    par.store("gdt", gdt);
    par.store("dt", dt);
    par.store("omegaX", omegaX);
    par.store("omegaY", omegaY);
    par.store("omegaZ", omegaZ);
    par.store("mass", mass);
    par.store("dx", dx);
    par.store("dy", dy);
    par.store("interaction", interaction);
    par.store("laser_power", laser_power);
    par.store("x", x);
    par.store("y", y);
    opr.store("V", V);
    opr.store("V_opt", V_opt);
    wave.store("Phi", Phi);
    opr.store("pAx_gpu", gpu1dpAx);
    opr.store("pAy_gpu", gpu1dpAy);
    wave.store("Phi_gpu", Phi_gpu);
    opr.store("EV", EV);
    //opr.store("V_gpu", V_gpu);
    //opr.store("K_gpu", K_gpu);
    opr.store("EV_opt", EV_opt);

    // getting data from Cuda class
    cupar.store("result", result);
    cupar.store("plan_1d", plan_1d);
    cupar.store("plan_2d", plan_2d);
    cupar.store("grid", grid);
*/

}

/*----------------------------------------------------------------------------//
* 3D
* Notes: In this case, we need to think about how to do the vortex tracking
*        Kicking will also be hard to do... Though not impossible, I suppose.
*-----------------------------------------------------------------------------*/

void evolve_3d(Wave &wave, Op &opr,
               cufftDoubleComplex *gpuParSum, int numSteps, Cuda &cupar,
               unsigned int gstate, Grid &par, 
               std::string buffer){

    // Re-establishing variables from parsed Grid class
    std::string data_dir = par.sval("data_dir");
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
    double dz = par.dval("dz");
    double interaction = par.dval("interaction");
    double laser_power = par.dval("laser_power");
    double gDenConst = par.dval("gDenConst");
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double *z = par.dsval("z");
    double *V = opr.dsval("V");
    double *V_opt = opr.dsval("V_opt");
    double *Phi = wave.dsval("Phi");
    double *gpu1dpAx = opr.dsval("pAx_gpu");
    double *gpu1dpAy = opr.dsval("pAy_gpu");
    double *gpu1dpAz = opr.dsval("pAz_gpu");
    double *Phi_gpu = wave.dsval("Phi_gpu");
    bool write_it = par.bval("write_it");
    bool graph = par.bval("graph");
    int N = par.ival("atoms");
    int printSteps = par.ival("printSteps");
    bool nonlin = par.bval("gpe");
    bool lz = par.bval("corotating");
    std::cout << "COROTATING IS: " << lz << '\n';
    bool ramp = par.bval("ramp");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int gridSize = xDim * yDim * zDim;
    cufftDoubleComplex *EV = opr.cufftDoubleComplexval("EV");
    cufftDoubleComplex *wfc = wave.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *EV_opt = opr.cufftDoubleComplexval("EV_opt");
    cufftDoubleComplex *gpuWfc = wave.cufftDoubleComplexval("wfc_gpu");
    cufftDoubleComplex *K_gpu =
        opr.cufftDoubleComplexval("K_gpu");
    cufftDoubleComplex *V_gpu =
        opr.cufftDoubleComplexval("V_gpu");

    std::cout << x[0] << '\t' << EV[0].x << '\t' << wfc[0].x << '\t'
              << EV_opt[0].x << '\t' << '\n';

    // getting data from Cuda class
    cufftResult result = cupar.cufftResultval("result");
    cufftHandle plan_1d = cupar.cufftHandleval("plan_1d");
    cufftHandle plan_3d = cupar.cufftHandleval("plan_3d");
    cufftHandle plan_dim2 = cupar.cufftHandleval("plan_dim2");
    cufftHandle plan_dim3 = cupar.cufftHandleval("plan_dim3");
    dim3 threads = cupar.dim3val("threads");
    dim3 grid = cupar.dim3val("grid");

    // Because no two operations are created equally. 
    // Multiplication is faster than divisions.
    double renorm_factor_3d=1.0/pow(gridSize,0.5);
    double renorm_factor_1d=1.0/pow(xDim,0.5);

    // outputting a bunch of variables just to check thigs out...
    std::cout << omega << '\t' << angle_sweep << '\t' << gdt << '\t'
              << dt << '\t' << omegaX << '\t' << omegaY << '\t' 
              << mass << '\t' << dx << '\t' << dy << '\t' << interaction << '\t'
              << laser_power << '\t' << N << '\t' << xDim << '\t' 
              << yDim << '\n';


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

    // ** ############################################################## ** //
    // **         HERE BE DRAGONS OF THE MOST DANGEROUS KIND!            ** //
    // ** ############################################################## ** //

    //std::cout << "numSteps is: " << numSteps << '\n';
    // Iterating through all of the steps in either g or esteps.
    for(int i=0; i < numSteps; ++i){
        if (ramp){
            //Adjusts omega for the appropriate trap frequency.
            omega_0=omegaX*((omega-0.39)*((double)i/(double)(numSteps)) + 0.39);
        }

        // Print-out at pre-determined rate.
        // Vortex & wfc analysis performed here also.
        if(i % printSteps == 0) { 
            // If the unit_test flag is on, we need a special case
            printf("Step: %d    Omega: %lf\n", i, omega_0 / omegaX);
            cudaMemcpy(wfc, gpuWfc, sizeof(cufftDoubleComplex)*xDim*yDim*zDim, 
                       cudaMemcpyDeviceToHost);

            // Printing out time of iteration
            end = clock();
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
            printf("Time spent: %lf\n", time_spent);
            std::string fileName = "";
            printf("ramp=%d        gstate=%d    rg=%d        \n", 
                   ramp, gstate, ramp | (gstate << 1));
            switch (ramp | (gstate << 1)) {
                case 0: //Groundstate solver, constant Omega value.
                    std::cout << "we are in case 0" << '\n';
                    fileName = "wfc_0_const";
                    break;
                case 1: //Groundstate solver, ramped Omega value.
                    std::cout << "we are in state 1" << '\n';
                    fileName = "wfc_0_ramp";
                    break;
                case 2: //Real-time evolution, constant Omega value.
                    // Note: In the case of 3d, we need to think about
                    //       vortex tracking in a new way.
                    //       It may be as simple as splitting the problem into
                    //       2D elements and working from there, but let's 
                    //       look into it when we need it in the future.
                    std::cout << "we are in case 2" << '\n';
                    fileName = "wfc_ev";
                    break;

                case 3:
                    fileName = "wfc_ev_ramp";
                    break;
                default:
                    break;
            }

            //std::cout << "writing" << '\n';
            if (write_it) {
                FileIO::writeOut(buffer, data_dir + fileName, 
                                 wfc, xDim*yDim*zDim, i);
            }
            //std::cout << "written" << '\n';
            //printf("Energy[t@%d]=%E\n",i,energy_angmom(V_gpu, 
            //       K_gpu, dx, dy, gpuWfc,gstate));
        }

        // No longer writing out

        // ** ########################################################## ** //
        // **                     More F'n' Dragons!                     ** //
        // ** ########################################################## ** //

        // U_r(dt/2)*wfc
        if(nonlin == 1){
            //std::cout << Dt << '\t' << mass << '\t' << omegaZ << '\t' 
            //          << gstate << '\t' << N*interaction << '\n';
            cMultDensity<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc,0.5*Dt,
                                           mass,gstate,N*interaction,gDenConst);
        }
        else {
            cMult<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc);
        }
                
        // U_p(dt)*fft2(wfc)
        result = cufftExecZ2Z(plan_3d,gpuWfc,gpuWfc,CUFFT_FORWARD);

        // Normalise
        scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_3d,gpuWfc);
        cMult<<<grid,threads>>>(K_gpu,gpuWfc,gpuWfc);
        result = cufftExecZ2Z(plan_3d,gpuWfc,gpuWfc,CUFFT_INVERSE);

        // Normalise
        scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_3d,gpuWfc);
        
        // U_r(dt/2)*wfc
        if(nonlin == 1){
            cMultDensity<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc,Dt*0.5,
                                           mass,gstate,N*interaction,gDenConst);
        }
        else {
            cMult<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc);
        }

        // Angular momentum pAy-pAx (if engaged)  //
        if(lz == true){
            // Multiplying by ramping factor if necessary
            // Note: using scalarPow to do the scaling inside of the exp
            if (ramp){
                scalarPow<<<grid,threads>>>((cufftDoubleComplex*) gpu1dpAy, 
                                            omega_0/(omega * omegaY),
                                            (cufftDoubleComplex*) gpu1dpAy);
                scalarPow<<<grid,threads>>>((cufftDoubleComplex*) gpu1dpAx, 
                                            omega_0/(omega * omegaX),
                                            (cufftDoubleComplex*) gpu1dpAx);
                scalarPow<<<grid,threads>>>((cufftDoubleComplex*) gpu1dpAz, 
                                            omega_0/(omega * omegaZ),
                                            (cufftDoubleComplex*) gpu1dpAz);
            }
            int size = xDim*zDim;
            switch(gstate){
                case 0: //Groundstate solver, even step

                    // 1d forward / mult by Az
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAz, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d, gpuWfc);

                    // loop to multiply by Ay
                    for (int i = 0; i < yDim; i++){
                        //size = xDim * zDim;
                        result = cufftExecZ2Z(plan_dim2,
                                 &gpuWfc[i*size],
                                 &gpuWfc[i*size],CUFFT_FORWARD);
                    }

                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAy, gpuWfc);

                    for (int i = 0; i < yDim; i++){
                        //size = xDim * zDim;
                        result = cufftExecZ2Z(plan_dim2,
                                 &gpuWfc[i*size],
                                 &gpuWfc[i*size],CUFFT_INVERSE);
                    }
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);

                    // 1D FFT to wfc_pAx
                    result = cufftExecZ2Z(plan_dim3,gpuWfc,gpuWfc,
                                          CUFFT_FORWARD);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAx, gpuWfc);
    
                    result = cufftExecZ2Z(plan_dim3,gpuWfc,gpuWfc,
                                          CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d, gpuWfc);

                    break; 
                
                case 1: //Real time evolution, even step
                    // 1d forward / mult by Az
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD); 
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAz, gpuWfc);
                    result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d, gpuWfc);

                    // loop to multiply by Ay
                    for (int i = 0; i < yDim; i++){
                        result = cufftExecZ2Z(plan_dim2,
                                 &gpuWfc[i*xDim*zDim],
                                 &gpuWfc[i*xDim*zDim],CUFFT_FORWARD);
                    }

                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAy, gpuWfc);

                    for (int i = 0; i < yDim; i++){
                        result = cufftExecZ2Z(plan_dim2,
                                 &gpuWfc[i*xDim*zDim],
                                 &gpuWfc[i*xDim*zDim],CUFFT_INVERSE);
                    }
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);


                    // 1D FFT to wfc_pAx
                    result = cufftExecZ2Z(plan_dim3,gpuWfc,gpuWfc,
                                          CUFFT_FORWARD);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d,gpuWfc);
                    cMult<<<grid,threads>>>(gpuWfc, 
                        (cufftDoubleComplex*) gpu1dpAx, gpuWfc);
    
                    result = cufftExecZ2Z(plan_dim3,gpuWfc,gpuWfc,
                                          CUFFT_INVERSE);
                    scalarMult<<<grid,threads>>>(gpuWfc,
                                                 renorm_factor_1d, gpuWfc);

                    break;
                
            }
        }
    
        if(gstate==0){
            parSum(gpuWfc, gpuParSum, par, cupar);
        }
    }

    // std::cout << "finished evolution" << '\n';
    // Storing wavefunctions for later
    //std::cout << gpuWfc[0].x << '\t' << gpuWfc[0].y << '\n';
    wave.store("wfc", wfc);
    wave.store("wfc_gpu", gpuWfc);
/*

    par.store("omega", omega);
    par.store("angle_sweep", angle_sweep);
    par.store("gdt", gdt);
    par.store("dt", dt);
    par.store("omegaX", omegaX);
    par.store("omegaY", omegaY);
    par.store("omegaZ", omegaZ);
    par.store("mass", mass);
    par.store("dx", dx);
    par.store("dy", dy);
    par.store("interaction", interaction);
    par.store("laser_power", laser_power);
    par.store("x", x);
    par.store("y", y);
    opr.store("V", V);
    opr.store("V_opt", V_opt);
    wave.store("Phi", Phi);
    opr.store("pAx_gpu", gpu1dpAx);
    opr.store("pAy_gpu", gpu1dpAy);
    wave.store("Phi_gpu", Phi_gpu);
    opr.store("EV", EV);
    //opr.store("V_gpu", V_gpu);
    //opr.store("K_gpu", K_gpu);
    opr.store("EV_opt", EV_opt);

    // getting data from Cuda class
    cupar.store("result", result);
    cupar.store("plan_1d", plan_1d);
    cupar.store("plan_2d", plan_2d);
    cupar.store("grid", grid);

*/
}
