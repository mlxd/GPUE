/*** operators.cc - GPUE: Split Operator based GPU solver for Nonlinear 
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

#include "../include/operators.h"

// Function for simple 2d rotation with i and j as the interators
double rotation_K(Grid &par, Op &opr, int i, int j, int k){
    double *xp = par.dsval("xp");
    double *yp = par.dsval("yp");
    double mass = par.dval("mass");
    return (HBAR*HBAR/(2*mass))*(xp[i]*xp[i] + yp[j]*yp[j]);
}

// Function for simple 2d rotation with i and j as the interators
double rotation_gauge_K(Grid &par, Op &opr, int i, int j, int k){
    double *xp = par.dsval("xp");
    double *yp = par.dsval("yp");
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double omega = par.dval("omega");
    double omegaX = par.dval("omegaX");
    double omega_0 = omega * omegaX;
    double mass = par.dval("mass");
    double p1 = HBAR*HBAR*(xp[i]*xp[i] + yp[j]*yp[j]);
    double p2 = mass*mass*omega_0*omega_0*(x[i]*x[i] + y[j]*y[j]);
    double p3 = 2*HBAR*mass*omega_0*(xp[i]*y[j] - yp[j]*x[i]);

    return (1/(2*mass))*(p1 + p2 + p3) *0.5;
}

// Function for simple 2d harmonic V with i and j as the iterators
double harmonic_V(Grid &par, Op &opr, int i , int j, int k){
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double gammaY = par.dval("gammaY");
    double yOffset = 0.0;
    double xOffset = 0.0;
    double mass = par.dval("mass");
    double V_x = omegaX*(x[i]+xOffset) - opr.Ax_fn(par.Afn)(par, opr, i, j, k);
    double V_y = gammaY*omegaY*(y[j]+yOffset) - 
                     opr.Ay_fn(par.Afn)(par, opr, i, j, k);
    return 0.5*mass*( V_x * V_x + V_y * V_y);

}

// Function for simple 2d harmonic V with i and j as iterators, gauge
double harmonic_gauge_V(Grid &par, Op &opr, int i , int j, int k){
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double omega = par.dval("omega");
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double gammaY = par.dval("gammaY");
    double omega_0 = omega * omegaX;
    double omega_1 = omega * omegaY;
    double yOffset = 0.0;
    double xOffset = 0.0;
    double mass = par.dval("mass");
    double ox = omegaX - omega_0;
    double oy = omegaY - omega_1;
    double v1 = ox * (x[i]+xOffset) * ox * (x[i]+xOffset) 
                + gammaY*oy*(y[j]+yOffset) * gammaY*oy*(y[j]+yOffset);
    return 0.5 * mass * (v1 );
    //return 0.5*mass*( pow(omegaX*(x[i]+xOffset) - omega_0,2) + 
    //       pow(gammaY*omegaY*(y[j]+yOffset) - omega_1,2) );

}

// Functions for pAx, y, z for rotation along the z axis
double rotation_pAx(Grid &par, Op &opr, int i, int j, int k){
    double *xp = par.dsval("xp");
    return opr.Ax_fn(par.Afn)(par, opr, i, j, k) * xp[i];
}

double rotation_pAy(Grid &par, Op &opr, int i, int j, int k){
    double *yp = par.dsval("yp");
    return opr.Ay_fn(par.Afn)(par, opr, i, j, k) * yp[j];
}

double rotation_Ax(Grid &par, Op &opr, int i, int j, int k){
    double *y = par.dsval("y");
    double omega = par.dval("omega");
    double omegaX = par.dval("omegaX");
    return -y[j] * omega * omegaX;
}

double rotation_Ay(Grid &par, Op &opr, int i, int j, int k){
    double *x = par.dsval("x");
    double omega = par.dval("omega");
    double omegaY = par.dval("omegaY");
    return x[i] * omega * omegaY;
}

// Fuinctions for pAx, y, z for rotation along the z axis
double rotation_squared_pAx(Grid &par, Op &opr, int i, int j, int k){
    double *y = par.dsval("y");
    double *xp = par.dsval("xp");
    double omega = par.dval("omega");
    double omegaX = par.dval("omegaX");
    //double yMax = par.dval("yMax");
    double val = -y[j]*y[j]*y[j]*xp[i]*xp[i]*xp[i] * omega * omegaX;
    //double val = 1/(y[j]*xp[i]) * omega * omegaX;
    if (abs(val) > 10){
        return 10;
    }
    else{
        return val;
    }
}

double rotation_squared_pAy(Grid &par, Op &opr, int i, int j, int k){
    double *x = par.dsval("x");
    double *yp = par.dsval("yp");
    double omega = par.dval("omega");
    double omegaY = par.dval("omegaY");
    //double xMax = par.dval("xMax");
    double val = x[i]*x[i]*x[i]*yp[j]*yp[j]*yp[j] * omega * omegaY;
    //double val = 1/(yp[j]*x[i]) * omega * omegaY;
    if (abs(val) > 10){
        return 10;
    }
    else{
        return val;
    }
}
