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
double rotation_squared_Ax(Grid &par, Op &opr, int i, int j, int k){
    double *y = par.dsval("y");
    double omega = par.dval("omega");
    double omegaX = par.dval("omegaX");
    //double yMax = par.dval("yMax");
    double val = -y[j]*y[j] * omega * omegaX;
    return val;
}

double rotation_squared_Ay(Grid &par, Op &opr, int i, int j, int k){
    double *x = par.dsval("x");
    double omega = par.dval("omega");
    double omegaY = par.dval("omegaY");
    //double xMax = par.dval("xMax");
    double val = x[i]*x[i] * omega * omegaY;
    return val;
}

double dynamic_Ax(Grid &par, Op &opr, int i, int j, int k){
    double val = 0;
    std::string equation = par.sval("Axstring");
    parse_equation(par, equation, val, i, j, k);
    // For debugging
    //exit(0);
    return val;
}

double dynamic_Ay(Grid &par, Op &opr, int i, int j, int k){
    double val = 0;
    std::string equation = par.sval("Aystring");
    parse_equation(par, equation, val, i, j, k);
    // For debugging
    //exit(0);
    return val;
}

double dynamic_Az(Grid &par, Op &opr, int i, int j, int k){
    double val = 0;
    std::string equation = par.sval("Aystring");
    parse_equation(par, equation, val, i, j, k);
    return val;
}

// This function will be used with the dynamic gauge fields for AX,y,z (above)
void parse_equation(Grid par, std::string &equation, double &val, 
                    int i, int j, int k){

    // boolean value iff first minus
    bool minus = false;

    // Because this will be called recursively, we need to return if the string
    // length is 0
    if (equation.length() == 0){
        std::cout << "There's nothing here!" << '\n';
        //return;
    }

    // vector of all possibe mathematical operators (not including functions)
    std::vector<std::string> moperators(4);
    moperators = {
        "-", "/", "*", "+"
    };

    // And another vector for brackets of various types which indicate recursive
    // parsing of the equation
    std::vector<std::string> mbrackets;
    mbrackets = {
        "(", "[", "]", ")"
    };

    // vector of all possible mathematical functions... more to come
    std::vector<std::string> mfunctions(5);
    mfunctions = {
        "sin", "cos", "exp", "tan", "erf", "sqrt"
    };

    // We also need a specific map for the functions above
    typedef double (*functionPtr)(double);
    std::unordered_map<std::string, functionPtr> mfunctions_map;
    mfunctions_map["sin"] = sin;
    mfunctions_map["cos"] = cos;
    mfunctions_map["tan"] = tan;
    mfunctions_map["exp"] = exp;
    mfunctions_map["erf"] = erf;
    mfunctions_map["sqrt"] = sqrt;

    // We will have values and operators, but some operators will need to 
    // recursively call this function (think exp(), sin(), cos())...
    // We now need to do some silly sorting to figure out which operator 
    // comes first and where it is
    size_t index = equation.length();
    std::string currmop = "";
    size_t moppos;
    for (auto &mop : moperators){
        moppos = equation.find(mop);
        if (moppos < equation.length()){
            if (moppos < index && moppos > 0){
                currmop = mop;
                index = moppos;
            }
            else if(moppos == 0){
                minus = true;
                equation = equation.substr(1,equation.size());
            }
            else{
                currmop = equation.length();
            }
        }
    }

    //std::cout << currmop << '\t' << index << '\n';

    // Now we do a similar thing for the mbrackets
    // Sharing moppos from above
    for (auto &mbra : mbrackets){
        moppos = equation.find(mbra);
        if (moppos < equation.length()){
            if (moppos < index){
                currmop = mbra;
                index = moppos;
            }
            else{
                currmop = equation.length();
            }
        }
    }

    // Now we need to get the string we are working with...
    std::string item = equation.substr(0,index);

    // now we need to find the string in either mfunctions or par
    // First, we'll check mfunctions

    // Now we need to check to see if the string is in mfunctions
    auto it = mfunctions_map.find(item);
    if (it != mfunctions_map.end()){
        int openbracket, closebracket;
        openbracket = index;
        closebracket = equation.find(equation[openbracket]);
        std::string ineqn = equation.substr(openbracket + 1, 
                                            closebracket - 1);
        double inval = 0;
        parse_equation(par, ineqn, inval, i, j, k);
        val = mfunctions_map[item](inval);
    }

    // Now we need to do a similar thing for all the maps in par.
    if (par.is_double(item)){
        val = par.dval(item);
    }
    else if (par.is_dstar(item)){
        if (item == "x" || item == "px"){
            val = par.dsval(item)[i];
        }
        if (item == "y" || item == "py"){
            val = par.dsval(item)[j];
        }
        if (item == "z" || item == "pz"){
            val = par.dsval(item)[k];
        }
    }
    else if (item.size() > 0){
        std::cout << "could not find string " << item << "! please use one of "
                  << "the following variables:" << '\n';
        par.print_map();
    }

    if (minus){
        val *= -1;
    }

    //std::cout << item << '\t' << currmop << '\n';

    // Now to deal with the operator at the end
    if (currmop == "+"){
        double inval = 0;
        std::string new_eqn = equation.substr(index+1,equation.size());
        //std::cout << new_eqn << '\n';
        parse_equation(par, new_eqn, inval, i, j, k);
        val += inval;
    }
    if (currmop == "-"){
        double inval = 0;
        std::string new_eqn = equation.substr(index+1,equation.size());
        //std::cout << new_eqn << '\n';
        parse_equation(par, new_eqn, inval, i, j, k);
        val -= inval;
    }
    if (currmop == "*"){
        double inval = 0;
        std::string new_eqn = equation.substr(index+1,equation.size());
        //std::cout << new_eqn << '\n';
        parse_equation(par, new_eqn, inval, i, j, k);
        val *= inval;
    }
    if (currmop == "/"){
        double inval = 0;
        std::string new_eqn = equation.substr(index+1,equation.size());
        //std::cout << new_eqn << '\n';
        parse_equation(par, new_eqn, inval, i, j, k);
        val /= inval;
    }
}
