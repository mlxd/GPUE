/*
* fileIO.h - GPUE: Split Operator based GPU solver for Nonlinear 
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


#ifndef FILEIO_H
#define FILEIO_H

/* Buffer and FILE for IO */
char buffer[100];


int readIn(char*, char*, int, int, double2 *arr);
void writeOut(char*, char*, double2*, int, int );
void writeOutDouble(char*, char*, double*, int, int);
int readState(char*);

#endif
