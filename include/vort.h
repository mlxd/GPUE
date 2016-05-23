/*** vort.h - GPUE: Split Operator based GPU solver for Nonlinear
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
//##############################################################################
/**
 *  @file    vort.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Class for keeping track of vortices. Implementation TBC
 *
 *  @section DESCRIPTION
 *  Each vortex is treated as a struct of position, winding, and least-squares
 *  calculated positions. UID, intervortex separations, and lifetime can be
 *  tracked and maintained.
 */
 //##############################################################################

#ifndef GPUE_1_VORT_H
#define GPUE_1_VORT_H

#include <memory>
#include <list>
#include<cuda.h>
#include<cuda_runtime.h>

//@todo Everything, as nothing is implemented in this class; struct Vortex is used though.
namespace Vtx {

	/**
	* Maintains vortex index in grid and least-squares calculated values
	* with winding/direction of vortex rotation.
	*/
    struct Vortex {
	    int2 coords;
	    double2 coordsD;
	    int wind;
    };
};

// Adding block comments in the least intuitive way because of doxygen comments
//##############################################################################
//	/**
//	* Vortex list for storing and retrieving vortices and associated values
//	*/
//	class VtxList {
//		private:
//			std::list< std::shared_ptr<Vortex> > vortices;
//		public:
//			/**
//			* @brief	Returns a reference to the vortex list
//			* @ingroup	vtx
//			* @return	Reference to vortex list
//			*/
//			std::list< std::shared_ptr<Vortex> > &getVortices();
//
//			/**
//			* @brief	Returns a reference to the vortex list
//			* @ingroup	vtx
//			* @return	Reference to vortex list
//			*/
//			std::shared_ptr<Vortex> get_Uid();
//
//			/**
//			* @brief	Returns a shared_ptr to the vortex by an index
//			* @ingroup	vtx
//			* @return	shared_ptr<Vortex> Shared pointer to vortex by index
//			*/
//			std::shared_ptr<Vortex> get_Idx();
//
//			/**
//			* @brief	Returns an index based upon a given vortex UID
//			* @ingroup	vtx
//			* @return	unsigned int Vortex index for UID
//			*/
//			unsigned int getIdx_Uid();
//
//			/**
//			* @brief	Returns the largest UID given
//			* @ingroup	vtx
//			* @return	unsigned int Largest vortex UID
//			*/
//			unsigned int getMax_Uid();
//
//			/**
//			* @brief	Returns index of vortex with shortest coordinate distance from current vortex
//			* @ingroup	vtx
//			* @return	unsigned int Vortex index for shortest distance
//			*/
//			unsigned int getIdx_MinDist();
//
//			/**
//			* @brief	In-place swap of the UID for the two given vortices
//			* @ingroup	vtx
//			* @param	v1 Shared pointer of Vortex v1
//			* @param	v2 Shared pointer of Vortex v2
//			*/
//			void swapUid(std::shared_ptr<Vortex> v1, std::shared_ptr<Vortex> v2);
//
//			/**
//			* @brief	Turns vortex activation off. Useful if vortex no longer exists in condensate.
//			* @ingroup	vtx
//			*/
//			void vortOff();
//	};
//}
//##############################################################################
#endif //GPUE_1_VORT_H
