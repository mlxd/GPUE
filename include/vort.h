//
// Created by Lee James O'Riordan on 12/08/15.
//

#ifndef GPUE_1_VORT_H
#define GPUE_1_VORT_H

#include <memory>
#include <list>
#include<cuda.h>
#include<cuda_runtime.h>

namespace Vtx {
    struct Vortex {
	    int2 coords;
	    double2 coordsD;
	    int wind;
    };

	class VtxList {
		private:
			std::list< std::shared_ptr<Vortex> > vortices;
		public:
			std::list< std::shared_ptr<Vortex> > &getVortices();
			std::shared_ptr<Vortex> get_Uid();
			std::shared_ptr<Vortex> get_Idx();
			unsigned int getIdx_Uid();
			unsigned int getMax_Uid();
			unsigned int getIdx_MinDist();

			void swapUid(std::shared_ptr<Vortex> v1, std::shared_ptr<Vortex> v2);
			void vortOff();
	};
}

#endif //GPUE_1_VORT_H