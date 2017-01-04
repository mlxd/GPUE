
#include "../include/minions.h"

namespace Minions{

    double psi2(double2 in){
        return in.x*in.x + in.y*in.y;
    }

    double maxValue(double* grid,int len){
        double max = grid[0];
        for (int i=1;i<len-1;++i){
            if(max<grid[i]){
                max=grid[i];
            }
        }
        return max;
    }

    double minValue(double* grid,int len){
        double min = grid[0];
        for (int i=1;i<len-1;++i){
            if(min>grid[i])
                min=grid[i];
        }
        return min;
    }

    double sumAvg(double* in, int len){
        double avg = 0.0;

        for (int i=0; i<len; ++i){
            avg += in[i];
        }
        return avg/len;
    }

    /**
     * Double precision fast inverse square-root. Useless, but necessary to have.
     */
    double fInvSqRt(double in){
        long long l;
        double in05, calc;

        in05 = in*0.5;
        calc=in;
        l = * (long long*) &calc;
        l = 0x5fe6eb50c7b537a9LL - (l >> 1);
        calc = *(double *) &l;
        calc = calc*( 1.5 - (in05*calc*calc) );

        return calc;
    }



     void coordSwap(struct Vtx::Vortex *vCoords, int src, int dest){
        struct Vtx::Vortex d = vCoords[dest];
        vCoords[dest] = vCoords[src];
        vCoords[src] = d;
    }

     double complexMag(double2 in){
        return sqrt(in.x*in.x + in.y*in.y);
    }

     double complexMag2(double2 in){
        return in.x*in.x + in.y*in.y;
    }

     double2 complexMult(double2 in1, double2 in2){
        double2 result;
        result.x = (in1.x*in2.x - in1.y*in2.y);
        result.y = (in1.x*in2.y + in1.y*in2.x);
        return result;
    }

     double2 complexScale(double2 comp, double scale){
        double2 result;
        result.x = comp.x*scale;
        result.y = comp.y*scale;
        return result;
    }

     double2 conj(double2 c){
        double2 result = c;
        result.y = -result.y;
        return result;
    }

     double2 complexDiv(double2 num, double2 den){
        double2 c = conj(den);
        return complexScale(complexMult(num,c),(1.0/complexMag2(den)));
    }

    void trans2x2(double *in, double *out){
        out[0] = in[0];
        out[1] = in[2];
        out[2] = in[1];
        out[3] = in[3];
    }

    void inv2x2(double *in, double *out){
        double det = 1.0/(in[0]*in[3] - in[1]*in[2]);
        out[0] = det*in[3];
        out[1] = -det*in[1];
        out[2] = -det*in[2];
        out[3] = det*in[0];
    }
}

    /*
    int qSort(int2 *vCoords, int *vCoordsP int index, int length){
        if(index < 2){
            return 0;
        }
        int2 pivot;
        int l = 0;
        int r = length - 1;
        while (l <= r){
            0;
        }
    }
*/
