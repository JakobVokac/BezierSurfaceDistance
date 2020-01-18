//
// Created by s3179222 on 12/9/19.
//
/*
 * This is the binary-search/bisection implementation for preprocessing.
 *
 * The 1D version uses basic bisection.
 *
 * The 2D version uses a binary search, which works better for surfaces. It splits the surface into two parts per parameter,
 * takes the median of each side and continues with the closer side to point P. It alternates between u and v.
 *
 * Along further work, the binary search should be cut out of this class and placed in its own class and a proper 2D
 * bisection method should be implemented here.
 */
#include "bisection.h"

OptState2D bisection::preprocess(surface &sur, const vec3d &P) {
    double ul = 0, ur = 1, vl = 0, vr = 1;

    double um, vm;
    for (int i = 0; i < 2*iterations; ++i) {
        if( i % 2 == 0 ) {
            um = (ul+ur)/2;
            vm = (vl+vr)/2;
            double u1 = (ul+um)/2;
            double u2 = (ur+um)/2;
            double dist1 = sur.sqDistTo(u1,vm,P);
            double dist2 = sur.sqDistTo(u2,vm,P);
            if(dist1 < dist2){
                ur = um;
            }else{
                ul = um;
            }
        }else{
            um = (ul+ur)/2;
            vm = (vl+vr)/2;
            double v1 = (vl + vm)/2;
            double v2 = (vr + vm)/2;
            double dist1 = sur.sqDistTo(um,v1,P);
            double dist2 = sur.sqDistTo(um,v2,P);
            if(dist1 < dist2){
                vr = vm;
            }else{
                vl = vm;
            }
        }
    }
    return {um, vm, sur.distTo(um, vm, P)};
}

OptState1D bisection::preprocess(curve &crv, const vec3d &P) {
    double t1 = 0, t2 = 1, tm, hm, ht;
    for (int i = 0; i < 8; i++) {
        tm = (t1+t2)/2;
        hm = crv.sqDistToDer1(tm,P);
        ht = crv.sqDistToDer1(t1,P);
        if(hm == 0)
            return {tm,crv.distTo(tm,P)};
        if((hm < 0 && ht < 0) || (hm > 0 && ht > 0)){
            t1 = tm;
        }else{
            t2 = tm;
        }
    }
    return {tm,crv.distTo(tm,P)};
}
