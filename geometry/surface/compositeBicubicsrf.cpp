//
// Created by s3179222 on 12/16/19.
//
/*
 * This class is an implementation of a Bezier surface, split into multiple patches. The main function of this class
 * is to take in a Bezier surface and continually split it until all parts have valid control nets (see source below).
 * Once the surface has been initialized, it can be used to check whether the whole Bezier surface has a closest point
 * to a point P on the inside or on an edge.
 *
 * The net and edge detection algorithms were taken and implemented from:
 *      Yingliang Ma & W. Terry Hewitt,
 *      "Point inversion and projection for NURBS curve: Control polygon approach",
 *      https://www.researchgate.net/publication/232638188_Point_inversion_and_projection_for_NURBS_curve_Control_polygon_approach
 */
#include "compositeBicubicsrf.h"


bool compositeBicubicsrf::closestPointInPatch(vec3d P) {
    if(hasValidControlNet()) {

        if(topLeft == nullptr)
            return srf.closestPointInPatch(P);

        return topLeft->closestPointInPatch(P) ||
               topRight->closestPointInPatch(P) ||
               bottomLeft->closestPointInPatch(P) ||
               bottomRight->closestPointInPatch(P);
    }
    return false;
}

bool compositeBicubicsrf::hasValidControlNet() {
    if(topLeft == nullptr)
        return srf.hasValidControlNet();

    return topLeft->hasValidControlNet() ||
           topRight->hasValidControlNet() ||
           bottomLeft->hasValidControlNet() ||
           bottomRight->hasValidControlNet();
}

vec3d compositeBicubicsrf::at(double u, double v) {
    if(topLeft == nullptr)
        return srf.at(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->at(u*2,v*2);
        }else{
            return topLeft->at(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->at((u-0.5)*2,v*2);
        }else{
            return topRight->at((u-0.5)*2,(v-0.5)*2);
        }
    }
}

vec3d compositeBicubicsrf::atDerU(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerU(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerU(u*2,v*2);
        }else{
            return topLeft->atDerU(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerU((u-0.5)*2,v*2);
        }else{
            return topRight->atDerU((u-0.5)*2,(v-0.5)*2);
        }
    }
}

vec3d compositeBicubicsrf::atDerV(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerV(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerV(u*2,v*2);
        }else{
            return topLeft->atDerV(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerV((u-0.5)*2,v*2);
        }else{
            return topRight->atDerV((u-0.5)*2,(v-0.5)*2);
        }
    }
}

vec3d compositeBicubicsrf::atDerUU(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerUU(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerUU(u*2,v*2);
        }else{
            return topLeft->atDerUU(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerUU((u-0.5)*2,v*2);
        }else{
            return topRight->atDerUU((u-0.5)*2,(v-0.5)*2);
        }
    }
}

vec3d compositeBicubicsrf::atDerVV(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerVV(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerVV(u*2,v*2);
        }else{
            return topLeft->atDerVV(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerVV((u-0.5)*2,v*2);
        }else{
            return topRight->atDerVV((u-0.5)*2,(v-0.5)*2);
        }
    }
}

vec3d compositeBicubicsrf::atDerUV(double u, double v) {
    if(topLeft == nullptr)
        return srf.atDerUV(u,v);
    if(u <= 0.5){
        if(v <= 0.5){
            return bottomLeft->atDerUV(u*2,v*2);
        }else{
            return topLeft->atDerUV(u*2,(v-0.5)*2);
        }
    }else{
        if(v <= 0.5){
            return bottomRight->atDerUV((u-0.5)*2,v*2);
        }else{
            return topRight->atDerUV((u-0.5)*2,(v-0.5)*2);
        }
    }
}



curve &compositeBicubicsrf::edgeU0() {
    return srf.edgeU0();
}

curve &compositeBicubicsrf::edgeU1() {
    return srf.edgeU1();
}

curve &compositeBicubicsrf::edgeV0() {
    return srf.edgeV0();
}

curve &compositeBicubicsrf::edgeV1() {
    return srf.edgeV1();
}

