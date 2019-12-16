//
// Created by s3179222 on 12/16/19.
//

#include "compositeBicubicsrf.h"

void compositeBicubicsrf::subdivide(BezierSurface *tl, BezierSurface *tr, BezierSurface *bl, BezierSurface *br) {
    tl = topLeft;
    tr = topRight;
    bl = bottomLeft;
    br = bottomRight;
}

bool compositeBicubicsrf::closestPointInPatch(vec3d P) {
    return  topLeft->closestPointInPatch(P) ||
            topRight->closestPointInPatch(P) ||
            bottomLeft->closestPointInPatch(P) ||
            bottomRight->closestPointInPatch(P);
}

bool compositeBicubicsrf::hasValidControlNet() {
    return topLeft->hasValidControlNet() ||
           topRight->hasValidControlNet() ||
           bottomLeft->hasValidControlNet() ||
           bottomRight->hasValidControlNet();
}

vec3d compositeBicubicsrf::at(double u, double v) {
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

void compositeBicubicsrf::recursiveControlNetSubdivide(BezierSurface *sur) {
    if(!sur->hasValidControlNet()){
        auto *temp = new compositeBicubicsrf(sur);
        delete sur;
        sur = temp;
        sur->controlNetSubdivide();
    }
}

void compositeBicubicsrf::controlNetSubdivide() {
    if(!hasValidControlNet()){
        recursiveControlNetSubdivide(topLeft);
        recursiveControlNetSubdivide(topRight);
        recursiveControlNetSubdivide(bottomLeft);
        recursiveControlNetSubdivide(bottomRight);
    }
}

curve &compositeBicubicsrf::edgeU0() {
    return *u0;
}

curve &compositeBicubicsrf::edgeU1() {
    return *u1;
}

curve &compositeBicubicsrf::edgeV0() {
    return *v0;
}

curve &compositeBicubicsrf::edgeV1() {
    return *v1;
}

