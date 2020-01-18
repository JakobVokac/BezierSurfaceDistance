//
// Created by s3179222 on 12/9/19.
//
/*
 * This is the interface class for preprocessing techniques for the optimizer, such as quadratic interpolation or bisection.
 *
 * The class requires any preprocessing method to be implemented for both 2D (surface search) and 1D (edge search).
 */
#ifndef HEARTVALVEMODEL_PREPROCESSOR_H
#define HEARTVALVEMODEL_PREPROCESSOR_H

#include "../states.h"
#include "../../geometry/surface/surface.h"

class preprocessor {
protected:
    int iterations;
public:
    explicit preprocessor(int iter){
        this->iterations = iter;
    };
    virtual OptState2D preprocess(surface &sur, const vec3d &P) = 0;
    virtual OptState1D preprocess(curve &crv, const vec3d &P) = 0;

};


#endif //HEARTVALVEMODEL_PREPROCESSOR_H
