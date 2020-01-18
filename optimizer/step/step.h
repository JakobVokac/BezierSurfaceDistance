//
// Created by s3179222 on 12/9/19.
//
/*
 * This is the interface class for the numerical optimization step, which is mean't to represent one step in an optimization
 * technique such as the Newton-Raphson method.
 *
 * The interface requires the method be implemented for both 2D (surface search) and 1D (edge search).
 */
#ifndef HEARTVALVEMODEL_STEP_H
#define HEARTVALVEMODEL_STEP_H


#include "../states.h"
#include "../../geometry/surface/surface.h"

class step {
public:
    virtual OptState2D doStep(surface &sur, vec3d &P, OptState2D &loc) = 0;
    virtual OptState1D doStep(curve &crv, vec3d &P, OptState1D &loc) = 0;
};


#endif //HEARTVALVEMODEL_STEP_H
