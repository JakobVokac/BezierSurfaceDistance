//
// Created by s3179222 on 12/9/19.
//

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
