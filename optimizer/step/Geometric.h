//
// Created by s3179222 on 12/10/19.
//

#ifndef HEARTVALVEMODEL_GEOMETRIC_H
#define HEARTVALVEMODEL_GEOMETRIC_H

#include "step.h"

class Geometric : public step {
    OptState2D doStep(surface &sur, vec3d &P, OptState2D &loc) override;
    OptState1D doStep(curve &crv, vec3d &P, OptState1D &loc) override;
};


#endif //HEARTVALVEMODEL_GEOMETRIC_H
