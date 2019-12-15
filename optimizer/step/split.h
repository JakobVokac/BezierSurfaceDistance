//
// Created by s3179222 on 12/15/19.
//

#ifndef HEARTVALVEMODEL_SPLIT_H
#define HEARTVALVEMODEL_SPLIT_H


#include "step.h"

class split : public step{
    OptState2D doStep(surface &sur, vec3d &P, OptState2D &loc) override;
    OptState1D doStep(curve &crv, vec3d &P, OptState1D &loc) override;
};


#endif //HEARTVALVEMODEL_SPLIT_H
