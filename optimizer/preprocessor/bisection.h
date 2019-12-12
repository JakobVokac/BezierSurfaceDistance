//
// Created by s3179222 on 12/9/19.
//

#ifndef HEARTVALVEMODEL_BISECTION_H
#define HEARTVALVEMODEL_BISECTION_H


#include "preprocessor.h"

class bisection : public preprocessor{
public:
    explicit bisection(int iter) : preprocessor(iter){};
    ~bisection() = default;
    OptState2D preprocess(surface &sur, const vec3d &P) override;
    OptState1D preprocess(curve &crv, const vec3d &P) override;
};


#endif //HEARTVALVEMODEL_BISECTION_H
