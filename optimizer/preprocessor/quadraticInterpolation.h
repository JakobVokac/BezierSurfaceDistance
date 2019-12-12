//
// Created by s3179222 on 12/10/19.
//

#ifndef HEARTVALVEMODEL_QUADRATICINTERPOLATION_H
#define HEARTVALVEMODEL_QUADRATICINTERPOLATION_H


#include "preprocessor.h"

class quadraticInterpolation : public preprocessor{
public:
    explicit quadraticInterpolation(int iter) : preprocessor(iter){};
    ~quadraticInterpolation() = default;
    OptState2D preprocess(surface &sur, const vec3d &P) override;
    OptState1D preprocess(curve &crv, const vec3d &P) override;
};


#endif //HEARTVALVEMODEL_QUADRATICINTERPOLATION_H
