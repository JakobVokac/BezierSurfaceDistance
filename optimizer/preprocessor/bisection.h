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
