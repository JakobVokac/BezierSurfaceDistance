//
// Created by s3179222 on 12/10/19.
//
/*
 * This class is the implementation of the Geometric method for the numerical optimizer.
 *
 * Source: Shi-min Hu and Johannes Wallner,
 *        "A second order algorithm for orthogonal projection onto curves and surfaces",
 *        http://www.geometrie.tugraz.at/wallner/sproj.pdf
 */
#ifndef HEARTVALVEMODEL_GEOMETRIC_H
#define HEARTVALVEMODEL_GEOMETRIC_H

#include "step.h"

class Geometric : public step {
    OptState2D doStep(surface &sur, vec3d &P, OptState2D &loc) override;
    OptState1D doStep(curve &crv, vec3d &P, OptState1D &loc) override;
};


#endif //HEARTVALVEMODEL_GEOMETRIC_H
