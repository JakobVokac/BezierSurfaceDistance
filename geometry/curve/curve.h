//
// Created by s3179222 on 12/8/19.
//
/*
 * This class is an interface for a curve in 3D space. The interface exists for future work, if any kind of curve
 * that is not a normal Bezier curve should be implemented (such as splines). The interface implements the distance functions
 * and leaves the actual curve functions to be implemented by concrete classes.
 *
 * It also half way implements an output function for testing and debugging.
 *
 * It also contains a virtual function to be implemented for checking whether a curve sits on a plane, which is useful
 * in certain optimization applications such as the Geometric search.
 */
#ifndef HEARTVALVEMODEL_CURVE_H
#define HEARTVALVEMODEL_CURVE_H


#include "../vector/vec3d.h"
#include <iostream>

class curve {
public:
    virtual vec3d f(double t) = 0;
    virtual vec3d df(double t) = 0;
    virtual vec3d ddf(double t) = 0;

    double sqDistTo(double t, vec3d A);
    double sqDistToDer1(double t, vec3d A);
    double sqDistToDer2(double t, vec3d A);

    double distTo(double t, vec3d A);
    double distToDer1(double t, vec3d A);
    double distToDer2(double t, vec3d A);

    virtual std::string asString() = 0;
    friend std::ostream & operator << (std::ostream &out, curve &c);

    virtual vec3d curvePlaneNormal() = 0;
};


#endif //HEARTVALVEMODEL_CURVE_H
