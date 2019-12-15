//
// Created by s3179222 on 12/8/19.
//

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

    virtual std::string asString() = 0;
    friend std::ostream & operator << (std::ostream &out, curve &c);

    virtual vec3d curvePlaneNormal() = 0;
};


#endif //HEARTVALVEMODEL_CURVE_H
