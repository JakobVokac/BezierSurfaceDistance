//
// Created by s3179222 on 12/8/19.
//
/*
 * This class implements the surface interface, which is used by the optimizer, such that it can use any type of surface.
 *
 * The interface imeplements all of the distance equations, leaving the actual surface function to be implemented by
 * concrete subclasses.
 *
 * Each class must also be able to provide all of its edges for edge searches.
 */
#ifndef HEARTVALVEMODEL_SURFACE_H
#define HEARTVALVEMODEL_SURFACE_H

#include "../curve/cubiccrv.h"

class surface {
public:
    virtual vec3d at(double u, double v) = 0;
    virtual vec3d atDerU(double u, double v) = 0;
    virtual vec3d atDerV(double u, double v) = 0;
    virtual vec3d atDerUU(double u, double v) = 0;
    virtual vec3d atDerVV(double u, double v) = 0;
    virtual vec3d atDerUV(double u, double v) = 0;

    double sqDistTo(double u, double v, vec3d A);
    double sqDistToDerU(double u, double v, vec3d A);
    double sqDistToDerV(double u, double v, vec3d A);
    double sqDistToDerUU(double u, double v, vec3d A);
    double sqDistToDerVV(double u, double v, vec3d A);
    double sqDistToDerUV(double u, double v, vec3d A);

    double distTo(double u, double v, vec3d A);
    double distToDerU(double u, double v, vec3d A);
    double distToDerV(double u, double v, vec3d A);
    double distToDerUU(double u, double v, vec3d A);
    double distToDerVV(double u, double v, vec3d A);
    double distToDerUV(double u, double v, vec3d A);

    virtual curve & edgeU0() = 0;
    virtual curve & edgeU1() = 0;
    virtual curve & edgeV0() = 0;
    virtual curve & edgeV1() = 0;
protected:
private:
};


#endif //HEARTVALVEMODEL_SURFACE_H
