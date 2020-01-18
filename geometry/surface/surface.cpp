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
#include "surface.h"


double surface::sqDistTo(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    double dist = top.sqdist(A);
    return dist;
}

double surface::sqDistToDerU(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer = atDerU(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}

double surface::sqDistToDerV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer = atDerV(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}

double surface::sqDistToDerUU(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer1 = atDerU(u,v);
    vec3d topDer2 = atDerUU(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}

double surface::sqDistToDerVV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDer1 = atDerV(u,v);
    vec3d topDer2 = atDerVV(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}

double surface::sqDistToDerUV(double u, double v, class vec3d A) {
    vec3d top = at(u,v);
    vec3d topDerU = atDerU(u,v);
    vec3d topDerV = atDerV(u,v);
    vec3d topDerUV = atDerUV(u,v);

    double dist = 2*(((top - A)*topDerUV).sum() + (topDerU*topDerV).sum());

    return dist;
}

double surface::distTo(double u, double v, class vec3d A) {
    return sqrt(sqDistTo(u,v,A));
}

double surface::distToDerU(double u, double v, vec3d A){

    return sqDistToDerU(u,v,A) / distTo(u,v,A);
}

double surface::distToDerV(double u, double v, vec3d A) {

    return sqDistToDerV(u, v, A) / distTo(u, v, A);
}

double surface::distToDerUU(double u, double v, vec3d A) {
    vec3d f = at(u,v);
    vec3d f1 = atDerU(u,v);
    vec3d f2 = atDerUU(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}

double surface::distToDerVV(double u, double v, vec3d A) {
    vec3d f = at(u,v);
    vec3d f1 = atDerV(u,v);
    vec3d f2 = atDerVV(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}

double surface::distToDerUV(double u, double v, vec3d A) {
    vec3d f = at(u, v);
    vec3d f1u = atDerU(u, v);
    vec3d f1v = atDerV(u, v);
    vec3d f2 = atDerUV(u, v);

    double dist = ((((f - A)*(f - A)).sum())*((f2*(f-A)).sum() + (f1u*f1v).sum()) - (f1u*(f-A)).sum()*(f1v*(f-A)).sum())/(pow(((f-A)*(f-A)).sum(),1.5));

    return dist;
}
