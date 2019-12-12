//
// Created by s3179222 on 12/8/19.
//

#include "curve.h"

double curve::sqDistTo(double t, class vec3d A){
    vec3d P = f(t);

    double dist = P.sqdist(A);
    return dist;
}

double curve::sqDistToDer1(double t, class vec3d A){
    vec3d P = f(t);
    vec3d Pd = df(t);
    double dist = 2*((P-A)*Pd).sum();
    return dist;
}

double curve::sqDistToDer2(double t, class vec3d A){
    vec3d P = f(t);
    vec3d Pd = df(t);
    vec3d Pdd = ddf(t);

    double dist = 2*(Pd*Pd + (P-A)*Pdd).sum();

    return dist;
}

double curve::distTo(double t, class vec3d A){
    return sqrt(sqDistTo(t,A));
}

std::ostream &operator<<(std::ostream &out, curve &c) {
    return out << c.asString();
}

