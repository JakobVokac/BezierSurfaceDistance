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

double curve::distToDer1(double t, vec3d A){

    return sqDistToDer1(t,A) / distTo(t,A);
}

double curve::distToDer2(double t, vec3d A) {
    vec3d P = f(t);
    vec3d f1 = df(t);
    vec3d f2 = ddf(t);

    double dist = ((P - A)*f2 + f1*f1).sum()/(sqrt(((P-A)*(P-A)).sum()) - pow(((P-A)*f1).sum(),2.0)/(pow(((P-A)*(P-A)).sum(),1.5)));

    return dist;
}

std::ostream &operator<<(std::ostream &out, curve &c) {
    return out << c.asString();
}

