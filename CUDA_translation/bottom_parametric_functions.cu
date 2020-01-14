#include "../model.h"
#include "vec3d.h"
#include "gpumodel.h"


__device__
vec3d at(double u, double v) {
    vec3d point = (1 - v) * bendCurve.f(u) + v * leafCurve.f(u) +
                  (1 - u) * symCurve.f(1 - v) + u * sinCurve.f(1 - v) -
                  ((1 - u) * (1 - v) * Q_b_sym + u * v * Q_l_sin +
                   u * (1 - v) * Q_b_sin + (1 - u) * v * Q_l_sym);

    return point;
}

__device__
vec3d atDerU(double u, double v) {
    vec3d point = (1 - v) * bendCurve.df(u) + v * leafCurve.df(u)
                  - symCurve.f(1 - v) + sinCurve.f(1 - v) -
                  (- (1 - v) * Q_b_sym + v * Q_l_sin +
                   (1 - v) * Q_b_sin - v * Q_l_sym);

    return point;
}

__device__
vec3d atDerV(double u, double v) {
    vec3d point = - bendCurve.f(u) + leafCurve.f(u) -
                  (1 - u) * symCurve.df(1 - v) - u * sinCurve.df(1 - v) -
                  (- (1 - u) * Q_b_sym + u * Q_l_sin
                   - u * Q_b_sin + (1 - u) * Q_l_sym);

    return point;
}

__device__
vec3d atDerUU(double u, double v) {
    vec3d point = (1 - v) * bendCurve.ddf(u) + v * leafCurve.ddf(u);

    return point;
}

__device__
vec3d atDerVV(double u, double v) {
    vec3d point = (1 - u) * symCurve.ddf(1 - v) + u * sinCurve.ddf(1 - v);

    return point;
}

__device__
vec3d atDerUV(double u, double v) {
    vec3d point = - bendCurve.df(u) + leafCurve.df(u)
                  + symCurve.df(1 - v) - sinCurve.df(1 - v) -
                  (Q_b_sym + Q_l_sin - Q_b_sin - Q_l_sym);

    return point;
}
