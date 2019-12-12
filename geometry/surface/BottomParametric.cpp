//
// Created by s3179222 on 12/8/19.
//

#include "BottomParametric.h"

BottomParametric::BottomParametric(cubiccrv sinCurve, cubiccrv symCurve, cubiccrv bendCurve, vec3d Q_r, vec3d Q_b_sin,
                                   vec3d Q_b_sym) {
    this->sinCurve = sinCurve;
    this->symCurve = symCurve;
    this->bendCurve = bendCurve;
    this->Q_r = Q_r;
    this->Q_b_sym = Q_b_sym;
    this->Q_b_sin = Q_b_sin;
    this->Q_corner = BottomCorner(Q_r);
}

vec3d BottomParametric::at(double u, double v) {
    vec3d res = v * Q_r + (1 - v) * bendCurve.f(u) + u * sinCurve.f(1 - v) + (1 - u) * symCurve.f(1 - v) -
                (u * v * Q_r + (1 - u) * (1 - v) * Q_b_sym + (1 - u) * v * Q_r + u * (1 - v) * Q_b_sin);
    return res;
}

vec3d BottomParametric::atDerU(double u, double v) {
    vec3d res = (1 - v) * bendCurve.df(u) + sinCurve.f(1 - v) - symCurve.f(1 - v) -
                (-(1 - v) * Q_b_sym + (1 - v) * Q_b_sin);
    return res;
}

vec3d BottomParametric::atDerV(double u, double v) {
    vec3d res = Q_r - bendCurve.f(u) - u * sinCurve.df(1 - v) - (1 - u) * symCurve.df(1 - v) -
                (u * Q_r - (1 - u) * Q_b_sym + (1 - u) * Q_r - u * Q_b_sin);
    return res;
}

vec3d BottomParametric::atDerUU(double u, double v) {
    vec3d res = (1 - v) * bendCurve.ddf(u);
    return res;
}

vec3d BottomParametric::atDerVV(double u, double v) {
    vec3d res = u * sinCurve.ddf(1 - v) + (1 - u) * symCurve.ddf(1 - v);
    return res;
}

vec3d BottomParametric::atDerUV(double u, double v) {
    vec3d res = - bendCurve.df(u) - sinCurve.f(1 - v) + symCurve.f(1 - v) -
                (Q_b_sym - Q_b_sin);
    return res;
}

curve & BottomParametric::edgeU0() {
    return symCurve;
}

curve & BottomParametric::edgeU1() {
    return sinCurve;
}

curve & BottomParametric::edgeV0() {
    return bendCurve;
}

curve & BottomParametric::edgeV1() {
    return Q_corner;
}

vec3d BottomParametric::BottomCorner::f(double t) {
    return Q_r;
}

vec3d BottomParametric::BottomCorner::df(double t) {
    return {1,1,1};
}

vec3d BottomParametric::BottomCorner::ddf(double t) {
    return {1,1,1};
}
