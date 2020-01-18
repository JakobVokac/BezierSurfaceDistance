//
// Created by s3179222 on 12/8/19.
//
/*
 * This is the implementation of the top part of each leaf in the model. The class contains all of the corners and edges
 * (see model.h), which include the bending curve, the leaf curve, the top half of the symmetry curve and the top half
 * of the sinusoidal curve.
 *
 * The surface itself is a bilinear interpolation of the 4 edges.
 */
#include "TopParametric.h"

TopParametric::TopParametric(cubiccrv leafCurve, cubiccrv bendCurve, cubiccrv symCurve, cubiccrv sinCurve,
                             vec3d Q_l_sin, vec3d Q_l_sym, vec3d Q_b_sin, vec3d Q_b_sym) {
    this->leafCurve = leafCurve;
    this->bendCurve = bendCurve;
    this->sinCurve = sinCurve;
    this->symCurve = symCurve;

    this->Q_l_sin = Q_l_sin;
    this->Q_l_sym = Q_l_sym;
    this->Q_b_sin = Q_b_sin;
    this->Q_b_sym = Q_b_sym;
}

vec3d TopParametric::at(double u, double v) {
    vec3d point = (1 - v) * bendCurve.f(u) + v * leafCurve.f(u) +
                  (1 - u) * symCurve.f(1 - v) + u * sinCurve.f(1 - v) -
                  ((1 - u) * (1 - v) * Q_b_sym + u * v * Q_l_sin +
                   u * (1 - v) * Q_b_sin + (1 - u) * v * Q_l_sym);

    return point;
}

vec3d TopParametric::atDerU(double u, double v) {
    vec3d point = (1 - v) * bendCurve.df(u) + v * leafCurve.df(u)
                  - symCurve.f(1 - v) + sinCurve.f(1 - v) -
                  (- (1 - v) * Q_b_sym + v * Q_l_sin +
                   (1 - v) * Q_b_sin - v * Q_l_sym);

    return point;
}

vec3d TopParametric::atDerV(double u, double v) {
    vec3d point = - bendCurve.f(u) + leafCurve.f(u) -
                  (1 - u) * symCurve.df(1 - v) - u * sinCurve.df(1 - v) -
                  (- (1 - u) * Q_b_sym + u * Q_l_sin
                   - u * Q_b_sin + (1 - u) * Q_l_sym);


    return point;
}

vec3d TopParametric::atDerUU(double u, double v) {
    vec3d point = (1 - v) * bendCurve.ddf(u) + v * leafCurve.ddf(u);

    return point;
}

vec3d TopParametric::atDerVV(double u, double v) {
    vec3d point = (1 - u) * symCurve.ddf(1 - v) + u * sinCurve.ddf(1 - v);

    return point;
}

vec3d TopParametric::atDerUV(double u, double v) {
    vec3d point = - bendCurve.df(u) + leafCurve.df(u)
                  + symCurve.df(1 - v) - sinCurve.df(1 - v) -
                  (Q_b_sym + Q_l_sin - Q_b_sin - Q_l_sym);

    return point;
}

curve & TopParametric::edgeU0() {
    return symCurve;
}

curve & TopParametric::edgeU1() {
    return sinCurve;
}

curve & TopParametric::edgeV0() {
    return bendCurve;
}

curve & TopParametric::edgeV1() {
    return leafCurve;
}


