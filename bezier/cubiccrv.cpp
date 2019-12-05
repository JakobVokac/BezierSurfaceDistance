//
// Created by s3179222 on 12/4/19.
//

#include "cubiccrv.h"

vec3d cubiccrv::f(double t) {
    return c0 * pow((1 - t), 3.0) + c1 * 3 * pow((1 - t), 2.0) * t +
           c2 * 3 * (1 - t) * pow(t, 2.0) + c3 * pow(t, 3.0);
};
vec3d cubiccrv::df(double t) {
    return (c1 - c0) * 3 * pow((1 - t), 2.0) + (c2 - c1) * 6 * t * (1 - t) + (c3 - c2) * 3 * pow(t, 2.0);
};
vec3d cubiccrv::ddf(double t){
    return (c2 - c1 * 2 + c0) * 6 * (1 - t) + (c3 - c2 * 2 + c1) * 6 * t;
};
