//
// Created by s3179222 on 12/4/19.
//

#include <limits>
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
}

cubiccrv cubiccrv::RotAxisZ(double rad) {
    return cubiccrv( c0.RotAxisZ(rad), c1.RotAxisZ(rad), c2.RotAxisZ(rad), c3.RotAxisZ(rad));
}

cubiccrv cubiccrv::MirrorY() {
    return cubiccrv(c0.MirrorY(), c1.MirrorY(), c2.MirrorY(), c3.MirrorY());
}

vec3d cubiccrv::getCtrlP(int i) {
    switch (i) {
        case 0:
            return c0;
        case 1:
            return c1;
        case 2:
            return c2;
        case 3:
            return c3;
        default:
            return vec3d();
    }
}

std::string cubiccrv::asString() {
    std::string s{};
    s += "cubiccrv; x y z\n";
    s += "C0: " + std::to_string(c0.getx()) + " " + std::to_string(c0.gety()) + " " + std::to_string(c0.getz()) + "\n";
    s += "C1: " + std::to_string(c1.getx()) + " " + std::to_string(c1.gety()) + " " + std::to_string(c1.getz()) + "\n";
    s += "C2: " + std::to_string(c2.getx()) + " " + std::to_string(c2.gety()) + " " + std::to_string(c2.getz()) + "\n";
    s += "C3: " + std::to_string(c3.getx()) + " " + std::to_string(c3.gety()) + " " + std::to_string(c3.getz()) + "\n";

    return s;
}

vec3d cubiccrv::curvePlaneNormal() {
    vec3d c0c1 = c1-c0;
    vec3d c0c2 = c2-c0;
    vec3d c0c3 = c3-c0;

    std::cout << " c0c1 " << c0c1 << std::endl;
    std::cout << " c0c2 " << c0c2 << std::endl;
    std::cout << " c0c3 " << c0c3 << std::endl;

    vec3d N = c0c1.cross(c0c2);

    std::cout << " N " << N << std::endl;

    N /= N.mag();

    double dot = N.dot(c0c3);

    std::cout << " N . c0c3 " << dot << std::endl;

    if(dot > 10*std::numeric_limits<double>::epsilon())
        return {0,0,0};

    return N;
}

bool cubiccrv::hasConvexPolygon() {
    vec3d n1 = c2-c0, n2 = c3-c1;
    double nn1 =n1.dot(n1), nn2 = n2.dot(n2);
    double t11 = n1.dot(c1-c0)/nn1, t12 = n1.dot(c3-c0)/nn1;
    double t21 = n2.dot(c2-c1)/nn2, t22 = n2.dot(c0-c1)/nn2;

    vec3d v11 = c0 + t11*n1, v12 = c0 + t12*n1;
    vec3d v21 = c1 + t21*n2, v22 = c1 + t22*n2;

    return !(v11.dot(v12) > 0 || v21.dot(v22) > 0);
}


