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

    vec3d l1 = c2 - c0, l2 = c3 - c1;

    double ll1 =l1.dot(l1), ll2 = l2.dot(l2);

    double t11 = l1.dot(c1 - c0) / ll1, t12 = l1.dot(c3 - c0) / ll1;
    double t21 = l2.dot(c2 - c1) / ll2, t22 = l2.dot(c0 - c1) / ll2;

    vec3d v11 = c0 + t11 * l1, v12 = c0 + t12 * l1;
    vec3d v21 = c1 + t21 * l2, v22 = c1 + t22 * l2;

    v11 = c1 - v11;
    v12 = c3 - v12;
    v21 = c2 - v21;
    v22 = c0 - v22;

    return !(v11.dot(v12) > 0 || v21.dot(v22) > 0);
}

bool cubiccrv::closestPointInCurve(vec3d P){

    double R1 = (P-c0).dot(c1-c0);
    double R2 = (c3-P).dot(c3-c2);
    double R3 = (c0-c3).dot(P-c3);
    double R4 = (c0-c3).dot(P-c0);

    return !((R1 < 0 || R2 < 0) && R3 * R4 > 0);
}

void cubiccrv::subdivide(double t, cubiccrv &crv1, cubiccrv &crv2) {
    vec3d c01{}, c12{}, c23{}, c02{}, c13{}, c03{};

    c01 = (1-t)*c0 + t*c1;
    c12 = (1-t)*c1 + t*c2;
    c23 = (1-t)*c2 + t*c3;
    c02 = (1-t)*c01 + t*c12;
    c13 = (1-t)*c12 + t*c23;
    c03 = (1-t)*c02 + t*c13;

    crv1 = {c0,c01,c02,c03};
    crv2 = {c03,c13,c23,c3};
}

