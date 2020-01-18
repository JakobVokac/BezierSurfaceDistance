//
// Created by s3179222 on 12/4/19.
//
/*
 * This is an implementatin of a 3rd degree Bezier curve. The implementation gives two constructors, one with 4 control
 * points (the usual way to instantiate it) and one with 3 control points (if you want to implement a 2nd degree Bezier
 * curve as a 3rd degree for a surface or optimization reason).
 *
 * Aside from the basic Bezier functions for curves, the class also contains functions for checking whether the control polygon
 * of the curve is simple and convex (when projected onto a plane containing the first, last and a middle control point).
 *
 * It also implements a function for checking whether the closest point on the curve to a point P sits on the inside or
 * on a corner.
 *
 *  The polygon and corner detection algorithms were taken and implemented from:
 *      Yingliang Ma & W. Terry Hewitt,
 *      "Point inversion and projection for NURBS curve: Control polygon approach",
 *      https://www.researchgate.net/publication/232638188_Point_inversion_and_projection_for_NURBS_curve_Control_polygon_approach
 */
#ifndef HEARTVALVEMODEL_CUBICCRV_H
#define HEARTVALVEMODEL_CUBICCRV_H

#include "curve.h"

class cubiccrv : public curve {
private:
    vec3d c0{},c1{},c2{},c3{};

public:
    cubiccrv() = default;

    cubiccrv(vec3d c0, vec3d c1, vec3d c2, vec3d c3){
        this->c0 = c0;
        this->c1 = c1;
        this->c2 = c2;
        this->c3 = c3;
    }

    cubiccrv(vec3d c0, vec3d c1, vec3d c2){
        this->c0 = c0;
        this->c1 = c0 * 1.0/3.0 + c1 * 2.0/3.0;
        this->c2 = c1 * 2.0/3.0 + c2 * 1.0/3.0;
        this->c3 = c2;
    }
    ~cubiccrv() = default;

    vec3d f(double t) override;
    vec3d df(double t) override;
    vec3d ddf(double t) override;

    vec3d curvePlaneNormal() override;
    vec3d getCtrlP(int i);

    cubiccrv RotAxisZ(double rad);
    cubiccrv MirrorY();

    std::string asString() override;

    bool hasConvexPolygon();
    bool closestPointInCurve(vec3d P);
    void subdivide(double t, cubiccrv &crv1, cubiccrv &crv2);
};



#endif //HEARTVALVEMODEL_CUBICCRV_H
