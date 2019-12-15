//
// Created by s3179222 on 12/4/19.
//

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
    void subdivideAt(double t, cubiccrv &crv1, cubiccrv &crv2);
};



#endif //HEARTVALVEMODEL_CUBICCRV_H
