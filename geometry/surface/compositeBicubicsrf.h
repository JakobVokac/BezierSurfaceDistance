//
// Created by s3179222 on 12/16/19.
//

#ifndef HEARTVALVEMODEL_COMPOSITEBICUBICSRF_H
#define HEARTVALVEMODEL_COMPOSITEBICUBICSRF_H


#include "BezierSurface.h"

class compositeBicubicsrf : public BezierSurface {
private:
    BezierSurface *topLeft, *topRight, *bottomLeft, *bottomRight;
    curve *v0{},*v1{},*u0{},*u1{};
public:
    compositeBicubicsrf() = default;
    explicit compositeBicubicsrf(BezierSurface *sur){
        sur->subdivide(topLeft,topRight,bottomLeft,bottomRight);
        *u0 = sur->edgeU0();
        *v0 = sur->edgeV0();
        *u1 = sur->edgeU1();
        *v1 = sur->edgeV1();
    }
    ~compositeBicubicsrf() override = default;

    vec3d at(double u, double v) override;
    vec3d atDerU(double u, double v) override;
    vec3d atDerV(double u, double v) override;
    vec3d atDerUU(double u, double v) override;
    vec3d atDerVV(double u, double v) override;
    vec3d atDerUV(double u, double v) override;

    curve & edgeU0() override;
    curve & edgeU1() override;
    curve & edgeV0() override;
    curve & edgeV1() override;


    bool hasValidControlNet() override;
    void controlNetSubdivide() override;
    bool closestPointInPatch(vec3d P) override;
    void subdivide(BezierSurface *tl,  BezierSurface *tr, BezierSurface *bl, BezierSurface *br) override;

    void recursiveControlNetSubdivide(BezierSurface *sur);
};

#endif //HEARTVALVEMODEL_COMPOSITEBICUBICSRF_H
