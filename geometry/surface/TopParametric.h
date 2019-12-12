//
// Created by s3179222 on 12/8/19.
//

#ifndef HEARTVALVEMODEL_TOPPARAMETRIC_H
#define HEARTVALVEMODEL_TOPPARAMETRIC_H


#include "surface.h"

class TopParametric : public surface {
public:
    TopParametric() = default;
    TopParametric(
            cubiccrv leafCurve,
            cubiccrv bendCurve,
            cubiccrv symCurve,
            cubiccrv sinCurve,
            vec3d Q_l_sin,
            vec3d Q_l_sym,
            vec3d Q_b_sin,
            vec3d Q_b_sym
    );

    ~TopParametric() = default;

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

private:
    cubiccrv leafCurve, bendCurve, symCurve, sinCurve;

    vec3d Q_l_sin{}, Q_l_sym{}, Q_b_sin{}, Q_b_sym{};
};


#endif //HEARTVALVEMODEL_TOPPARAMETRIC_H
