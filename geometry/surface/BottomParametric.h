//
// Created by s3179222 on 12/8/19.
//
/*
 * This is the implementation of the top part of each leaf in the model. The class contains all of the corners and edges
 * (see model.h), which include the bending curve, the bottom half of the symmetry curve and the bottom half
 * of the sinusoidal curve.
 *
 * The surface is an interpolation between the 3 edges. It also has a subclass such that the bottom corner can be passed
 * as an edge, which is required by the surface interface. The bottom corner then outputs the same value every time.
 */

#ifndef HEARTVALVEMODEL_BOTTOMPARAMETRIC_H
#define HEARTVALVEMODEL_BOTTOMPARAMETRIC_H


#include "surface.h"

class BottomParametric : public surface {
public:
    BottomParametric() = default;
    BottomParametric(
            cubiccrv sinCurve,
            cubiccrv symCurve,
            cubiccrv bendCurve,
            vec3d Q_r,
            vec3d Q_b_sin,
            vec3d Q_b_sym
            );
    ~BottomParametric() = default;

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
    cubiccrv sinCurve, symCurve, bendCurve;
    vec3d Q_r{}, Q_b_sin{}, Q_b_sym{};

    class BottomCorner : public curve{
    public:
        explicit BottomCorner(vec3d Q_r){
            this->Q_r = Q_r;
        }
        ~BottomCorner() = default;
        vec3d f(double t) override;
        vec3d df(double t) override;
        vec3d ddf(double t) override;
        std::string asString() override{
            return "\nCorner "
            + std::to_string(Q_r.getx()) + " "
            + std::to_string(Q_r.gety()) + " "
            + std::to_string(Q_r.getz()) + "\n";
        };
        vec3d curvePlaneNormal() override{
            return {0,0,0};
        };

    private:
        vec3d Q_r{};
    };

    BottomCorner Q_corner = BottomCorner(vec3d());
};


#endif //HEARTVALVEMODEL_BOTTOMPARAMETRIC_H
