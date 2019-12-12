//
// Created by s3179222 on 10/12/19.
//

#ifndef HEARTVALVEMODEL_MODEL_H
#define HEARTVALVEMODEL_MODEL_H

#include "old_code/bezier.h"
#include <cmath>
#include "geometry/maths/vec3d.h"
#include "geometry/curve/cubiccrv.h"
#include "geometry/surface/surface.h"
#include "geometry/surface/TopParametric.h"
#include "geometry/surface/BottomParametric.h"
#include "geometry/surface/bicubicsrf.h"

class Model {
public:

    double root_radius;
    double leaflet_angle;
    double leaflet_height;
    double leaflet_radius_out;
    double leaflet_radius_in;
    double leaflet_power;

    double bending_angle;
    double bending_height;
    double bending_radius_out;
    double bending_radius_in;
    double bending_power;

    double s_sin;
    double s_sym;
    double xQ1_sym;
    double zQ1_sym;

    double leaflet_end_x;
    double leaflet_end_y;
    double bending_end_x;
    double bending_end_y;

    double eps;
    double t0;
    double t3;

    cubiccrv leafCurve, bendCurve, symCurveTop, symCurveBottom, sinCurveTop, sinCurveBottom;

    vec3d Q_l_sin{}, Q_l_sym{}, Q_b_sin{}, Q_b_sym{}, Q_r{};
    explicit Model(
            double root_radius = 12.0,
            double leaflet_angle = 59.5/180 * M_PI,
            double leaflet_height = 11.4,
            double leaflet_radius_out = 14.4,
            double leaflet_radius_in = 3.5,
            double leaflet_power = 1.35,
            double bending_angle = 50.0/180 * M_PI,
            double bending_height = 0.5,
            double bending_radius_out = 16.8,
            double bending_radius_in = 15.5,
            double bending_power = 1.35,
            double xQ1_sym = -0.2,
            double zQ1_sym = -0.2,
            double eps = 0.01,
            double t0 = 1.0,
            double t3 = 6.5
                    );

    double leafCurveFunc(double y);

    double bendCurveFunc(double y);

    cubiccrv leafCurveBezier();

    cubiccrv bendCurveBezier();

    cubiccrv sinCurveBottomBezier();

    cubiccrv sinCurveTopBezier();

    cubiccrv symCurveBottomBezier();

    cubiccrv symCurveTopBezier();

    vec3d fillTop(double u, double v);

    vec3d fillTopDerU(double u, double v);
    vec3d fillTopDerUU(double u, double v);
    vec3d fillTopDerUV(double u, double v);

    vec3d fillTopDerV(double u, double v);
    vec3d fillTopDerVV(double u, double v);

    vec3d fillBottom(double u, double v);

    double distanceToTopPoint(double u, double v, vec3d A);


    double completeDistanceTopDerU(double u, double v, vec3d A);
    double completeDistanceTopDerV(double u, double v, vec3d A);


    //old_code functions
    double completeDistanceTopDerUU(double u, double v, vec3d A);
    double completeDistanceTopDerVV(double u, double v, vec3d A);
    double completeDistanceTopDerUV(double u, double v, vec3d A);


    double distanceToBottomPoint(double u, double v, vec3d A);

    static Model getPart(Model model, int idx);

    double squaredTopDist(double u, double v, vec3d A);
    double squaredTopDistDerU(double u, double v, vec3d A);
    double squaredTopDistDerV(double u, double v, vec3d A);
    double squaredTopDistDerUU(double u, double v, vec3d A);
    double squaredTopDistDerVV(double u, double v, vec3d A);
    double squaredTopDistDerUV(double u, double v, vec3d A);

    TopParametric getTopParametric();
    BottomParametric getBottomParametric();
    bicubicsrf getTopBezier();
};


#endif
