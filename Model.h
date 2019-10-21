//
// Created by s3179222 on 10/12/19.
//

#ifndef HEARTVALVEMODEL_MODEL_H
#define HEARTVALVEMODEL_MODEL_H

#include "bezier.h"
#include <cmath>

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

    CBezier leafCurve, bendCurve, symCurveTop, symCurveBottom;
    QBezier sinCurveTop, sinCurveBottom;

    vector<double> Q_l_sin, Q_l_sym, Q_b_sin, Q_b_sym, Q_r;
    Model(
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

    CBezier leafCurveBezier(double eps, double t0, double t3);

    CBezier bendCurveBezier(double eps, double t0, double t3);

    QBezier sinCurveBottomBezier();

    QBezier sinCurveTopBezier();

    CBezier symCurveBottomBezier();

    CBezier symCurveTopBezier();

    vector<double> fillTop(double u, double v);

    vector<double> fillTopDerU(double u, double v);
    vector<double> fillTopDerUU(double u, double v);
    vector<double> fillTopDerUV(double u, double v);

    vector<double> fillTopDerV(double u, double v);
    vector<double> fillTopDerVV(double u, double v);
    vector<double> fillTopDerVU(double u, double v);

    vector<double> fillBottom(double u, double v);

    double distanceToTopPoint(double u, double v, vector<double> A);

    double distanceToTopPointDerU(double u, double v, vector<double> A);
    double distanceToTopPointDerUU(double u, double v, vector<double> A);
    double distanceToTopPointDerUV(double u, double v, vector<double> A);

    double distanceToTopPointDerV(double u, double v, vector<double> A);
    double distanceToTopPointDerVV(double u, double v, vector<double> A);
    double distanceToTopPointDerVU(double u, double v, vector<double> A);

    double completeDistanceTopDerU(double u, double v, vector<double> A);
    double completeDistanceTopDerUU(double u, double v, vector<double> A);
    double completeDistanceTopDerV(double u, double v, vector<double> A);
    double completeDistanceTopDerVV(double u, double v, vector<double> A);
    double completeDistanceTopDerUV(double u, double v, vector<double> A);



    //test functions
    double completeDistanceTopDerUU2(double u, double v, vector<double> A);
    double completeDistanceTopDerVV2(double u, double v, vector<double> A);
    double completeDistanceTopDerUV2(double u, double v, vector<double> A);





    double distranceToBottomPoint(double u, double v, vector<double> A);

    static Model getPart(Model model, int idx);
};


#endif
