//
// Created by s3179222 on 10/12/19.
//

#include "model.h"

Model::Model(
        double root_radius,
        double leaflet_angle,
        double leaflet_height,
        double leaflet_radius_out,
        double leaflet_radius_in,
        double leaflet_power,
        double bending_angle,
        double bending_height,
        double bending_radius_out,
        double bending_radius_in,
        double bending_power,
        double xQ1_sym_var,
        double zQ1_sym_var,
        double eps,
        double t0,
        double t3
){

    this->root_radius = root_radius;
    this->leaflet_angle = leaflet_angle;
    this->leaflet_height = leaflet_height;
    this->leaflet_radius_out = leaflet_radius_out;
    this->leaflet_radius_in = leaflet_radius_in;
    this->leaflet_power = leaflet_power;
    this->bending_angle = bending_angle;
    this->bending_height = bending_height;
    this->bending_radius_out = bending_radius_out;
    this->bending_radius_in = bending_radius_in;
    this->bending_power = bending_power;
    this->xQ1_sym = xQ1_sym_var;
    this->zQ1_sym = zQ1_sym_var;
    this->s_sin = bending_height / leaflet_height;
    this->s_sym = bending_height / leaflet_height;

    this->leaflet_end_x = leaflet_radius_out*cos(leaflet_angle);
    this->leaflet_end_y = leaflet_radius_out*sin(leaflet_angle);
    this->bending_end_x = bending_radius_out*cos(bending_angle);
    this->bending_end_y = bending_radius_out*sin(bending_angle);

    this->Q_l_sin = {leaflet_end_x, leaflet_end_y, leaflet_height};
    this->Q_l_sym = {leaflet_radius_in, 0.0, leaflet_height};
    this->Q_b_sin = {bending_end_x, bending_end_y, bending_height};
    this->Q_b_sym = {bending_radius_in, 0.0, bending_height};
    this->Q_r = {root_radius,0.0,0.0};

    this->eps = eps;
    this->t0 = t0;
    this->t3 = t3;

    this->leafCurve = leafCurveBezier(eps, t0, t3);
    this->bendCurve = bendCurveBezier(eps, t0, t3);
    this->symCurveBottom = symCurveBottomBezier();
    this->symCurveTop = symCurveTopBezier();
    this->sinCurveBottom = sinCurveBottomBezier();
    this->sinCurveTop = sinCurveTopBezier();


}
double Model::leafCurveFunc(double y){
    return leaflet_radius_in + (leaflet_end_x - leaflet_radius_in) * pow(y/leaflet_end_y, leaflet_power);
}

double Model::bendCurveFunc(double y){
    return bending_radius_in + (bending_end_x - bending_radius_in) * pow(y/bending_end_y, bending_power);
}

CBezier Model::leafCurveBezier(double eps, double t0, double t3){
    CBezier curve = CBezier();
    curve.C0 = {leaflet_radius_in,0,leaflet_height};
    curve.C3 = {leaflet_end_x, leaflet_end_y, leaflet_height};
    curve.C1 = {leaflet_radius_in + t0 * (leafCurveFunc(eps) - leaflet_radius_in)/eps, 0 + t0, leaflet_height};
    curve.C2 = {leaflet_end_x + t3 * (leafCurveFunc(leaflet_end_y - eps) - leaflet_end_x)/eps, leaflet_end_y - t3, leaflet_height};

    return curve;
}

CBezier Model::bendCurveBezier(double eps, double t0, double t3){
    CBezier curve = CBezier();
    curve.C0 = {bending_radius_in,0,bending_height};
    curve.C1 = {bending_radius_in + t0 * (bendCurveFunc(eps) - bending_radius_in)/eps, 0 + t0, bending_height};
    curve.C2 = {bending_end_x + t3 * (bendCurveFunc(bending_end_y - eps) - bending_end_x)/eps, bending_end_y - t3, bending_height};
    curve.C3 = {bending_end_x, bending_end_y, bending_height};

    return curve;
}

QBezier Model::sinCurveBottomBezier() {
    CBezier temp = CBezier();
    QBezier curve = QBezier();

    temp.C0 = {root_radius, 0.0, 0.0};
    temp.C2 = {leaflet_end_x, leaflet_end_y, leaflet_height};
    temp.C3 = {bending_end_x, bending_end_y, bending_height};
    for (int i = 0; i < 3; i++) {
        temp.C1[i] = (temp.C3[i] - pow((1 - s_sin), 2.0) * temp.C0[i] - pow(s_sin, 2.0) * temp.C2[i]) /
                     (2 * (1 - s_sin) * s_sin);
    }

    curve.C0 = temp.C0;
    curve.C2 = temp.C3;
    for (int i = 0; i < 3; i++) {
        curve.C1[i] = (1-s_sin)*temp.C0[i] + s_sin*temp.C1[i];
    }
    return curve;
}

QBezier Model::sinCurveTopBezier() {
    CBezier temp = CBezier();
    QBezier curve = QBezier();

    temp.C0 = {root_radius, 0.0, 0.0};
    temp.C2 = {leaflet_end_x, leaflet_end_y, leaflet_height};
    temp.C3 = {bending_end_x, bending_end_y, bending_height};
    for (int i = 0; i < 3; i++) {
        temp.C1[i] = (temp.C3[i] - pow((1 - s_sin), 2.0) * temp.C0[i] - pow(s_sin, 2.0) * temp.C2[i]) /
                     (2 * (1 - s_sin) * s_sin);
    }

    curve.C0 = temp.C2;
    curve.C2 = temp.C3;
    for (int i = 0; i < 3; i++) {
        curve.C1[i] = (1-s_sin)*temp.C1[i] + s_sin*temp.C2[i];
    }
    return curve;
}

CBezier Model::symCurveBottomBezier(){
    CBezier temp = CBezier();
    CBezier curve = CBezier();
    vector<double> B = vector<double>(3);
    vector<double> R12 = vector<double>(3);

    temp.C0 = {root_radius,0.0,0.0};
    temp.C3 = {leaflet_radius_in,0.0,leaflet_height};
    temp.C1 = {(root_radius + bending_radius_in)/2 + xQ1_sym*(root_radius - bending_radius_in), 0.0, bending_height/2 + zQ1_sym*bending_height};
    B = {bending_radius_in,0,bending_height};
    for (int i = 0; i < 3; ++i) {
        temp.C2[i] = (B[i] - pow((1 - s_sym), 3.0) * temp.C0[i] - 3 * s_sym * pow((1 - s_sym), 2.0) * temp.C1[i] -
                      pow(s_sym, 3.0) * temp.C3[i]) / (3 * (1 - s_sym) * pow(s_sym, 2.0));
    }

    curve.C0 = temp.C0;
    curve.C3 = B;
    for (int i = 0; i < 3; ++i) {
        curve.C1[i] = (1-s_sym)*temp.C0[i] + s_sym*temp.C1[i];
        R12[i] = (1-s_sym)*temp.C1[i] + s_sym*temp.C2[i];
        curve.C2[i] = (1-s_sym)*curve.C1[i] + s_sym*R12[i];
    }

    return curve;
}

CBezier Model::symCurveTopBezier(){
    CBezier temp = CBezier();
    CBezier curve = CBezier();
    vector<double> B = vector<double>(3);
    vector<double> R12 = vector<double>(3);

    temp.C0 = {root_radius,0,0};
    temp.C3 = {leaflet_radius_in,0,leaflet_height};
    temp.C1 = {(root_radius + bending_radius_in)/2 + xQ1_sym*(root_radius - bending_radius_in), 0, bending_height/2 + zQ1_sym*bending_height};
    B = {bending_radius_in,0,bending_height};
    for (int i = 0; i < 3; ++i) {
        temp.C2[i] = (B[i] - pow((1 - s_sym), 3) * temp.C0[i] - 3 * s_sym * pow((1 - s_sym), 2.0) * temp.C1[i] -
                      pow(s_sym, 3.0) * temp.C3[i]) / (3 * (1 - s_sym) * pow(s_sym, 2.0));
    }

    curve.C0 = temp.C3;
    curve.C3 = B;
    for (int i = 0; i < 3; ++i) {
        curve.C1[i] = (1-s_sym)*temp.C2[i] + s_sym*temp.C3[i];
        R12[i] = (1-s_sym)*temp.C1[i] + s_sym*temp.C2[i];
        curve.C2[i] = (1-s_sym)*R12[i] + s_sym*curve.C1[i];
    }

    return curve;
}

vector<double> Model::fillTop(double u, double v){
    vector<double> point = vector<double>(3);

    for (int i = 0; i < 3; ++i) {
        point[i] = (1 - v) * cubicBezier(u, i, bendCurve) + v * cubicBezier(u, i, leafCurve) +
                   (1 - u) * cubicBezier(1 - v, i, symCurveTop) + u * quadBezier(1 - v, i, sinCurveTop) -
                   ((1 - u) * (1 - v) * Q_b_sym[i] + u * v * Q_l_sin[i] +
                    u * (1 - v) * Q_b_sin[i] + (1 - u) * v * Q_l_sym[i]);
    }

    return point;
}

vector<double> Model::fillTopDerU(double u, double v) {
    vector<double> point = vector<double>(3);

    for (int i = 0; i < 3; ++i) {
        point[i] = (1 - v) * cubicBezierDer1(u, i, bendCurve) + v * cubicBezierDer1(u, i, leafCurve)
                 - cubicBezier(1 - v, i, symCurveTop) + quadBezier(1 - v, i, sinCurveTop) -
                   (- (1 - v) * Q_b_sym[i] + v * Q_l_sin[i] +
                   (1 - v) * Q_b_sin[i] - v * Q_l_sym[i]);
    }

    return point;
}

vector<double> Model::fillTopDerUU(double u, double v){
    vector<double> point = vector<double>(3);

    for (int i = 0; i < 3; ++i) {
        point[i] = (1 - v) * cubicBezierDer2(u, i, bendCurve) + v * cubicBezierDer2(u, i, leafCurve);
    }

    return point;
}

vector<double> Model::fillTopDerUV(double u, double v) {
    vector<double> point = vector<double>(3);

    for (int i = 0; i < 3; ++i) {
        point[i] = - cubicBezierDer1(u, i, bendCurve) + cubicBezierDer1(u, i, leafCurve)
                   + cubicBezierDer1(1 - v, i, symCurveTop) - quadBezierDer1(1 - v, i, sinCurveTop) -
                   ( + Q_b_sym[i] + Q_l_sin[i] - Q_b_sin[i] - Q_l_sym[i]);
    }

    return point;
}


vector<double> Model::fillTopDerV(double u, double v){
    vector<double> point = vector<double>(3);

    for (int i = 0; i < 3; ++i) {
        point[i] = - cubicBezier(u, i, bendCurve) + cubicBezier(u, i, leafCurve) -
                   (1 - u) * cubicBezierDer1(1 - v, i, symCurveTop) - u * quadBezierDer1(1 - v, i, sinCurveTop) -
                   (- (1 - u) * Q_b_sym[i] + u * Q_l_sin[i]
                    - u * Q_b_sin[i] + (1 - u) * Q_l_sym[i]);
    }

    return point;
}

vector<double> Model::fillTopDerVV(double u, double v){
    vector<double> point = vector<double>(3);

    for (int i = 0; i < 3; ++i) {
        point[i] = (1 - u) * cubicBezierDer2(1 - v, i, symCurveTop) + u * quadBezierDer2(1 - v, i, sinCurveTop);
    }

    return point;
}

vector<double> Model::fillTopDerVU(double u, double v){
    vector<double> point = vector<double>(3);

    for (int i = 0; i < 3; ++i) {
        point[i] = - cubicBezierDer1(u, i, bendCurve) + cubicBezierDer1(u, i, leafCurve)
                   + cubicBezierDer1(1 - v, i, symCurveTop) - quadBezierDer1(1 - v, i, sinCurveTop) -
                   ( + Q_b_sym[i] + Q_l_sin[i] - Q_b_sin[i] - Q_l_sym[i]);
    }

    return point;
}



vector<double> Model::fillBottom(double u, double v){
    vector<double> point = vector<double>(3);

    for (int i = 0; i < 3; ++i) {
        point[i] = v * Q_r[i] + (1 - v) * cubicBezier(u, i, bendCurve) +
                   u * quadBezier(1 - v, i, sinCurveBottom) + (1 - u) * cubicBezier(1 - v, i, symCurveBottom) -
                   (u * v * Q_r[i] + (1 - u) * (1 - v) * Q_b_sym[i] +
                    (1 - u) * v * Q_r[i] + u * (1 - v) * Q_b_sin[i]);
    }

    return point;
}

double Model::squaredTopDist(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);

    double distance = pow(top[0] - A[0],2.0) + pow(top[1] - A[1],2.0) + pow(top[2] - A[2],2.0);
    return distance;
}

double Model::squaredTopDistDerU(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDer = fillTopDerU(u,v);
    double distance = 2*((top[0] - A[0])*topDer[0] + (top[1] - A[1])*topDer[1] + (top[2] - A[2])*topDer[2]);
    return distance;
}
double Model::squaredTopDistDerV(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDer = fillTopDerV(u,v);
    double distance = 2*((top[0] - A[0])*topDer[0] + (top[1] - A[1])*topDer[1] + (top[2] - A[2])*topDer[2]);
    return distance;
}
double Model::squaredTopDistDerUU(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDer1 = fillTopDerU(u,v);
    vector<double> topDer2 = fillTopDerUU(u,v);

    double distance = 2*(pow((topDer1[0]),2.0) + (top[0] + A[0])*topDer2[0]
                      + pow((topDer1[1]),2.0) + (top[1] + A[1])*topDer2[1]
                      + pow((topDer1[2]),2.0) + (top[2] + A[2])*topDer2[2]);
    return distance;
}
double Model::squaredTopDistDerVV(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDer1 = fillTopDerV(u,v);
    vector<double> topDer2 = fillTopDerVV(u,v);

    double distance = 2*(pow((topDer1[0]),2.0) + (top[0] + A[0])*topDer2[0]
                         + pow((topDer1[1]),2.0) + (top[1] + A[1])*topDer2[1]
                         + pow((topDer1[2]),2.0) + (top[2] + A[2])*topDer2[2]);
    return distance;
}
double Model::squaredTopDistDerUV(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDerU = fillTopDerU(u,v);
    vector<double> topDerV = fillTopDerV(u,v);
    vector<double> topDerUV = fillTopDerUV(u,v);

    double distance = 2*((top[0] - A[0])*topDerUV[0] + topDerV[0]*topDerU[0]
                      + (top[1] - A[1])*topDerUV[1] + topDerV[1]*topDerU[1]
                      + (top[2] - A[2])*topDerUV[2] + topDerV[2]*topDerU[2]);
    return distance;
}

double Model::distanceToTopPoint(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);

    double distance = sqrt(pow(top[0] - A[0],2.0) + pow(top[1] - A[1],2.0) + pow(top[2] - A[2],2.0));
    return distance;
}


double Model::distanceToTopPointDerU(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDer = fillTopDerU(u,v);
    double distance = (top[0] - A[0])*topDer[0] + (top[1] - A[1])*topDer[1] + (top[2] - A[2])*topDer[2];
    return distance;
}

double Model::distanceToTopPointDerUU(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDer1 = fillTopDerU(u,v);
    vector<double> topDer2 = fillTopDerUU(u,v);

    double distance = pow((topDer1[0]),2.0) + (top[0] + A[0])*topDer2[0]
                    + pow((topDer1[1]),2.0) + (top[1] + A[1])*topDer2[1]
                    + pow((topDer1[2]),2.0) + (top[2] + A[2])*topDer2[2];
    return distance;
}


double Model::distanceToTopPointDerUV(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDerU = fillTopDerU(u,v);
    vector<double> topDerV = fillTopDerV(u,v);
    vector<double> topDerUV = fillTopDerUV(u,v);

    double distance = (top[0] - A[0])*topDerUV[0] + topDerV[0]*topDerU[0]
                      + (top[1] - A[1])*topDerUV[1] + topDerV[1]*topDerU[1]
                      + (top[2] - A[2])*topDerUV[2] + topDerV[2]*topDerU[2];
    return distance;
}


double Model::distanceToTopPointDerV(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDer = fillTopDerV(u,v);
    double distance = (top[0] - A[0])*topDer[0] + (top[1] - A[1])*topDer[1] + (top[2] - A[2])*topDer[2];
    return distance;
}

double Model::distanceToTopPointDerVV(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDer1 = fillTopDerV(u,v);
    vector<double> topDer2 = fillTopDerVV(u,v);

    double distance = pow((topDer1[0]),2.0) + (top[0] + A[0])*topDer2[0]
                    + pow((topDer1[1]),2.0) + (top[1] + A[1])*topDer2[1]
                    + pow((topDer1[2]),2.0) + (top[2] + A[2])*topDer2[2];
    return distance;
}

double Model::distanceToTopPointDerVU(double u, double v, vector<double> A){
    vector<double> top = fillTop(u,v);
    vector<double> topDerV = fillTopDerV(u,v);
    vector<double> topDerU = fillTopDerU(u,v);
    vector<double> topDerVU = fillTopDerVU(u,v);

    double distance = (top[0] - A[0])*topDerVU[0] + topDerV[0]*topDerU[0]
                    + (top[1] - A[1])*topDerVU[1] + topDerV[1]*topDerU[1]
                    + (top[2] - A[2])*topDerVU[2] + topDerV[2]*topDerU[2];
    return distance;

    // d(u,v)/dv = (f(u,v) - A)*f(u,v)/dv
    // d(u,v)/dvdu = f(u,v)/du * f(u,v)/dv + ((f(u,v) - A)*f(u,v)/dvdu
}

double Model::completeDistanceTopDerU(double u, double v, vector<double> A){

    return distanceToTopPointDerU(u,v,A)/distanceToTopPoint(u,v,A);
}
double Model::completeDistanceTopDerV(double u, double v, vector<double> A){

    return distanceToTopPointDerV(u,v,A)/distanceToTopPoint(u,v,A);
}
double Model::completeDistanceTopDerUU(double u, double v, vector<double> A){
    double f = distanceToTopPoint(u,v,A),
           fodu = distanceToTopPointDerU(u,v,A),
           fodudu = distanceToTopPointDerUU(u,v,A);

    return (f*fodudu - fodu/f*fodu)/(f*f);
}

double Model::completeDistanceTopDerUU2(double u, double v, vector<double> A) {
    vector<double> f = fillTop(u,v);
    vector<double> f1 = fillTopDerU(u,v);
    vector<double> f2 = fillTopDerUU(u,v);

    return (2*(f[0] - A[0])*f2[0] + 2*f1[0]*f1[0] + 2*(f[1] - A[1])*f2[1] + 2*f1[1]*f1[1] + 2*(f[2] - A[2])*f2[2] + 2*f1[2]*f1[2]) /
            (2*sqrt((f[0] - A[0])*(f[0] - A[0]) + (f[1] -  A[1])*(f[1] -  A[1]) + (f[2] - A[2])*(f[2] - A[2]))) -
            pow((2*(f[0] - A[0])*f1[0] + 2*(f[1] -  A[1])*f1[1] + 2*(f[2] - A[2])*f1[2]), 2.0)/
                    (4 * pow((f[0] - A[0])*(f[0] - A[0]) + (f[1] -  A[1])*(f[1] -  A[1]) + (f[2] - A[2])*(f[2] - A[2]),1.5));
}

double Model::completeDistanceTopDerVV2(double u, double v, vector<double> A) {
    vector<double> f = fillTop(u,v);
    vector<double> f1 = fillTopDerV(u,v);
    vector<double> f2 = fillTopDerVV(u,v);

    return (2*(f[0] - A[0])*f2[0] + 2*f1[0]*f1[0] + 2*(f[1] - A[1])*f2[1] + 2*f1[1]*f1[1] + 2*(f[2] - A[2])*f2[2] + 2*f1[2]*f1[2]) /
           (2*sqrt((f[0] - A[0])*(f[0] - A[0]) + (f[1] -  A[1])*(f[1] -  A[1]) + (f[2] - A[2])*(f[2] - A[2]))) -
           pow((2*(f[0] - A[0])*f1[0] + 2*(f[1] -  A[1])*f1[1] + 2*(f[2] - A[2])*f1[2]), 2.0)/
           (4 * pow((f[0] - A[0])*(f[0] - A[0]) + (f[1] -  A[1])*(f[1] -  A[1]) + (f[2] - A[2])*(f[2] - A[2]),1.5));
}

double Model::completeDistanceTopDerUV2(double u, double v, vector<double> A) {
    vector<double> f = fillTop(u, v);
    vector<double> f1u = fillTopDerU(u, v);
    vector<double> f1v = fillTopDerV(u, v);
    vector<double> f2 = fillTopDerUV(u, v);

    return (2 * (pow((f[0] - A[0]), 2.0) + pow((f[1] - A[1]), 2.0) + pow((f[2] - A[2]), 2.0)) *
             (f1v[0] * f1u[0] + (f[0] - A[0]) * f2[0] + f1v[1] * f1u[1] + (f[1] -A[1]) * f2[1] + f1v[2] * f1u[2] +
              (f[2] - A[2]) * f2[2])
             - 2 * ((f[0] - A[0]) * f1v[0] + (f[1] - A[1]) * f1v[1] + (f[2] - A[2]) * f1v[2]) *
               ((f[0] - A[0]) * f1u[0] + (f[1] - A[1]) * f1u[1] + (f[2] - A[2]) * f1u[2])) /
            (2 * pow((pow((f[0] - A[0]), 2.0) + pow((f[1] - A[1]), 2.0) + pow((f[2] - A[2]), 2.0)), (3.0 / 2.0)));

}
double Model::completeDistanceTopDerVV(double u, double v, vector<double> A){
    double f = distanceToTopPoint(u,v,A),
           fodu = distanceToTopPointDerV(u,v,A),
           fodudu = distanceToTopPointDerVV(u,v,A);

    return (f*fodudu - (fodu/f)*fodu)/(f*f);
}
double Model::completeDistanceTopDerUV(double u, double v, vector<double> A){
    double f = distanceToTopPoint(u,v,A),
            fodu = distanceToTopPointDerU(u,v,A),
            fodv = distanceToTopPointDerV(u,v,A),
            fodudv = distanceToTopPointDerUV(u,v,A);

    return (f*fodudv - fodv/f*fodu)/(f*f);
}

double Model::distranceToBottomPoint(double u, double v, vector<double> A){
    vector<double> bottom = fillBottom(u, v);

    double distance = sqrt(pow(bottom[0] - A[0], 2.0) + pow(bottom[1] - A[1], 2.0) + pow(bottom[2] - A[2], 2.0));
    return distance;
}

Model Model::getPart(Model model, int idx){
    if(idx % 2 != 0){
        model.leafCurve = MirrorY(model.leafCurve);
        model.bendCurve = MirrorY(model.bendCurve);
        model.sinCurveBottom = MirrorY(model.sinCurveBottom);
        model.sinCurveTop = MirrorY(model.sinCurveTop);
        model.Q_b_sym[1] *= -1;
        model.Q_l_sym[1] *= -1;
        model.Q_b_sin[1] *= -1;
        model.Q_l_sin[1] *= -1;
        model.Q_r[1] *= -1;
        idx--;
    }
    idx /= 2;
    for(int i = 0; i < idx; i++){
        model.leafCurve = RotAxisZ(120.0/180*M_PI,model.leafCurve);
        model.bendCurve = RotAxisZ(120.0/180*M_PI,model.bendCurve);
        model.symCurveBottom = RotAxisZ(120.0/180*M_PI,model.symCurveBottom);
        model.symCurveTop = RotAxisZ(120.0/180*M_PI,model.symCurveTop);
        model.sinCurveBottom = RotAxisZ(120.0/180*M_PI,model.sinCurveBottom);
        model.sinCurveTop = RotAxisZ(120.0/180*M_PI,model.sinCurveTop);
        model.Q_b_sym = RotAxisZ(120.0/180*M_PI,model.Q_b_sym);
        model.Q_l_sym = RotAxisZ(120.0/180*M_PI,model.Q_l_sym);
        model.Q_b_sin = RotAxisZ(120.0/180*M_PI,model.Q_b_sin);
        model.Q_l_sin = RotAxisZ(120.0/180*M_PI,model.Q_l_sin);
        model.Q_r = RotAxisZ(120.0/180*M_PI,model.Q_r);
    }
    return model;
}