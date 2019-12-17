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
    this->leafCurve = leafCurveBezier();
    this->bendCurve = bendCurveBezier();
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

cubiccrv Model::leafCurveBezier(){
    vec3d C0 = {leaflet_radius_in,0,leaflet_height};
    vec3d C3 = {leaflet_end_x, leaflet_end_y, leaflet_height};
    vec3d C1 = {leaflet_radius_in + t0 * (leafCurveFunc(eps) - leaflet_radius_in)/eps, 0 + t0, leaflet_height};
    vec3d C2 = {leaflet_end_x + t3 * (leafCurveFunc(leaflet_end_y - eps) - leaflet_end_x)/eps, leaflet_end_y - t3, leaflet_height};
    cubiccrv curve = cubiccrv(C0,C1,C2,C3);

    return curve;
}

cubiccrv Model::bendCurveBezier(){
    vec3d C0 = {bending_radius_in,0,bending_height};
    vec3d C1 = {bending_radius_in + t0 * (bendCurveFunc(eps) - bending_radius_in)/eps, 0 + t0, bending_height};
    vec3d C2 = {bending_end_x + t3 * (bendCurveFunc(bending_end_y - eps) - bending_end_x)/eps, bending_end_y - t3, bending_height};
    vec3d C3 = {bending_end_x, bending_end_y, bending_height};
    cubiccrv curve = cubiccrv(C0,C1,C2,C3);

    return curve;
}

cubiccrv Model::sinCurveBottomBezier() {
    vec3d tC0 = {root_radius, 0.0, 0.0};
    vec3d tC2 = {leaflet_end_x, leaflet_end_y, leaflet_height};
    vec3d tC3 = {bending_end_x, bending_end_y, bending_height};
    vec3d tC1 = (tC3 - pow((1 - s_sin), 2.0) * tC0 - pow(s_sin, 2.0) * tC2) / (2 * (1 - s_sin) * s_sin);

    vec3d C0 = tC0;
    vec3d C2 = tC3;
    vec3d C1 = (1-s_sin)*tC0 + s_sin*tC1;

    cubiccrv curve = cubiccrv(C0,C1,C2);
    return curve;
}

cubiccrv Model::sinCurveTopBezier() {
    vec3d tC0 = {root_radius, 0.0, 0.0};
    vec3d tC2 = {leaflet_end_x, leaflet_end_y, leaflet_height};
    vec3d tC3 = {bending_end_x, bending_end_y, bending_height};
    vec3d tC1 = (tC3 - pow((1 - s_sin), 2.0) * tC0 - pow(s_sin, 2.0) * tC2) / (2 * (1 - s_sin) * s_sin);


    vec3d C0 = tC2;
    vec3d C2 = tC3;
    vec3d C1 = (1-s_sin)*tC1 + s_sin*tC2;

    cubiccrv curve = cubiccrv(C0,C1,C2);
    return curve;
}

cubiccrv Model::symCurveBottomBezier(){
    vec3d tC0 = {root_radius,0.0,0.0};
    vec3d tC3 = {leaflet_radius_in,0.0,leaflet_height};
    vec3d tC1 = {(root_radius + bending_radius_in)/2 + xQ1_sym*(root_radius - bending_radius_in), 0.0, bending_height/2 + zQ1_sym*bending_height};
    vec3d B = {bending_radius_in,0,bending_height};
    vec3d tC2 = (B - pow((1 - s_sym), 3.0) * tC0 - 3 * s_sym * pow((1 - s_sym), 2.0) * tC1 -
            pow(s_sym, 3.0) * tC3) / (3 * (1 - s_sym) * pow(s_sym, 2.0));

    vec3d C0 = tC0;
    vec3d C3 = B;

    vec3d C1 = (1-s_sym)*tC0 + s_sym*tC1;
    vec3d R12 = (1-s_sym)*tC1 + s_sym*tC2;
    vec3d C2 = (1-s_sym)*C1 + s_sym*R12;

    cubiccrv curve = cubiccrv(C0,C1,C2,C3);

    return curve;
}

cubiccrv Model::symCurveTopBezier(){

    vec3d tC0 = {root_radius,0,0};
    vec3d tC3 = {leaflet_radius_in,0,leaflet_height};
    vec3d tC1 = {(root_radius + bending_radius_in)/2 + xQ1_sym*(root_radius - bending_radius_in), 0, bending_height/2 + zQ1_sym*bending_height};
    vec3d B = {bending_radius_in,0,bending_height};

    vec3d tC2 = (B - pow((1 - s_sym), 3) * tC0 - 3 * s_sym * pow((1 - s_sym), 2.0) * tC1 -
                  pow(s_sym, 3.0) * tC3) / (3 * (1 - s_sym) * pow(s_sym, 2.0));


    vec3d C0 = tC3;
    vec3d C3 = B;

    vec3d C1 = (1-s_sym)*tC2 + s_sym*tC3;
    vec3d R12 = (1-s_sym)*tC1 + s_sym*tC2;
    vec3d C2 = (1-s_sym)*R12 + s_sym*C1;

    cubiccrv curve = cubiccrv(C0,C1,C2,C3);
    return curve;
}

vec3d Model::fillTop(double u, double v){



    vec3d point = (1 - v) * bendCurve.f(u) + v * leafCurve.f(u) +
            (1 - u) * symCurveTop.f(1-v) + u * sinCurveTop.f(1-v) -
            ((1 - u) * (1 - v) * Q_b_sym + u * v * Q_l_sin +
            u * (1 - v) * Q_b_sin + (1 - u) * v * Q_l_sym);


    return point;
}

vec3d Model::fillTopDerU(double u, double v) {


    vec3d point = (1 - v) * bendCurve.df(u) + v * leafCurve.df(u)
             - symCurveTop.f(1-v) + sinCurveTop.f(1-v) -
               (- (1 - v) * Q_b_sym + v * Q_l_sin +
               (1 - v) * Q_b_sin - v * Q_l_sym);


    return point;
}

vec3d Model::fillTopDerUU(double u, double v){

    vec3d point = (1 - v) * bendCurve.ddf(u) + v * leafCurve.ddf(u);

    return point;
}

vec3d Model::fillTopDerUV(double u, double v) {

    vec3d point = - bendCurve.df(u) + leafCurve.df(u)
                   + symCurveTop.df(1-v) - sinCurveTop.df(1-v) -
                   (Q_b_sym + Q_l_sin - Q_b_sin - Q_l_sym);

    return point;
}


vec3d Model::fillTopDerV(double u, double v){

    vec3d point = - bendCurve.f(u) + leafCurve.f(u) -
               (1 - u) * symCurveTop.df(1-v) - u * sinCurveTop.df(1-v) -
               (- (1 - u) * Q_b_sym + u * Q_l_sin
                - u * Q_b_sin + (1 - u) * Q_l_sym);


    return point;
}

vec3d Model::fillTopDerVV(double u, double v){

    vec3d point = (1 - u) * symCurveTop.ddf(1-v) + u * sinCurveTop.ddf(1-v);

    return point;
}


vec3d Model::fillBottom(double u, double v){

    vec3d point = v * Q_r + (1 - v) * bendCurve.f(u) +
               u * sinCurveBottom.f(1-v) + (1 - u) * symCurveBottom.f(1-v) -
               (u * v * Q_r + (1 - u) * (1 - v) * Q_b_sym +
                (1 - u) * v * Q_r + u * (1 - v) * Q_b_sin);

    return point;
}

double Model::squaredTopDist(double u, double v, vec3d A){
    vec3d top = fillTop(u,v);
    double dist = top.sqdist(A);
    return dist;
}

double Model::squaredTopDistDerU(double u, double v, vec3d A){
    vec3d top = fillTop(u,v);
    vec3d topDer = fillTopDerU(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}
double Model::squaredTopDistDerV(double u, double v, vec3d A){
    vec3d top = fillTop(u,v);
    vec3d topDer = fillTopDerV(u,v);
    double dist = 2*((top - A)*topDer).sum();
    return dist;
}
double Model::squaredTopDistDerUU(double u, double v, vec3d A){
    vec3d top = fillTop(u,v);
    vec3d topDer1 = fillTopDerU(u,v);
    vec3d topDer2 = fillTopDerUU(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}
double Model::squaredTopDistDerVV(double u, double v, vec3d A){
    vec3d top = fillTop(u,v);
    vec3d topDer1 = fillTopDerV(u,v);
    vec3d topDer2 = fillTopDerVV(u,v);

    double dist = 2*((topDer1*topDer1).sum() + ((top - A)*topDer2).sum());

    return dist;
}
double Model::squaredTopDistDerUV(double u, double v, vec3d A){
    vec3d top = fillTop(u,v);
    vec3d topDerU = fillTopDerU(u,v);
    vec3d topDerV = fillTopDerV(u,v);
    vec3d topDerUV = fillTopDerUV(u,v);

    double dist = 2*(((top - A)*topDerUV).sum() + (topDerU*topDerV).sum());

    return dist;
}

double Model::distanceToTopPoint(double u, double v, vec3d A){

    return sqrt(squaredTopDist(u,v,A));
}

double Model::completeDistanceTopDerU(double u, double v, vec3d A){

    return squaredTopDistDerU(u,v,A) / distanceToTopPoint(u,v,A);
}
double Model::completeDistanceTopDerV(double u, double v, vec3d A) {

    return squaredTopDistDerV(u, v, A) / distanceToTopPoint(u, v, A);
}

double Model::completeDistanceTopDerUU(double u, double v, vec3d A) {
    vec3d f = fillTop(u,v);
    vec3d f1 = fillTopDerU(u,v);
    vec3d f2 = fillTopDerUU(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}

double Model::completeDistanceTopDerVV(double u, double v, vec3d A) {
    vec3d f = fillTop(u,v);
    vec3d f1 = fillTopDerV(u,v);
    vec3d f2 = fillTopDerVV(u,v);

    double dist = ((f - A)*f2 + f1*f1).sum()/(sqrt(((f-A)*(f-A)).sum()) - pow(((f-A)*f1).sum(),2.0)/(pow(((f-A)*(f-A)).sum(),1.5)));

    return dist;
}

double Model::completeDistanceTopDerUV(double u, double v, vec3d A) {
    vec3d f = fillTop(u, v);
    vec3d f1u = fillTopDerU(u, v);
    vec3d f1v = fillTopDerV(u, v);
    vec3d f2 = fillTopDerUV(u, v);

    double dist = ((((f - A)*(f - A)).sum())*((f2*(f-A)).sum() + (f1u*f1v).sum()) - (f1u*(f-A)).sum()*(f1v*(f-A)).sum())/(pow(((f-A)*(f-A)).sum(),1.5));

    return dist;
}

double Model::distanceToBottomPoint(double u, double v, vec3d A){
    vec3d bottom = fillBottom(u, v);
    double dist = ((bottom - A)*(bottom - A)).sum();

    return sqrt(dist);
}

Model Model::getPart(Model model, int idx){
    if(idx % 2 != 0){
        model.leafCurve = model.leafCurve.MirrorY();
        model.bendCurve = model.bendCurve.MirrorY();
        model.sinCurveBottom = model.sinCurveBottom.MirrorY();
        model.sinCurveTop = model.sinCurveTop.MirrorY();
        model.Q_b_sym = model.Q_b_sym.MirrorY();
        model.Q_l_sym = model.Q_l_sym.MirrorY();
        model.Q_b_sin = model.Q_b_sin.MirrorY();
        model.Q_l_sin = model.Q_l_sin.MirrorY();
        model.Q_r = model.Q_r.MirrorY();

        idx--;
    }
    idx /= 2;
    for(int i = 0; i < idx; i++){
        model.leafCurve = model.leafCurve.RotAxisZ(120.0/180*M_PI);
        model.bendCurve = model.bendCurve.RotAxisZ(120.0/180*M_PI);
        model.symCurveBottom = model.symCurveBottom.RotAxisZ(120.0/180*M_PI);
        model.symCurveTop = model.symCurveTop.RotAxisZ(120.0/180*M_PI);
        model.sinCurveBottom = model.sinCurveBottom.RotAxisZ(120.0/180*M_PI);
        model.sinCurveTop = model.sinCurveTop.RotAxisZ(120.0/180*M_PI);
        model.Q_b_sym = model.Q_b_sym.RotAxisZ(120.0/180*M_PI);
        model.Q_l_sym = model.Q_l_sym.RotAxisZ(120.0/180*M_PI);
        model.Q_b_sin = model.Q_b_sin.RotAxisZ(120.0/180*M_PI);
        model.Q_l_sin = model.Q_l_sin.RotAxisZ(120.0/180*M_PI);
        model.Q_r = model.Q_r.RotAxisZ(120.0/180*M_PI);
    }
    return model;
}

TopParametric Model::getTopParametric() {
    return TopParametric(
            leafCurveBezier(),
            bendCurveBezier(),
            symCurveTopBezier(),
            sinCurveTopBezier(),
            Q_l_sin,
            Q_l_sym,
            Q_b_sin,
            Q_b_sym
            );
}

BottomParametric Model::getBottomParametric() {
    return BottomParametric(
            sinCurveBottom,
            symCurveBottom,
            bendCurve,
            Q_r,
            Q_b_sin,
            Q_b_sym
    );
}

bicubicsrf Model::getTopBezier() {
    TopParametric top = getTopParametric();
    cubiccrv lC = leafCurveBezier();
    cubiccrv bC = bendCurveBezier();
    cubiccrv symC = symCurveTopBezier();
    cubiccrv sinC = sinCurveTopBezier();

    vec3d k00 = bC.getCtrlP(0);
    vec3d k01 = symC.getCtrlP(2);
    vec3d k10 = bC.getCtrlP(1);
    vec3d dUV00 = top.atDerUV(0,0);
    vec3d k11 = k01 + k10 - k00;// - dUV00/100;

    vec3d k30 = bC.getCtrlP(3);
    vec3d k31 = sinC.getCtrlP(2);
    vec3d k20 = bC.getCtrlP(2);
    vec3d dUV10 = top.atDerUV(1,0);
    vec3d k21 = k31 + k20 - k30;// - dUV10/100;

    vec3d k03 = lC.getCtrlP(0);
    vec3d k02 = symC.getCtrlP(1);
    vec3d k13 = lC.getCtrlP(1);
    vec3d dUV01 = top.atDerUV(0,1);
    vec3d k12 = k02 + k13 - k03;// - dUV01/100;

    vec3d k33 = lC.getCtrlP(3);
    vec3d k32 = sinC.getCtrlP(1);
    vec3d k23 = lC.getCtrlP(2);
    vec3d dUV11 = top.atDerUV(1,1);
    vec3d k22 = k32 + k23 - k33;// - dUV11/100;

    vec3d ctrls[16] ={k00,k01,k02,k03,
                      k10,k11,k12,k13,
                      k20,k21,k22,k23,
                      k30,k31,k32,k33};

    return bicubicsrf(ctrls);
}
