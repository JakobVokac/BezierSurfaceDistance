//
// Created by s3179222 on 10/14/19.
//

#include "bezier.h"
vector<double> cubicBezier(double t, CBezier b){
    vector<double> point = vector<double>(3);
    for(int i = 0; i<3; i++) {
        point[i] = cubicBezier(t, i, b);
    }
    return point;
}

double cubicBezier(double t, int i, CBezier b){

    return pow((1 - t), 3.0) * b.C0[i] + 3 * pow((1 - t), 2.0) * t * b.C1[i] +
           3 * (1 - t) * pow(t, 2.0) * b.C2[i] + pow(t, 3.0) * b.C3[i];

}

vector<double> cubicBezierDer1(double t, CBezier b){
    vector<double> point = vector<double>(3);
    for(int i = 0; i<3; i++) {
        point[i] = cubicBezierDer1(t, i, b);
    }
    return point;
}

double cubicBezierDer1(double t, int i, CBezier b) {

    return 3 * pow((1 - t), 2.0) * (b.C1[i] - b.C0[i]) + 6 * t * (1 - t) * (b.C2[i] - b.C1[i]) +
           3 * pow(t, 2.0) * (b.C3[i] - b.C2[i]);

}

vector<double> cubicBezierDer2(double t, CBezier b){
    vector<double> point = vector<double>(3);
    for(int i = 0; i<3; i++) {
        point[i] = cubicBezierDer2(t, i, b);
    }
    return point;
}

double cubicBezierDer2(double t, int i, CBezier b) {

    return 6 * (1 - t) * (b.C2[i] - 2 * b.C1[i] + b.C0[i]) + 6 * t * (b.C3[i] - 2 * b.C2[i] + b.C1[i]);

}

vector<double> quadBezier(double t, QBezier b){
    vector<double> point = vector<double>(3);
    for(int i = 0; i<3; i++) {
        point[i] = quadBezier(t, i, b);
    }
    return point;
}

double quadBezier(double t, int i, QBezier b){

    return pow((1 - t), 2.0) * b.C0[i] + 2 * (1 - t) * t * b.C1[i] +
           pow(t, 2.0) * b.C2[i];

}

vector<double> quadBezierDer1(double t, QBezier b){
    vector<double> point = vector<double>(3);
    for(int i = 0; i<3; i++) {
        point[i] = quadBezierDer1(t, i, b);
    }
    return point;
}

double quadBezierDer1(double t, int i, QBezier b) {

    return 2 * (1 - t) * (b.C1[i] - b.C0[i]) + 2 * t * (b.C2[i] - b.C1[i]);

}

vector<double> quadBezierDer2(double t, QBezier b){
    vector<double> point = vector<double>(3);
    for(int i = 0; i<3; i++) {
        point[i] = quadBezierDer2(t, i, b);
    }
    return point;
}

double quadBezierDer2(double t, int i, QBezier b) {

    return 2*b.C0[i] - 4*b.C1[i] + 2*b.C2[i];

}

double sqDistBezier(CBezier c, double t, vector<double> A){
    vector<double> P = cubicBezier(t,c);

    double distance = pow(P[0] - A[0], 2.0) + pow(P[1] - A[1], 2.0) + pow(P[2] - A[2], 2.0);
    return distance;
}

double sqDistBezierDer1(CBezier c, double t, vector<double> A){
    vector<double> P = cubicBezier(t,c);
    vector<double> Pd = cubicBezierDer1(t,c);
    double distance = 2*((P[0] - A[0])*Pd[0] + (P[1] - A[1])*Pd[1] + (P[2] - A[2])*Pd[2]);
    return distance;
}

double sqDistBezierDer2(CBezier c, double t, vector<double> A){
    vector<double> P = cubicBezier(t,c);
    vector<double> Pd = cubicBezierDer1(t,c);
    vector<double> Pdd = cubicBezierDer2(t,c);

    double distance = 2*(pow((Pd[0]),2.0) + (P[0] + A[0])*Pdd[0]
                         + pow((Pd[1]),2.0) + (P[1] + A[1])*Pdd[1]
                         + pow((Pd[2]),2.0) + (P[2] + A[2])*Pdd[2]);
    return distance;
}

vector<double> RotAxisZ(double rad, vector<double> P){
    vector<double> Pn = {0,0,0};

    Pn[0] = P[0] * cos(rad) - P[1] * sin(rad) + 0   ;
    Pn[1] = P[0] * sin(rad) + P[1] * cos(rad) + 0   ;
    Pn[2] = 0               + 0               + P[2];

    return Pn;
}

QBezier RotAxisZ(double rad, QBezier b){
    vector<double> C0 = {0,0,0}, C1 = {0,0,0}, C2 = {0,0,0};

    C0[0] = b.C0[0] * cos(rad) - b.C0[1] * sin(rad) + 0      ;
    C0[1] = b.C0[0] * sin(rad) + b.C0[1] * cos(rad) + 0      ;
    C0[2] = 0                  + 0                  + b.C0[2];

    C1[0] = b.C1[0] * cos(rad) - b.C1[1] * sin(rad) + 0      ;
    C1[1] = b.C1[0] * sin(rad) + b.C1[1] * cos(rad) + 0      ;
    C1[2] = 0                  + 0                  + b.C1[2];

    C2[0] = b.C2[0] * cos(rad) - b.C2[1] * sin(rad) + 0      ;
    C2[1] = b.C2[0] * sin(rad) + b.C2[1] * cos(rad) + 0      ;
    C2[2] = 0                  + 0                  + b.C2[2];

    QBezier bNew = {C0,C1,C2};
    bNew.C0 = C0;
    bNew.C1 = C1;
    bNew.C2 = C2;
    return bNew;
}

CBezier RotAxisZ(double rad, CBezier b){
    vector<double> C0 = {0,0,0}, C1 = {0,0,0}, C2 = {0,0,0}, C3 = {0,0,0};

    C0[0] = b.C0[0] * cos(rad) - b.C0[1] * sin(rad) + 0      ;
    C0[1] = b.C0[0] * sin(rad) + b.C0[1] * cos(rad) + 0      ;
    C0[2] = 0                  + 0                  + b.C0[2];

    C1[0] = b.C1[0] * cos(rad) - b.C1[1] * sin(rad) + 0      ;
    C1[1] = b.C1[0] * sin(rad) + b.C1[1] * cos(rad) + 0      ;
    C1[2] = 0                  + 0                  + b.C1[2];

    C2[0] = b.C2[0] * cos(rad) - b.C2[1] * sin(rad) + 0      ;
    C2[1] = b.C2[0] * sin(rad) + b.C2[1] * cos(rad) + 0      ;
    C2[2] = 0                  + 0                  + b.C2[2];

    C3[0] = b.C3[0] * cos(rad) - b.C3[1] * sin(rad) + 0      ;
    C3[1] = b.C3[0] * sin(rad) + b.C3[1] * cos(rad) + 0      ;
    C3[2] = 0                  + 0                  + b.C3[2];

    CBezier bNew = {C0,C1,C2,C3};
    bNew.C0 = C0;
    bNew.C1 = C1;
    bNew.C2 = C2;
    bNew.C3 = C3;
    return bNew;
}

QBezier MirrorY(QBezier b){
    b.C0[1] *= -1;
    b.C1[1] *= -1;
    b.C2[1] *= -1;
    return b;
}

CBezier MirrorY(CBezier b){
    b.C0[1] *= -1;
    b.C1[1] *= -1;
    b.C2[1] *= -1;
    b.C3[1] *= -1;
    return b;
}