//
// Created by s3179222 on 10/14/19.
//

#ifndef HEARTVALVEMODEL_BEZIER_H
#define HEARTVALVEMODEL_BEZIER_H

#include <vector>
#include <cmath>

using namespace std;

struct QBezier{
    vector<double> C0 = vector<double>(3);
    vector<double> C1 = vector<double>(3);
    vector<double> C2 = vector<double>(3);
};

struct CBezier{
    vector<double> C0 = vector<double>(3);
    vector<double> C1 = vector<double>(3);
    vector<double> C2 = vector<double>(3);
    vector<double> C3 = vector<double>(3);
};

vector<double> cubicBezier(double t, CBezier b);
double cubicBezier(double t, int i, CBezier b);

vector<double> cubicBezierDer1(double t, CBezier b);
double cubicBezierDer1(double t, int i, CBezier b);

vector<double> cubicBezierDer2(double t, CBezier b);
double cubicBezierDer2(double t, int i, CBezier b);

vector<double> quadBezier(double t, QBezier b);
double quadBezier(double t, int i, QBezier b);

vector<double> quadBezierDer1(double t, QBezier b);
double quadBezierDer1(double t, int i, QBezier b);

vector<double> quadBezierDer2(double t, QBezier b);
double quadBezierDer2(double t, int i, QBezier b);

QBezier RotAxisZ(double rad, QBezier b);
CBezier RotAxisZ(double rad, CBezier b);
QBezier MirrorY(QBezier b);
CBezier MirrorY(CBezier b);
vector<double> RotAxisZ(double rad, vector<double> P);

#endif //HEARTVALVEMODEL_BEZIER_H
