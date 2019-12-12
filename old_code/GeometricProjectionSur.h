//
// Created by s3179222 on 12/3/19.
//

#ifndef HEARTVALVEMODEL_GEOMETRICPROJECTIONSUR_H
#define HEARTVALVEMODEL_GEOMETRICPROJECTIONSUR_H

#include <cmath>
#include <vector>
#include "plotting.h"
#include "vectorMath.h"

using namespace std;

namespace plt = matplotlibcpp;

int factorial(int n);
int binom(int n, int i);
vector<double> BSur(double u, double v, vector<vector<double>> k, int n, int m);
vector<double> BSurDu(double u, double v, vector<vector<double>> k, int n, int m);
vector<double> BSurDv(double u, double v, vector<vector<double>> k, int n, int m);
vector<double> BSurDuu(double u, double v, vector<vector<double>> k, int n, int m);
vector<double> BSurDuv(double u, double v, vector<vector<double>> k, int n, int m);
vector<double> BSurDvv(double u, double v, vector<vector<double>> k, int n, int m);
vector<double> BSurParam(double u0, double v0, double t, vector<vector<double>> k, int n, int m, double lamu, double lamv);
vector<double> BSurParamDt(double u0, double v0, double t, vector<vector<double>> k, int n, int m, double lamu, double lamv);
vector<double> BSurParamDtt(double u0, double v0, double t, vector<vector<double>> k, int n, int m, double lamu, double lamv);
void drawSur(vector<vector<double>> C, int n, int m);
void drawSurWithLine(vector<vector<double>> C, int n, int m, vector<double> line);
vector<double> orthogonalProjectionTestSur3D();


#endif //HEARTVALVEMODEL_GEOMETRICPROJECTIONSUR_H
