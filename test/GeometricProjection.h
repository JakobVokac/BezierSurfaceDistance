//
// Created by s3179222 on 12/3/19.
//

#ifndef HEARTVALVEMODEL_GEOMETRICPROJECTION_H
#define HEARTVALVEMODEL_GEOMETRICPROJECTION_H

#include <cmath>
#include <vector>
#include "../plotting.h"
#include "../vectorMath.h"

using namespace std;

namespace plt = matplotlibcpp;

double B(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, int dim);
double B1(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, int dim);
double B2(vector<double> Q0, vector<double> Q1, vector<double> Q2, int dim);
double C(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, vector<double> Q3, int dim);
double C1(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, vector<double> Q3, int dim);
double C2(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, vector<double> Q3, int dim);
void drawB(vector<double> Q0, vector<double> Q1, vector<double> Q2);
void drawC(vector<double> Q0, vector<double> Q1, vector<double> Q2, vector<double> Q3);
void drawLine(vector<double> P0, vector<double> P1);
vector<double> geometricProjectionTest();


#endif //HEARTVALVEMODEL_GEOMETRICPROJECTION_H
