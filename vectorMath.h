//
// Created by s3179222 on 12/3/19.
//

#ifndef HEARTVALVEMODEL_VECTORMATH_H
#define HEARTVALVEMODEL_VECTORMATH_H

#include <vector>
#include <cmath>
#include <assert.h>

using namespace std;

vector<double> add(vector<double> P, vector<double> Q);
vector<double> sub(vector<double> P, vector<double> Q);
vector<double> mul(vector<double> P, double q);
vector<double> div(vector<double> P, double q);
double det(vector<double> P, vector<double> Q);
double sign2D(vector<double> P, vector<double> Q);
double magnitude(vector<double> v);
double dot(vector<double> a, vector<double> b);
vector<double> cross(vector<double> a, vector<double> b);
double curvature3D(vector<double> c1, vector<double> c2);
vector<double> circleCenterDir3D(vector<double> c1, vector<double> c2);
vector<double> solve3Dlinear(vector<double> x, vector<double> y, vector<double> z);
double sign3D(vector<double> P, vector<double> Q);
vector<double> circleCenterDir2D(vector<double> c1, vector<double> c2);
double curvature2D(vector<double> c1, vector<double> c2);

#endif //HEARTVALVEMODEL_VECTORMATH_H
