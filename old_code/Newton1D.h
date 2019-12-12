//
// Created by s3179222 on 12/3/19.
//

#ifndef HEARTVALVEMODEL_NEWTON1D_H
#define HEARTVALVEMODEL_NEWTON1D_H

#include <cmath>
#include <vector>
#include "plotting.h"

using namespace std;

namespace plt = matplotlibcpp;

double func(double x);
double funcdx(double x);
double funcdxdx(double x);
double NewtonMethodTest1D(double x, double eps, int plot);

#endif //HEARTVALVEMODEL_NEWTON1D_H
