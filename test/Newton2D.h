//
// Created by s3179222 on 12/3/19.
//

#ifndef HEARTVALVEMODEL_NEWTON2D_H
#define HEARTVALVEMODEL_NEWTON2D_H

#include <cmath>
#include <vector>
#include "../plotting.h"

using namespace std;

namespace plt = matplotlibcpp;

double f(double x, double y);
double fdx(double x, double y);
double fdy(double x, double y);
double fdxdx(double x, double y);
double fdydy(double x, double y);
double fdxdy(double x, double y);
double NewtonMethodTest2D(double x, double y, double eps, double sigma, int plot);


#endif //HEARTVALVEMODEL_NEWTON2D_H
