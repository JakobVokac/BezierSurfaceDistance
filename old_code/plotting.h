//
// Created by s3179222 on 12/3/19.
//

#ifndef HEARTVALVEMODEL_PLOTTING_H
#define HEARTVALVEMODEL_PLOTTING_H

#include "../model.h"
#include "../matplotlibcpp.h"
#include "vectorMath.h"
#include "../model.h"
#include <vector>

namespace plt = matplotlibcpp;

void DrawDist(Model &model, const vec3d &P);
double constraintFunc(double u, double v);
void plotSurface(vector<double> P, Model m, double ul, double vl, double ur, double vr, double dU, double dV, int withConst);
void drawFillTopWithLine(Model m, vector<double> line, vector<double> axes);
void plotEdgeDistSq(curve &c, vec3d P);
void drawPart(Model &model, vector<vector<double>> &x, vector<vector<double>> &y, vector<vector<double>> &z, int top, int bottom);

#endif //HEARTVALVEMODEL_PLOTTING_H