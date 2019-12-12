//
// Created by s3179222 on 12/12/19.
//

#ifndef HEARTVALVEMODEL_MEASUREMENTS_H
#define HEARTVALVEMODEL_MEASUREMENTS_H

#include "../model.h"
#include <chrono>
#include "../optimizer/optimizer.h"
#include "../optimizer/preprocessor/bisection.h"
#include "../optimizer/preprocessor/quadraticInterpolation.h"
#include "../optimizer/step/Geometric.h"
#include "../optimizer/step/Newton.h"

void randomPointFromEdgeOrSurfaceNormal(surface &sur, int iterations, int i, OptState2D &trueLoc, vec3d &P, double minDist,
                                        double distVariation);
void randomPointFromSurfaceNormal(surface &sur, OptState2D &trueLoc, vec3d &P, double minDist, double distVariation);
void randomPointInSpace(surface &sur, OptState2D &trueLoc, vec3d &P, double minDist, double distVariation);
void roughCurveSearch(const vec3d &P, curve *c, const vec3d &dir, double &Pdist, double &t);
void roughGridSearch(surface &sur, const vec3d &P, double dist, double &u, double &v);
void plotSurfaceDistance(const vec3d &P, surface &sur);
void TestOptimizerPerformance(optimizer &opt, int seed, int testType, double minDist, double distVar, int iterations = 1000,
                              double eps = 0.0000001, int plot = 0);

#endif //HEARTVALVEMODEL_MEASUREMENTS_H
