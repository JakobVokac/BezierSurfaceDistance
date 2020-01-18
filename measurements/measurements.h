//
// Created by s3179222 on 12/12/19.
//
/*
 * This class implements different types of measurements for the optimizer. The 3 main types of measurements all take
 * a collection of randomly produced points (with some specifications) for which the distance is known and tests the
 * optimizer on them. They keep track of how many times the optimizer found a solution with a lower margin of error than
 * the specified eps (its reliability) and measures by how much it missed the solution when it does (average error).
 * They also keep track of how many times a wrong solution was produced by an edge search (the optimizer wrongly decided
 * that the closest point is on an edge) or corner search (... on a corner).
 *
 * The 3 types of tests differ by how they produce points:
 *      - the first takes a random point on the surface (not on an edge) and calculates the surface normal. It then
 *        decides on a distance and calculates a point in the direction of the surface normal for that distance.
 *        The distance is then known and correct, since the normal represents the points projection, which is the closest
 *        distance to the surface (keeping in mind that this works since the surface's curvature is too small for any better
 *        solutions on another point or in other words there is only one minimum in the distance function in the domain of
 *        the search for any point),
 *      - the second takes points from the first method, but also adds points from the same kind of projection from edges,
 *        such that the algorithm must find both surface and edge solution,
 *      - the third takes a collection of random points in space and makes a rough grid search to find a rough estimate
 *        of the solution (this test is only useful to check whether the optimizer is heavily unreliable)
 *
 * The are also tests for the splitting algorithm in the same fashion and also for detecting edge solutions.
 */
#ifndef HEARTVALVEMODEL_MEASUREMENTS_H
#define HEARTVALVEMODEL_MEASUREMENTS_H

#include "../model.h"
#include <chrono>
#include "../optimizer/optimizer.h"
#include "../optimizer/preprocessor/bisection.h"
#include "../optimizer/preprocessor/quadraticInterpolation.h"
#include "../optimizer/step/Geometric.h"
#include "../optimizer/step/Newton.h"
#include "../splittingAlgorithm/splittingAlgorithm.h"
#include "../geometry/surface/compositeBicubicsrf.h"

using namespace std;

void TestEdgeDistTestValidity(surface &sur, double minDist, double distVar, int iterations);

void TestEdgeSolutionDetection(TopParametric &sur, compositeBicubicsrf &bez, int seed, double minDist, double distVar, int iterations);

void
randomPointFromEdgeOrSurfaceNormal(surface &sur, int iterations, int i, OptState2D &trueLoc, vec3d &P, double minDist,
                                   double distVariation, default_random_engine &generator,
                                   uniform_real_distribution<double> &distribution);
void randomPointFromSurfaceNormal(surface &sur, OptState2D &trueLoc, vec3d &P, double minDist, double distVariation,
                                  default_random_engine &generator, uniform_real_distribution<double> &distribution);
void randomPointInSpace(surface &sur, OptState2D &trueLoc, vec3d &P, double minDist, double distVariation,
                        default_random_engine &generator, uniform_real_distribution<double> &distribution);
void roughCurveSearch(const vec3d &P, curve *c, const vec3d &dir, double &Pdist, double &t);
void roughGridSearch(surface &sur, const vec3d &P, double &dist, double &u, double &v);
void plotSurfaceDistance(const vec3d &P, surface &sur);
void TestOptimizerPerformance(optimizer &opt, int seed, int testType, double minDist, double distVar, int iterations = 1000,
                              double eps = 0.0000001, int plot = 0);
void TestSplitterPerformance(bicubicsrf &sur, int seed, int testType, double minDist, double distVar, int iterations, double eps, int plot);

#endif //HEARTVALVEMODEL_MEASUREMENTS_H
