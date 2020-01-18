//
// Created by s3179222 on 12/16/19.
//
/*
 * This class is an implementation of a recursive Bezier surface splitting algorithm for finding a point's closest distance
 * to the surface.
 *
 * The algorithm takes in a Bezier surface and a point (vec3d), uses functions from the surface class to check whether the point's
 * closest distance is on the edge of the surface or inside.
 *
 * First, it takes in a surface and checks whether the closest distance to the point is on an edge. If so, the algorithm
 * continues with the edge. It first checks whether the closest distance sits on a corner. If so, it found the solution.
 * If not, it splits the edge recursively. If an edge is split, one part will have a closest distance on the inside and
 * one on a corner. The algorithm continues with the first edge. It recursively splits it until the resulting edge's is
 * length is less than the specified margin of error (eps). Any distance taken from that edge is then a good enough
 * approximation for the true distance.
 *
 * If the closest distance is to the inside of the surface, then the algorithm continually splits the surface much like
 * the edge. The main differences here are that each iteration changes the direction of the split such that each surface
 * is split into 4 equally sized patches in 2 iterations. This process is continued until the algorithm reaches a patch,
 * whose diagonal length is less than the margin of error. The center of that patch is taken as the closest point to the
 * original point.
 */

#ifndef HEARTVALVEMODEL_SPLITTINGALGORITHM_H
#define HEARTVALVEMODEL_SPLITTINGALGORITHM_H

#include "../optimizer/states.h"
#include "../geometry/surface/bicubicsrf.h"

class splittingAlgorithm {
public:
    static OptState2D optimize(bicubicsrf s, vec3d P, double eps, int recursiveLimit = 100);
    static void surfaceSplitting(int dir, double ul, double ur, double vl, double vr, double &uf, double &vf, bicubicsrf s, bicubicsrf &sf, vec3d P, double eps, int i = 0, int limit = 100);
    static void curveSplitting(double tl, double tr, double &tf, cubiccrv c, cubiccrv &cf, vec3d P, double eps, int i = 0, int limit = 100);
};


#endif //HEARTVALVEMODEL_SPLITTINGALGORITHM_H
