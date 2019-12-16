//
// Created by s3179222 on 12/16/19.
//

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
