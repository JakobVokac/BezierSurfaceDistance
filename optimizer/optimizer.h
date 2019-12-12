//
// Created by s3179222 on 12/9/19.
//

#ifndef HEARTVALVEMODEL_OPTIMIZER_H
#define HEARTVALVEMODEL_OPTIMIZER_H


#include "../geometry/surface/surface.h"
#include "states.h"
#include "preprocessor/preprocessor.h"
#include "step/step.h"
#include <cstdlib>
#include "../old_code/plotting.h"

class optimizer {
protected:
    preprocessor &prep2D;
    preprocessor &prep1D;
    step &step2D;
    step &step1D;
    surface &sur;
    vec3d P{};
    double eps;
    double distOld = 0;
    int iter = 0;
    int iterMax;
    int rangecount = 0;
    int rcMax;
    OptState2D loc{};
    OptState1D eLoc{};
    bool edgeLoc = false;
    bool cornerLoc = false;
    bool cornerSearched = false;
public:
    optimizer(
            preprocessor &prep2D,
            preprocessor &prep1D,
            step &step2D,
            step &step1D,
            surface &sur,
            vec3d P,
            double eps,
            int iterMax,
            int outOfRangeMax
    )
            : prep2D(prep2D), prep1D(prep1D), step2D(step2D), step1D(step1D), sur(sur) {
        this->eps = eps;
        this->iterMax = iterMax;
        this->rcMax = outOfRangeMax;
        this->sur = sur;
        this->P = P;
        this->loc = {0.5, 0.5, distOld};
    };

    OptState2D optimize();

    OptState2D optimizeForPoint(vec3d P);

    OptState1D optimizeEdge(curve &crv);

    ~optimizer() = default;

    void checkEdges();

    void plot();

    void plotEdge(curve &crv);

    void cornerSearch(double u, double v);

    void checkCorner();

    void initialize();

    surface &getSurface();

    int getIterations();

    bool edgeSolution();

    bool cornerSolution();

    bool cornerSearchSolution();
};

#endif //HEARTVALVEMODEL_OPTIMIZER_H
