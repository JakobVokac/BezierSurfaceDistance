//
// Created by s3179222 on 12/9/19.
//
/*
 * This is the main numerical optimizer class. It implements the skeleton of the numerical optimization, which includes:
 *  - calling the preprocessing method for the surface search (prep2D)
 *  - iterating through the surface search with a specified step (step2D),
 *  - checking if the optimization has been outsize of the domain ((u,v) e [0,1]^2) for 2 or more iterations (rangecount)
 *    in which case an edge or corner search is called,
 *  - calling the preprocessing method for the edge search (prep1D)
 *  - iterating through the edge with a specified step (step1D)
 *  - choosing the most appropirate solution between the surface and edge searches,
 *  - keeping track of the number of iterations
 */
#ifndef HEARTVALVEMODEL_OPTIMIZER_H
#define HEARTVALVEMODEL_OPTIMIZER_H


#include "../geometry/surface/surface.h"
#include "states.h"
#include "preprocessor/preprocessor.h"
#include "step/step.h"
#include <cstdlib>
#include "../plotting/plotting.h"

class optimizer {
protected:
    preprocessor &prep2D;
    preprocessor &prep1D;
    step &step2D;
    step &step1D;
    surface *sur;
    vec3d P{};
    double eps;
    double distOld = 0;
    int iter = 0;
    int iterMax;
    int rangecount = 0;
    int rcMax;
    OptState2D loc{};
    OptState1D eLoc{};
    OptState2D tLoc{};
    bool edgeLoc = false;
    bool cornerLoc = false;
    bool cornerSearched = false;
    bool edgesSearched = false;
public:
    optimizer(
            preprocessor &prep2D,
            preprocessor &prep1D,
            step &step2D,
            step &step1D,
            surface *sur,
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

    void setSur(surface *sur){
        this->sur = sur;
    }
};

#endif //HEARTVALVEMODEL_OPTIMIZER_H
