//
// Created by s3179222 on 12/9/19.
//

#include "optimizer.h"

using namespace std;

OptState2D optimizer::optimize() {
    initialize();
    loc = prep2D.preprocess(sur, P);
    do{
        distOld = loc.dist;
        loc = step2D.doStep(sur,P,loc);
        checkEdges();
        iter++;
    }while(abs(distOld-loc.dist) > eps && iter < iterMax && !edgeLoc);

    return loc;
}

OptState1D optimizer::optimizeEdge(curve &crv) {
    eLoc = prep1D.preprocess(crv, P);
    eLoc.dist = crv.distTo(eLoc.t,P);
    do{
        distOld = eLoc.dist;
        eLoc = step1D.doStep(crv, P, eLoc);
        checkCorner();
        iter++;
    }while(abs(distOld - eLoc.dist) > eps && iter < iterMax && !cornerLoc);

    return eLoc;
}

void optimizer::checkEdges() {

    if (loc.u < 0 && rangecount > rcMax) {

        if (loc.v < 0) {
            cornerSearch(0, 0);
        }else if (loc.v > 1) {
            cornerSearch(0, 1);
        }else {
            OptState1D res = optimizeEdge(sur.edgeU0());
            loc = {0, res.t, res.dist};

            plotEdge(sur.edgeU0());
        }

        edgeLoc = true;

    } else
        rangecount++;

    if (loc.u > 1 && rangecount > rcMax) {

        if(loc.v < 0) {
            cornerSearch(1, 0);
        }else if (loc.v > 1) {
            cornerSearch(1, 1);
        }else {
            OptState1D res = optimizeEdge(sur.edgeU1());
            loc = {1, res.t, res.dist};

            plotEdge(sur.edgeU1());
        }
        edgeLoc = true;

    } else
        rangecount++;

    if (loc.v < 0 && rangecount > rcMax) {

        OptState1D res = optimizeEdge(sur.edgeV0());
        loc = {res.t, 0, res.dist};

        plotEdge(sur.edgeV0());

        edgeLoc = true;

    } else
        rangecount++;

    if (loc.v > 1 && rangecount > rcMax) {

        OptState1D res = optimizeEdge(sur.edgeV1());
        loc = {res.t, 1, res.dist};

        plotEdge(sur.edgeV1());

        edgeLoc = true;

    } else
        rangecount++;
}

void optimizer::cornerSearch(double u, double v) {
    OptState1D res1{}, res2{};

    res1 = optimizeEdge(sur.edgeU0());
    res2 = optimizeEdge(sur.edgeV0());

    plotEdge(sur.edgeU0());
    plotEdge(sur.edgeV0());

    if (res1.dist < res2.dist)
        loc = {res1.t, v, res1.dist};
    else
        loc = {u, res2.t, res2.dist};

    edgeLoc = true;
    cornerSearched = true;
}

void optimizer::plotEdge(curve &crv) {
//    cout << "loc: u: " << loc.u << " v: " << loc.v << endl;
//    plotEdgeDistSq(crv, P);
}

void optimizer::plot() {}

void optimizer::checkCorner() {
    if(eLoc.t < 0){
        eLoc.t = 0;
        cornerLoc = true;
    }else if(eLoc.t > 1){
        eLoc.t = 1;
        cornerLoc = true;
    }
}

void optimizer::initialize() {
    edgeLoc = false;
    cornerLoc = false;
    cornerSearched = false;
    iter = 0;
    rangecount = 0;
}

OptState2D optimizer::optimizeForPoint(vec3d P) {
    this->P = P;

    return optimize();
}

surface &optimizer::getSurface() {
    return sur;
}

bool optimizer::edgeSolution() {
    return edgeLoc;
}

bool optimizer::cornerSearchSolution() {
    return cornerSearched;
}


bool optimizer::cornerSolution() {
    return cornerLoc;
}

int optimizer::getIterations() {
    return iter;
}
