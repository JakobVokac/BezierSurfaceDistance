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
#include "optimizer.h"

using namespace std;

OptState2D optimizer::optimize() {
    initialize();
    tLoc = loc = prep2D.preprocess(*sur, P);

    do{
        distOld = loc.dist;
        loc = step2D.doStep(*sur,P,loc);
        if(!edgesSearched) {
            checkEdges();
        }
        iter++;


    }while(abs(distOld-loc.dist) > eps && iter < iterMax);

    if(loc.u < 0 || loc.v < 0 || loc.u > 1 || loc.v > 1)
        return tLoc;

    if(tLoc.dist < loc.dist)
        return tLoc;

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

    }while(abs(distOld - eLoc.dist) > eps && iter < iterMax);

    OptState1D lc, rc;
    lc = {0,crv.distTo(0,P)};
    rc = {1,crv.distTo(1,P)};
    if(eLoc.t < 0 || eLoc.t > 1)
        return (lc.dist < rc.dist ? lc : rc);

    if(lc.dist < eLoc.dist)
        eLoc = lc;
    if(rc.dist < eLoc.dist)
        eLoc = rc;

    return eLoc;
}

void optimizer::checkEdges() {

    if (loc.u < 0 && rangecount > rcMax) {

        if (loc.v < 0) {
            cornerSearch(0, 0);
        }else if (loc.v > 1) {
            cornerSearch(0, 1);
        }else {
            OptState1D res = optimizeEdge(sur->edgeU0());

            tLoc = {0, res.t, res.dist};

            plotEdge(sur->edgeU0());
        }

        edgeLoc = true;
        edgesSearched = true;

    } else
        rangecount++;

    if (loc.u > 1 && rangecount > rcMax) {

        if(loc.v < 0) {
            cornerSearch(1, 0);
        }else if (loc.v > 1) {
            cornerSearch(1, 1);
        }else {
            OptState1D res = optimizeEdge(sur->edgeU1());

            tLoc = {1, res.t, res.dist};

            plotEdge(sur->edgeU1());
        }
        edgeLoc = true;
        edgesSearched = true;

    } else
        rangecount++;

    if (loc.v < 0 && rangecount > rcMax) {

        OptState1D res = optimizeEdge(sur->edgeV0());

        tLoc = {res.t, 0, res.dist};

        plotEdge(sur->edgeV0());

        edgeLoc = true;
        edgesSearched = true;

    } else
        rangecount++;

    if (loc.v > 1 && rangecount > rcMax) {

        OptState1D res = optimizeEdge(sur->edgeV1());

        tLoc = {res.t, 1, res.dist};

        plotEdge(sur->edgeV1());

        edgeLoc = true;
        edgesSearched = true;

    } else
        rangecount++;
}

void optimizer::cornerSearch(double u, double v) {
    OptState1D res1{}, res2{};

    if(u == 0)
        res1 = optimizeEdge(sur->edgeU0());
    else
        res1 = optimizeEdge(sur->edgeU1());
    if(v == 0)
        res2 = optimizeEdge(sur->edgeV0());
    else
        res2 = optimizeEdge(sur->edgeV1());

    if (res1.dist < res2.dist)
        tLoc = {res1.t, v, res1.dist};
    else
        tLoc = {u, res2.t, res2.dist};

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
    edgesSearched = false;
    iter = 0;
    rangecount = 0;
}

OptState2D optimizer::optimizeForPoint(vec3d P) {
    this->P = P;

    return optimize();
}

surface &optimizer::getSurface() {
    return *sur;
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
