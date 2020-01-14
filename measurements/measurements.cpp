//
// Created by s3179222 on 12/12/19.
//

#include "measurements.h"


void TestEdgeDistTestValidity(surface &sur, double minDist, double distVar, int iterations){
    vec3d P;
    OptState2D trueLoc;
    OptState2D gridLoc;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    bool edgeLoc;
    int reliability = 0;

    for (int i = 0; i < iterations; i++) {
        randomPointFromEdgeOrSurfaceNormal(sur, iterations, i, trueLoc, P, minDist, distVar, generator,
                                           distribution);

        roughGridSearch(sur,P,trueLoc.dist,gridLoc.u,gridLoc.v);
        gridLoc.dist = sur.distTo(gridLoc.u,gridLoc.v,P);

        if(trueLoc.dist < gridLoc.dist && i/(iterations/5) != 4)
            reliability++;

        cout << i << endl;
    }
    plotSurface(sur);
    cout << "Edge test reliability: " << double(reliability)/((iterations*4)/5) << endl;
}

bool onEdge(OptState2D s){
    if(s.u < 1.0001 && s.u > 0.9999)
        return true;
    if(s.v < 1.0001 && s.v > 0.9999)
        return true;
    if(s.u < 0.0001 && s.u > -0.0001)
        return true;
    if(s.v < 0.0001 && s.v > -0.0001)
        return true;
    return false;
}

void TestEdgeSolutionDetection(TopParametric &sur, compositeBicubicsrf &bez, int seed, double minDist, double distVar, int iterations) {


    vec3d P;
    OptState2D trueLoc;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    bool edgeLoc;
    int reliability = 0;

    for (int i = 0; i < iterations; i++) {
        randomPointInSpace(sur,trueLoc,P,  minDist, distVar, generator,
                                           distribution);



        if(!bez.closestPointInPatch(P) == onEdge(trueLoc))
            reliability++;

        if(i % 20 == 0){
            cout << "iteration: " << i << endl;
            cout << "point: " << P << endl;
            cout << "edgeLoc: " << onEdge(trueLoc) << " u: " << trueLoc.u << " v: " << trueLoc.v << endl;
            cout << "edge solution: " << !bez.closestPointInPatch(P) << endl;
        }
    }
    plotSurface(sur);
    cout << "Edge detection from Bezier surface reliability: " << double(reliability)/iterations << endl;
}

void TestSplitterPerformance(bicubicsrf &sur, int seed, int testType, double minDist, double distVar, int iterations, double eps, int plot){
    cout << "Testing: ";
    int reliability = 0;
    double avgError = 0;
    int sumIters = 0;
    srand(seed);
    time_t avgTime = 0;
    OptState2D trueLoc;
    OptState2D loc;
    int onedimcalls = 0;
    int count1D = 0;
    int cornercalls = 0;
    int countCor = 0;
    bool edgeSolution, cornerSolution;

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for (int i = 0; i < iterations; i++) {
        cout << "iteration: " << i << endl;
        edgeSolution = false, cornerSolution = false;
        vec3d P;
        if(testType == 0)
            randomPointFromEdgeOrSurfaceNormal(sur, iterations, i, trueLoc, P, minDist, distVar, generator,
                                               distribution);
        else if(testType == 1)
            randomPointFromSurfaceNormal(sur, trueLoc, P, minDist, distVar, generator, distribution);
        else{
            randomPointInSpace(sur, trueLoc, P, minDist, distVar, generator, distribution);
            cout << "Random points test: distances are approximated by a grid search; not useful for precise measurements";
        }

        t_start = chrono::high_resolution_clock::now();
        loc = splittingAlgorithm::optimize(sur,P,eps);
        t_stop = chrono::high_resolution_clock::now();

        if(loc.u == 0 || loc.u == 1){
            if(loc.v == 0 || loc.v == 1){
                cornerSolution = true;
            }
            edgeSolution = true;
        }
        if(loc.v == 0 || loc.v == 1){
            edgeSolution = true;
        }
        t_duration = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
        avgTime += t_duration.count();

        onedimcalls += edgeSolution;
        cornercalls += cornerSolution;

        if(plot)
            cout << "u: " << loc.u << " v: " << loc.v << " V: " << loc.dist << endl;

        if(loc.dist <= trueLoc.dist + eps) {
            reliability++;
        }else{
            avgError += (loc.dist - trueLoc.dist)/trueLoc.dist;
            if(edgeSolution){
                count1D++;
            }
            if(cornerSolution){
                countCor++;
            }
            if(plot) {
                cout << "Point: " << P.getx() << " " << P.gety() << " " << P.getz() << endl;
                cout << "True location: u: " << trueLoc.u << " v: " << trueLoc.v << " True distance: " << trueLoc.dist << endl;

                plotSurfaceDistance(P, sur);
            }
        }
    }

    cout << "reliability: " << double(reliability)/iterations << " average error: " << avgError/iterations << " average time: " << avgTime/iterations << " microseconds " << endl;
    cout << "% onedimcalls: " << double(onedimcalls)/iterations << " error by 1D search %: " << double(count1D)/onedimcalls << endl;
    cout << "% cornercalls: " << double(cornercalls)/iterations << " error by corner search %: " << double(countCor)/onedimcalls << endl;

}

void TestOptimizerPerformance(optimizer &opt, int seed, int testType, double minDist, double distVar, int iterations, double eps, int plot){
    cout << "Testing: ";
    int reliability = 0;
    double avgError = 0;
    double totalDistance = 0;
    double trueTotalDist = 0;
    int sumIters = 0;
    srand(seed);
    time_t avgTime = 0;
    OptState2D trueLoc;
    OptState2D loc;
    int onedimcalls = 0;
    int count1D = 0;
    int cornercalls = 0;
    int countCor = 0;

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    if(testType == 2){
        cout << "Random points test: distances are approximated by a grid search; not useful for precise measurements" << endl;
    }
    for (int i = 0; i < iterations; i++) {
        vec3d P;
        if(testType == 0)
            randomPointFromEdgeOrSurfaceNormal(opt.getSurface(), iterations, i, trueLoc, P, minDist, distVar,
                                               generator, distribution);
        else if(testType == 1)
            randomPointFromSurfaceNormal(opt.getSurface(), trueLoc, P, minDist, distVar, generator, distribution);
        else{
            randomPointInSpace(opt.getSurface(), trueLoc, P, minDist, distVar, generator, distribution);
        }

        t_start = chrono::high_resolution_clock::now();
        loc = opt.optimizeForPoint(P);
        t_stop = chrono::high_resolution_clock::now();

        t_duration = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
        avgTime += t_duration.count();
        onedimcalls += opt.edgeSolution();
        cornercalls += opt.cornerSearchSolution();

        if(plot)
            cout << "u: " << loc.u << " v: " << loc.v << " V: " << loc.dist << " iterations: " << opt.getIterations() << endl;

        if(isnan(loc.dist)){
            cout << "Iteration: " << i << " u: " << loc.u << " v: " << loc.v << endl;
        }
        
        totalDistance += loc.dist;
        trueTotalDist += trueLoc.dist;

        if(loc.dist <= trueLoc.dist + eps) {
            reliability++;
        }else{
            avgError += (loc.dist - trueLoc.dist)/trueLoc.dist;
            if(opt.edgeSolution()){
                count1D++;
            }
            if(opt.cornerSearchSolution()){
                countCor++;
            }
            if(plot) {
                cout << "Point: " << P.getx() << " " << P.gety() << " " << P.getz() << endl;
                cout << "True location: u: " << trueLoc.u << " v: " << trueLoc.v << " True distance: " << trueLoc.dist << endl;
                surface &sur = opt.getSurface();

                plotSurfaceDistance(P, sur);
            }
        }

        sumIters += opt.getIterations();

    }
    cout << "total error: " << avgError << endl;
    cout << "reliability: " << double(reliability)/iterations << " average error: " << avgError/iterations << " average iterations: " << double(sumIters)/iterations << " average time: " << avgTime/iterations << " microseconds " << endl;
    cout << "% onedimcalls: " << double(onedimcalls)/iterations << " error by 1D search %: " << double(count1D)/onedimcalls << endl;
    cout << "% cornercalls: " << double(cornercalls)/iterations << " error by corner search %: " << double(countCor)/onedimcalls << endl;
    cout << "Total distance: " << totalDistance << ", true total distance: " << trueTotalDist << endl;
}

void plotSurfaceDistance(const vec3d &P, surface &sur) {
    vector<vector<double>> u, v, dist;
    double min_dist = 100000;
    for (double i = 0; i <= 1; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j <= 1; j += 0.05) {
            uRow.push_back(j);
            vRow.push_back(i);
            double tDist = sur.sqDistTo(j, i, P);
            distRow.push_back(tDist);
            if (min_dist > tDist) {
                min_dist = tDist;
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    vector<double> x2, y2, z2, u2, v2, w2;
    for (double i = 0; i <= 1; i += 0.05) {
        for (double j = 0; j <= 1; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = sur.sqDistToDerU(j, i, P);
            tV = sur.sqDistToDerV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);
    plt::show();
}

void plotSurface(const vec3d &P, surface &sur, int derivative = 0) {
    vector<vector<double>> u, v, dist;
    double min_dist = 100000;
    for (double i = 0; i <= 1; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j <= 1; j += 0.05) {
            vec3d p{};
            switch (derivative){
                case 0:
                    p = sur.at(j,i);
                    break;
                case 1:
                    p = sur.atDerU(j,i);
                    break;
                case 2:
                    p = sur.atDerV(j,i);
                    break;
                case 3:
                    p = sur.atDerUU(j,i);
                    break;
                case 4:
                    p = sur.atDerVV(j,i);
                    break;
                case 5:
                    p = sur.atDerUV(j,i);
                    break;
                default:
                    break;
            }
            uRow.push_back(p.getx());
            vRow.push_back(p.gety());
            distRow.push_back(p.getz());
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
}

void randomPointInSpace(surface &sur, OptState2D &trueLoc, vec3d &P, double minDist, double distVariation,
                        default_random_engine &generator, uniform_real_distribution<double> &distribution) {

    double x = distribution(generator) * distVariation + minDist;
    double y = distribution(generator) * distVariation + minDist;
    double z = distribution(generator) * distVariation + minDist;

    double u = distribution(generator), v = distribution(generator);
    P = sur.at(u,v);
    P += {x,y,z};
    double dist = x+y+z;
    roughGridSearch(sur, P, dist, u, v);
    trueLoc = {u,v,dist};
}

void roughGridSearch(surface &sur, const vec3d &P, double &dist, double &u, double &v) {
    for(double i = 0; i <= 1.000001; i += 0.01){
        for(double j = 0; j <= 1.000001; j += 0.01){
            double tDist = sur.distTo(i,j,P);
            if(tDist < dist){
                u = i;
                v = j;
                dist = tDist;
            }
        }
    }
}

void randomPointFromSurfaceNormal(surface &sur, OptState2D &trueLoc, vec3d &P, double minDist, double distVariation,
                                  default_random_engine &generator, uniform_real_distribution<double> &distribution) {

    double Pdist = distribution(generator) * distVariation + minDist;
    double u = distribution(generator), v = distribution(generator);
    P = sur.at(u,v);
    vec3d Pu = sur.atDerU(u,v);
    vec3d Pv = sur.atDerV(u,v);
    vec3d N = Pu.cross(Pv);
    N /= N.mag();
    P += N*Pdist;
    trueLoc = {u,v,Pdist};
}


void
randomPointFromEdgeOrSurfaceNormal(surface &sur, int iterations, int i, OptState2D &trueLoc, vec3d &P, double minDist,
                                   double distVariation, default_random_engine &generator,
                                   uniform_real_distribution<double> &distribution) {

    double Pdist= distribution(generator) * distVariation + minDist;
    curve *c;
    vec3d dir{};
    double t;
    vec3d P1{}, P2{};

    switch (i/(iterations/5)){
        case 0:
            c = &sur.edgeV0();
            break;
        case 1:
            c = &sur.edgeV1();
            break;
        case 2:
            c = &sur.edgeU0();
            break;
        case 3:
            c = &sur.edgeU1();
            break;
        case 4:
            randomPointFromSurfaceNormal(sur, trueLoc, P, minDist, distVariation, generator, distribution);
            return;
        default:
            break;
    }

    dir = c->curvePlaneNormal();
    t = distribution(generator);
    P1 = c->f(t) + dir * Pdist;
    P2 = c->f(t) - dir * Pdist;

    switch (i/(iterations/5)){
        case 0:
            if(sur.distTo(t,0.1,P1) > sur.distTo(t,0.1,P2))
                P = P1;
            else
                P = P2;

            roughCurveSearch(P, c, dir, Pdist, t);

            trueLoc = {t,0,Pdist};
            break;
        case 1:
            if(sur.distTo(t,0.9,P1) > sur.distTo(t,0.9,P2))
                P = P1;
            else
                P = P2;

            roughCurveSearch(P, c, dir, Pdist, t);

            trueLoc = {t,1,Pdist};
            break;
        case 2:
            if(sur.distTo(0.1,t,P1) > sur.distTo(0.1,t,P2))
                P = P1;
            else
                P = P2;

            roughCurveSearch(P, c, dir, Pdist, t);

            trueLoc = {0,t,Pdist};
            break;
        case 3:
            if(sur.distTo(0.9,t,P1) > sur.distTo(0.9,t,P2))
                P = P1;
            else
                P = P2;

            roughCurveSearch(P, c, dir, Pdist, t);

            trueLoc = {1,t,Pdist};
            break;
        default:
            break;
    }
}

void roughCurveSearch(const vec3d &P, curve *c, const vec3d &dir, double &Pdist, double &t) {
    if(dir == vec3d{0, 0, 0}){
//        cout << "true dist for edge can't be found analytically; true distance for this point is now a rough estimate" << endl;
        for (double j = 0; j <= 1; j+=0.001) {
            double tDist = c->distTo(j,P);
            if(tDist < Pdist){
                t = j;
                Pdist = tDist;
            }
        }
    }
}


