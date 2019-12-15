//
// Created by s3179222 on 12/12/19.
//

#include "measurements.h"



void TestOptimizerPerformance(optimizer &opt, int seed, int testType, double minDist, double distVar, int iterations, double eps, int plot){
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

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration;


    for (int i = 0; i < iterations; i++) {
        vec3d P;
        if(testType == 0)
            randomPointFromEdgeOrSurfaceNormal(opt.getSurface(), iterations, i, trueLoc, P, minDist, distVar);
        else if(testType == 1)
            randomPointFromSurfaceNormal(opt.getSurface(), trueLoc, P, minDist, distVar);
        else{
            randomPointInSpace(opt.getSurface(), trueLoc, P, minDist, distVar);
            cout << "Random points test: distances are approximated by a grid search; not useful for precise measurements";
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

    cout << "reliability: " << double(reliability)/iterations << " average error: " << avgError/iterations << " average iterations: " << double(sumIters)/iterations << " average time: " << avgTime/iterations << " microseconds " << endl;
    cout << "% onedimcalls: " << double(onedimcalls)/iterations << " error by 1D search %: " << double(count1D)/onedimcalls << endl;
    cout << "% cornercalls: " << double(cornercalls)/iterations << " error by corner search %: " << double(countCor)/onedimcalls << endl;

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

void randomPointInSpace(surface &sur, OptState2D &trueLoc, vec3d &P, double minDist,
                                  double distVariation) {
    double x = (double(rand()) / RAND_MAX) * distVariation + minDist;
    double y = (double(rand()) / RAND_MAX) * distVariation + minDist;
    double z = (double(rand()) / RAND_MAX) * distVariation + minDist;

    double u = double(rand())/RAND_MAX, v = double(rand())/RAND_MAX;
    P = sur.at(u,v);
    P += {x,y,z};
    double dist = x+y+z;
    roughGridSearch(sur, P, dist, u, v);
    trueLoc = {u,v,x+y+z};
}

void roughGridSearch(surface &sur, const vec3d &P, double dist, double &u, double &v) {
    for(double i = 0; i <= 1; i += 0.01){
        for(double j = 0; j <= 1; j += 0.01){
            double tDist = sur.distTo(i,j,P);
            if(tDist < dist){
                u = i;
                v = j;
                dist = tDist;
            }
        }
    }
}

void randomPointFromSurfaceNormal(surface &sur, OptState2D &trueLoc, vec3d &P, double minDist,
                               double distVariation) {
    double Pdist = (double(rand()) / RAND_MAX) * distVariation + minDist;
    double u = double(rand())/RAND_MAX, v = double(rand())/RAND_MAX;
    P = sur.at(u,v);
    vec3d Pu = sur.atDerU(u,v);
    vec3d Pv = sur.atDerV(u,v);
    vec3d N = Pu.cross(Pv);
    N /= N.mag();
    P += N*Pdist;
    trueLoc = {u,v,Pdist};
}


void randomPointFromEdgeOrSurfaceNormal(surface &sur, int iterations, int i, OptState2D &trueLoc, vec3d &P, double minDist,
                                        double distVariation) {
    double Pdist= (double(rand()) / RAND_MAX) * distVariation + minDist;
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
            randomPointFromSurfaceNormal(sur,trueLoc,P,minDist,distVariation);
            return;
        default:
            break;
    }

    dir = c->curvePlaneNormal();

    t = double(rand())/RAND_MAX;
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
        cout << "true dist for edge can't be found analytically; true distance for this point is now a rough estimate" << endl;
        for (double j = 0; j <= 1; j+=0.001) {
            double tDist = c->distTo(j,P);
            if(tDist < Pdist){
                t = j;
                Pdist = tDist;
            }
        }
    }
}

