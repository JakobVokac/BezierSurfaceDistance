#include <iostream>
#include <vector>
#include <cmath>
#include "model.h"
#include "optimizer/optimizer.h"
#include "optimizer/preprocessor/bisection.h"
#include "optimizer/step/Newton.h"
#include "optimizer/step/Geometric.h"
#include "optimizer/preprocessor/quadraticInterpolation.h"
#include "measurements/measurements.h"
#include "splittingAlgorithm/splittingAlgorithm.h"
#include "geometry/surface/compositeBicubicsrf.h"


double compareSurfaces(TopParametric top, bicubicsrf bez){
    double totalDist = 0;
    for (double i = 0; i <= 1; i += 0.01) {
        for (double j = 0; j <= 1; j += 0.01) {
            vec3d p1, p2;
            p1 = top.at(j,i);
            p2 = bez.at(j,i);
            cout << "param: " << p1 << endl;
            cout << "bez: " << p2 << endl;
            totalDist += p1.dist(p2);
        }
    }
    return totalDist/10000;
}

int main() {

    Model model = Model(
            12.0,
            59.5/180 * M_PI,
            11.4,
            14.4,
            10,
            1.2,
            50.0/180*M_PI,
            7.2,
            16.8,
            3.5,
            1.35,
            -0.2,
            -0.2,
            0.01,
            1.0,
            6.5);

    Model p1 = Model::getPart(model,0),
          p2 = Model::getPart(model,1),
          p3 = Model::getPart(model,2),
          p4 = Model::getPart(model,3),
          p5 = Model::getPart(model,4),
          p6 = Model::getPart(model,5);

    preprocessor *bisect;
    preprocessor *quad;
    step *newton;
    step *geometric;
    optimizer *opt;
    quad = dynamic_cast<preprocessor*>(new quadraticInterpolation(8));
    bisect = dynamic_cast<preprocessor*>(new bisection(6));
    geometric = dynamic_cast<step*>(new Geometric());
    newton = dynamic_cast<step*>(new Newton(1));

    TopParametric top = model.getTopParametric();
    BottomParametric bot = model.getBottomParametric();
    bicubicsrf bez = model.getTopBezier();
    surface *sur = dynamic_cast<surface*>(&top);
    opt = new optimizer(
            *bisect,
            *quad,
            *newton,
            *newton,
            *sur,
            {0,0,0},
            0.00000001,
            20,
            2
    );

    //TestOptimizerPerformance(*opt,3,2,1,5,1000,0.00000001,0);

    vec3d P = {1,1.5,1};

    vec3d ctrls[16] = {
            {0,0,0},{1,1.5,0},{2,1.5,0},{3,0,0},
            {0,0,1},{1,2.5,1},{2,2.5,1},{3,0,1},
            {0,0,2},{1,2.5,2},{2,2.5,2},{3,0,2},
            {0,0,3},{1,1.5,3},{2,1.5,3},{3,0,3}
    };

    bicubicsrf s(ctrls);

//    s.subdivideInDir(false,0.5,s1,s2);
//    plotSurface(s,0.1);
//    plotSurface(s1,0.1);
//    plotSurface(s2,0.1);
//
//    cout << s.closestPointInPatch(P) << endl;
//    cout << s1.closestPointInPatch(P) << endl;
//    cout << s2.closestPointInPatch(P) << endl;

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration;

    t_start = chrono::high_resolution_clock::now();
    cout << bez.hasValidControlNet() << endl;



    OptState2D loc = splittingAlgorithm::optimize(bez,P,0.00000001);
    t_stop = chrono::high_resolution_clock::now();

    cout << "loc: u: " << loc.u << " v: " << loc.v << " dist: " << loc.dist << endl;

    t_duration = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);

    cout << t_duration.count() << " microseconds" << endl;

    TestOptimizerPerformance(*opt,3,0,1,5,1000,0.00000001,0);
    //TestSplitterPerformance(bez, 3, 0, 3.0, 1, 1000, 0.00000001, 0);

//    cout << "param vs bez: " << compareSurfaces(top,bez) << endl;
//
//    compositeBicubicsrf b(bez);
//
//    cout << b.hasValidControlNet() << endl;

    delete quad;
    delete geometric;
    delete opt;
    delete bisect;
    delete newton;
}



