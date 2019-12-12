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

    BottomParametric b = p1.getBottomParametric();

    vector<vector<double>> x,y,z;
    for (double i = 0; i <= 1; i += 0.05) {
        vector<double> xs, ys, zs;
        for (double j = 0; j <= 1; j += 0.05) {
            vec3d p = p1.fillBottom(j,i);
            xs.push_back(p.getx());
            ys.push_back(p.gety());
            zs.push_back(p.getz());
        }
        x.push_back(xs);
        y.push_back(ys);
        z.push_back(zs);
    }
    plt::plot_surface(x,y,z);

    x.clear(),y.clear(),z.clear();
    for (double i = 0; i <= 1; i += 0.05) {
        vector<double> xs, ys, zs;
        for (double j = 0; j <= 1; j += 0.05) {
            vec3d p = b.at(j,i);
            xs.push_back(p.getx());
            ys.push_back(p.gety());
            zs.push_back(p.getz());
        }
        x.push_back(xs);
        y.push_back(ys);
        z.push_back(zs);
    }
    plt::plot_surface(x,y,z);
    plt::show();

    cout << p1.symCurveBottom;

    cout << p1.sinCurveBottom;

    cout << p1.Q_r;

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
            *bisect,
            *newton,
            *newton,
            *sur,
            {0,0,0},
            0.00000001,
            20,
            2
    );

    TestOptimizerPerformance(*opt,3,2,1,5,1000,0.00000001,0);


}



