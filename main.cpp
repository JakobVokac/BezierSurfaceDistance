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


void plotSurface(surface &s, double eps = 0.1){
    if(eps >= 1){
        cout << "surface plotting: eps needs to be a fraction of 1" << endl;
        return;
    }
    vector<vector<double>> x,y,z;
    for (double i = 0; i <= 1; i+=eps) {
        vector<double> xr,yr,zr;
        for (double j = 0; j <= 1; j+=eps) {
            vec3d p = s.at(j,i);
            xr.push_back(p.getx());
            yr.push_back(p.gety());
            zr.push_back(p.getz());
        }
        x.push_back(xr);
        y.push_back(yr);
        z.push_back(zr);
    }
//    plt::plot_surface_with_line_and_axes_set(x,y,z,{-1,5,-1,5,-1,5},{},{},{},{},{},{},0);
    plt::plot_surface(x,y,z);
    plt::show();
}

void curveSplitting(double tl, double tr, double &tf, cubiccrv c, cubiccrv &cf, vec3d P, double eps, int i = 0, int limit = 100){
    bool valid = c.hasConvexPolygon();
    if(valid){
        if(!c.closestPointInCurve(P)){
            return;
        }else{
            vec3d c0c3 = c.getCtrlP(3) - c.getCtrlP(0);
            if(c0c3.mag() < eps || i > limit){
                cf = c;
                tf = (tl+tr)/2;
                return;
            }
        }
    }
    cubiccrv c1,c2;
    c.subdivideAt(0.5,c1,c2);
    curveSplitting(tl, (tl+tr)/2, tf, c1, cf, P, eps, i + 1, limit);
    curveSplitting((tl+tr)/2, tr, tf, c2, cf, P, eps, i + 1, limit);
}

void surfaceSplitting(int dir, double ul, double ur, double vl, double vr, double &uf, double &vf, bicubicsrf s, bicubicsrf &sf, vec3d P, double eps, int i = 0, int limit = 100){
    bool valid = s.hasValidControlNet();
    if(valid){
        if(!s.closestPointInPatch(P)){
            return;
        }else{
            vec3d diag = s.ctrlP(0) - s.ctrlP(15);
            if(diag.mag() < eps || i > limit){
                sf = s;
                uf = (ul+ur)/2;
                vf = (vl+vr)/2;
                return;
            }
        }
    }
    bicubicsrf s1,s2;

    s.subdivide(dir,0.5,s1,s2);

//    cout << "uvlr: " << ul << " " << ur << " " << vl << " " << vr << " i: " << i << endl;
//    plotSurface(s1,0.1);
//    plotSurface(s2,0.1);
    if(!dir){
        surfaceSplitting((dir + 1)%2,ul,(ul+ur)/2,vl,vr,uf,vf,s1,sf,P,eps,i+1,limit);
        surfaceSplitting((dir + 1)%2,(ul+ur)/2,ur,vl,vr,uf,vf,s2,sf,P,eps,i+1,limit);
    }else{
        surfaceSplitting((dir + 1)%2,ul,ur,vl,(vl+vr)/2,uf,vf,s1,sf,P,eps,i+1,limit);
        surfaceSplitting((dir + 1)%2,ul,ur,(vl+vr)/2,vr,uf,vf,s2,sf,P,eps,i+1,limit);
    }
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
            *bisect,
            *newton,
            *newton,
            *sur,
            {0,0,0},
            0.00000001,
            20,
            2
    );

    //TestOptimizerPerformance(*opt,3,2,1,5,1000,0.00000001,0);
    cubiccrv C =  {{0,0,0},{1,1.5,0},{2,1.5,0},{3,0,0}};
    cubiccrv Cf;
    vec3d P = {1,1.5,1};
    double tf = -1;
    curveSplitting(0,1,tf,C,Cf,P,0.00000001,0,100);
    cout << Cf << endl;
    cout << tf << endl;
    vec3d v = C.f(tf);
    cout << v.dist(P) << endl;
    cout << Cf.f(0.5).dist(P) << endl;
    cout << C.df(tf).dot(C.getCtrlP(1)-C.f(tf)) << endl;

    bicubicsrf s1, s2;
    vec3d ctrls[16] = {
            {0,0,0},{1,1.5,0},{2,1.5,0},{3,0,0},
            {0,0,1},{1,2.5,1},{2,2.5,1},{3,0,1},
            {0,0,2},{1,2.5,2},{2,2.5,2},{3,0,2},
            {0,0,3},{1,1.5,3},{2,1.5,3},{3,0,3}
    };

    bicubicsrf s(ctrls);

//    s.subdivide(false,0.5,s1,s2);
    plotSurface(s,0.1);
//    plotSurface(s1,0.1);
//    plotSurface(s2,0.1);
//
//    cout << s.closestPointInPatch(P) << endl;
//    cout << s1.closestPointInPatch(P) << endl;
//    cout << s2.closestPointInPatch(P) << endl;

    double uf,vf;
    bicubicsrf sf;
    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration;

    t_start = chrono::high_resolution_clock::now();
    surfaceSplitting(1,0,1,0,1,uf,vf,s,sf,P,0.000000001,0,100);
    t_stop = chrono::high_resolution_clock::now();

    t_duration = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);

    cout << t_duration.count() << " microseconds" << endl;
    cout << "uv: " << uf << " " << vf << endl;
    plotSurface(sf,0.1);
}



