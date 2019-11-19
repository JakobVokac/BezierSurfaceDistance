#include <iostream>
#include <vector>
#include <chrono>
#include <float.h>
#include "Model.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;



double bisection(Model model, vector<double> A, double min_u, double min_v, double step, double eps=0.01){
    double min_distance = 100;
    double distance_old = 105;

    int loop = 1;

    while(distance_old - min_distance > eps && loop < 6){
        double u_start = max(0.0, min_u - step);
        double u_end = max(0.0, min_u + step);
        double v_start = max(0.0, min_v - step);
        double v_end = max(0.0, min_v + step);
        step /= 2.0;

        distance_old = min_distance;
        double u = u_start;
        while(u <= u_end){
            double v = v_start;
            while(v <= v_end){
                double distance = model.distanceToTopPoint(u, v, A);
                if(distance < min_distance){
                    min_distance = distance;
                    min_u = u;
                    min_v = v;
                }
                v += step;
            }
            u += step;
        }
        loop++;
    }

    return min_distance;
}
vector<double> computeInitialUV(Model model, vector<double> A, double step){
    double min_distance = 100;
    double u = 0;
    double u_min, v_min;
    while(u < 1.0001){
        double v = 0;
        while(v < 1.0001){
            double distance = model.distanceToTopPoint(u, v, A);
            if(distance < min_distance){
                min_distance = distance;
                u_min = u;
                v_min = v;
            }
            v += step;
        }
        u += step;
    }
    return vector<double>({u_min,v_min});
}
double distanceToCusp(Model model, vector<double> A, double eps = 0.01){
    double step = 0.025;

    vector<double> initUV = computeInitialUV(model,A,step);
    return model.distanceToTopPoint(initUV[0], initUV[1], A);
}
void drawPart(Model &model, vector<vector<double>> &x, vector<vector<double>> &y, vector<vector<double>> &z) {
    for (double i = 1; i >= 0; i -= 0.1) {
        vector<double> x_row, y_row, z_row;
        for (double j = 0; j <= 1; j += 0.1) {
            vector<double> p = model.fillTop(j, i);
//            cout << p[0] << " " << p[1] << " " << p[2] << endl;
            x_row.push_back(p[0]);
            y_row.push_back(p[1]);
            z_row.push_back(p[2]);
        }
        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }

    for (double i = 0.1; i <= 1;  i += 0.1) {
        vector<double> x_row, y_row, z_row;
        for (double j = 0; j <= 1; j += 0.1) {
            vector<double> p = model.fillBottom(j, i);
//            cout << p[0] << " " << p[1] << " " << p[2] << endl;
            x_row.push_back(p[0]);
            y_row.push_back(p[1]);
            z_row.push_back(p[2]);
        }
        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }
}

double f(double x, double y){
    return sin(x)*x*x + sin(y)*y*y;
}
double fdx(double x, double y){
    return x*x*cos(x) + 2*x*sin(x);
}
double fdy(double x, double y){
    return y*y*cos(y) + 2*y*sin(y);
}
double fdxdx(double x, double y){
    return -x*x*sin(x) + 4*x*cos(x) + 2*sin(x);
}
double fdydy(double x, double y){
    return -y*y*sin(y) + 4*y*cos(y) + 2*sin(y);
}
double fdxdy(double x, double y){
    return 0;
}
double NewtonMethodTest(double x, double y, double eps=0.01, double sigma = 0.1){

    double V = f(x,y);
    double V_old;
    int iter = 0;

    vector<double> xs, ys;
    xs.push_back(x);
    ys.push_back(y);
    vector<double> grad(2);
    vector<vector<double>> H = {{0, 0}, {0, 0}}, Ht;
    do{
        iter++;
        grad[0] = fdx(x,y);
        grad[1] = fdy(x,y);

        H[0][0] = fdxdx(x,y);
        H[0][1] = fdxdy(x,y);
        H[1][0] = fdxdy(x,y);
        H[1][1] = fdydy(x,y);

//        cout << "Derivatives:" << endl;
//        cout << grad[0] << " " << grad[1] << endl;
//        cout << H[0][0] << " " << H[0][1] << " " << H[1][0] << " " << H[1][1] << endl;

        double invConst = 1/(H[0][0] * H[1][1] - H[1][0] * H[0][1]);

        Ht = H;
        H[0][0] = invConst * Ht[1][1];
        H[0][1] = invConst * (-Ht[0][1]);
        H[1][0] = invConst * (-Ht[1][0]);
        H[1][1] = invConst * Ht[0][0];

//        cout << "Delta x: " << -sigma*(H[0][0] * grad[0] + H[0][1] * grad[1]) << endl;
//        cout << "Delta y: " << -sigma*(H[1][0] * grad[0] + H[1][1] * grad[1]) << endl;
        x -= sigma*(H[0][0] * grad[0] + H[0][1] * grad[1]);
        y -= sigma*(H[1][0] * grad[0] + H[1][1] * grad[1]);

        V_old = V;
        V = f(x,y);
//        cout << "V: " << V << endl;
        xs.push_back(x);
        ys.push_back(y);
    }while(abs(V_old - V) > eps);
    plt::plot(xs,ys);
    return V;
}
double func(double x){
    return pow((x-2.5),2.0) + 2*x + pow(x, 3.0);
}
double funcdx(double x) {
    return 2 * (x-2.5) + 2 + 3*pow(x, 2.0);
}
double funcdxdx(double x){
    return 2 + 6*x;
};
double NewtonMethodSimpleFunc(double x, double eps){
    double V = func(x);
    double V_old;
    double t = 0;

    vector<double> xs,ys, ys2, ys3;
    cout << "Fletcher Reeves Newton Test" << endl;
    for(double i = 0; i < 1; i+= 0.01 ){
        xs.push_back(i);
        ys.push_back(func(i));
        ys2.push_back(funcdx(i));
        ys3.push_back(funcdxdx(i));
    }
    plt::plot(xs,ys,"g");
    plt::plot(xs,ys2,"b");
    plt::plot(xs,ys3,"r");

    plt::show();
    do{
        cout << x << endl;
        x -= (funcdx(x))
             /(funcdxdx(x));
//        cout << "t: " << t << endl;
        V_old = V;
        V = func(x);
    }while(abs(V-V_old) > eps);

    return x;
}


double constraintFunc(double u, double v){
    return pow((u < 0 ? -u : 0),3.0)
           + pow((u > 1 ? u-1 : 0),3.0)
           + pow((v < 0 ? -v : 0),3.0)
           + pow((v > 1 ? v-1 : 0),3.0);
}
double constraintFuncdU(double u, double v){
    return (u > 1 ? 3*pow((u-1),2.0) : 0) + (u < 0 ? -3*pow(u,2.0) : 0);
}
double constraintFuncdV(double u, double v){
    return (v > 1 ? 3 * pow((v - 1), 2.0) : 0) + (v < 0 ? -3 * pow(v, 2.0) : 0);
}
double constraintFuncdUdU(double u, double v){
    return (u > 1 ? 6*(u-1) : 0) + (u < 0 ? -6*u : 0);
}
double constraintFuncdVdV(double u, double v){
    return (u > 1 ? 6*(u-1) : 0) + (u < 0 ? -6*u : 0);
}
double constraintFuncdUdV(double u, double v){
    return 0;
}

vector<double> gradientDescent(Model model, vector<double> A, double u, double v, double eps=0.01, double gamma = 0.1, double lambda = 100){
    double dist = model.distanceToTopPoint(u,v,A) + lambda*constraintFunc(u,v);
    double u0 = u, v0 = v;
    double dist_old;

    int iter = 0;

    vector<double> us, vs;
    us.push_back(u);
    vs.push_back(v);

    do{
//        cout << "u: " << u << " v: " << v << " A: x:" << A[0] << " y:" << A[1] << " z:" << A[2] << endl;
//        cout << "Iter: " << iter << " Delta u: " << model.distanceToTopPointDerU(u,v,A) << " Delta v: " << model.distanceToTopPointDerV(u,v,A) << endl;
        u -= gamma * (model.completeDistanceTopDerU(u,v,A) + lambda*constraintFuncdU(u,v));
        v -= gamma * (model.completeDistanceTopDerV(u,v,A) + lambda*constraintFuncdV(u,v));


        dist_old = dist;
        dist = model.distanceToTopPoint(u,v,A) + lambda*constraintFunc(u,v);
        us.push_back(u);
        vs.push_back(v);
        iter++;
    }while(abs(dist_old - dist) > eps);
    plt::plot(us,vs);
//    cout << "Gradient descent iterations: " << iter << " u: " << u << " v: " << v << endl;
    return {u,v,dist,u0,v0,static_cast<double>(iter)};
}

vector<double> Newton1DV(Model model, vector<double> A, double u, double v, double u0, double v0, int iter, double eps=0.01){
    double dist = model.distanceToTopPoint(u,v,A);
    double dist_old;
    if(v < 0){
        v = 0;
    }else if(v > 1){
        v = 1;
    }
    do{
        iter++;
        v -= model.squaredTopDistDerV(u,v,A)/model.squaredTopDistDerVV(u, v, A);
        dist_old = dist;
        dist = model.distanceToTopPoint(u,v,A);
        if(v < 0){
            return {u,0,model.distanceToTopPoint(u,0,A),u0,v0,static_cast<double>(iter)};
        }else if(v > 1){
            return {u,1,model.distanceToTopPoint(u,1,A),u0,v0,static_cast<double>(iter)};
        }
    }while(abs( model.squaredTopDistDerV(u,v,A)) > eps && iter < 100);

    return {u,v,dist,u0,v0,static_cast<double>(iter)};
}

vector<double> Newton1DU(Model model, vector<double> A, double u, double v, double u0, double v0, int iter, double eps=0.01){
    double dist = model.distanceToTopPoint(u,v,A);
    double dist_old;
    if(u < 0){
        u = 0;
    }else if(u > 1){
        u = 1;
    }
    do{
        iter++;
        u -= model.squaredTopDistDerU(u,v,A)/model.squaredTopDistDerUU(u, v, A);
        dist_old = dist;
        dist = model.distanceToTopPoint(u,v,A);
        if(u < 0){
            return {0,v,model.distanceToTopPoint(0,v,A),u0,v0,static_cast<double>(iter)};
        }else if(u > 1){
            return {1,v,model.distanceToTopPoint(1,v,A),u0,v0,static_cast<double>(iter)};
        }
    }while(abs(model.squaredTopDistDerU(u,v,A)) > eps && iter < 100);

    return {u,v,dist,u0,v0,static_cast<double>(iter)};
}

vector<double> NewtonMethodProjected(Model model, vector<double> A, double u, double v, double eps=0.01, double sigma = 0.1, double HModifier = 1600){

    double dist = model.distanceToTopPoint(u,v,A);
    double u0 = u, v0 = v;
    double dist_old;
    int iter = 0;

    vector<double> us, vs;
    us.push_back(u);
    vs.push_back(v);
    vector<double> grad(2);
    vector<vector<double>> H = {{0, 0}, {0, 0}}, Ht;
    do {
//        cout << "Distance at iteration " << iter << ": " << dist << " u: " << u << " v: " << v << endl;
        iter++;
//        cout << "u: " << u << " v: " << v << endl;
        grad[0] = model.squaredTopDistDerU(u, v, A);
        grad[1] = model.squaredTopDistDerV(u, v, A);
//        cout << "du: " << grad[0] << " dv: " << grad[1] << endl;
        H[0][0] = model.squaredTopDistDerUU(u, v, A) + HModifier;
        H[0][1] = model.squaredTopDistDerUV(u, v, A);
        H[1][0] = model.squaredTopDistDerUV(u, v, A);
        H[1][1] = model.squaredTopDistDerVV(u, v, A) + HModifier;

//        cout << "dudu: " << H[0][0] << " dudv: " << H[0][1] << " dvdv: " << H[1][1] << endl;
        double invConst = 1 / (H[0][0] * H[1][1] - H[1][0] * H[0][1]);
//        cout << "determinant: " << invConst << endl;
        Ht = H;
        H[0][0] = invConst * Ht[1][1];
        H[0][1] = invConst * (-Ht[0][1]);
        H[1][0] = invConst * (-Ht[1][0]);
        H[1][1] = invConst * Ht[0][0];
//        cout << "dudu: " << H[0][0] << " dudv: " << H[0][1] << " dvdv: " << H[1][1] << endl;

//        cout << "Delta u: " << -sigma*(H[0][0] * grad[0] + H[0][1] * grad[1]) << endl;
//        cout << "Delta v: " << -sigma*(H[1][0] * grad[0] + H[1][1] * grad[1]) << endl;

        u -= sigma * (H[0][0] * grad[0] + H[0][1] * grad[1]);
        v -= sigma * (H[1][0] * grad[0] + H[1][1] * grad[1]);
        us.push_back(u);
        vs.push_back(v);
        if (u < 0) {
//            plt::plot(us,vs);
//            plt::show();
            return Newton1DV(model, A, 0, v, u0, v0, iter, eps);
        }
        if (u > 1) {
//            plt::plot(us,vs);
//            plt::show();
            return Newton1DV(model, A, 1, v, u0, v0, iter, eps);
        }
        if (v < 0) {
//            plt::plot(us,vs);
//            plt::show();
            return Newton1DU(model, A, u, 0, u0, v0, iter, eps);
        }
        if (v > 1) {
//            plt::plot(us,vs);
//            plt::show();
            return Newton1DU(model, A, u, 1, u0, v0, iter, eps);
        }
        dist_old = dist;
        dist = model.distanceToTopPoint(u, v, A);
        HModifier /= 2;
    }while(abs(grad[0]) + abs(grad[1]) > eps && iter < 10);

//    cout << "Newton's method iterations: " << iter << " u: " << u << " v: " << v << endl;
    return {u,v,dist,u0,v0,static_cast<double>(iter)};
}

vector<double> NewtonMethod(Model model, vector<double> A, double u, double v, double eps=0.01, double sigma = 0.1, double lambda = 100, double HModifier = 1600){

    double dist = model.distanceToTopPoint(u,v,A) + lambda*constraintFunc(u,v);
    double u0 = u, v0 = v;
    double dist_old;
    int iter = 0;

    vector<double> us, vs;
    us.push_back(u);
    vs.push_back(v);
    vector<double> grad(2);
    vector<vector<double>> H = {{0, 0}, {0, 0}}, Ht;
    do{
//        cout << "Distance at iteration " << iter << ": " << dist << " u: " << u << " v: " << v << endl;
        iter++;
//        cout << "u: " << u << " v: " << v << endl;
        grad[0] = model.completeDistanceTopDerU(u,v,A) + lambda*constraintFuncdU(u,v);
        grad[1] = model.completeDistanceTopDerV(u,v,A) + lambda*constraintFuncdV(u,v);
//        cout << "du: " << grad[0] << " dv: " << grad[1] << endl;
        H[0][0] = model.completeDistanceTopDerUU2(u, v, A) + lambda*constraintFuncdUdU(u,v);
        H[0][1] = model.completeDistanceTopDerUV2(u, v, A) + lambda*constraintFuncdUdU(u,v);
        H[1][0] = model.completeDistanceTopDerUV2(u, v, A) + lambda*constraintFuncdUdU(u,v);
        H[1][1] = model.completeDistanceTopDerVV2(u, v, A) + lambda*constraintFuncdUdU(u,v);

//        cout << "dudu: " << H[0][0] << " dudv: " << H[0][1] << " dvdv: " << H[1][1] << endl;
        double invConst = 1/(H[0][0] * H[1][1] - H[1][0] * H[0][1]);
//        cout << "determinant: " << invConst << endl;
        Ht = H;
        H[0][0] = invConst * Ht[1][1];
        H[0][1] = invConst * (-Ht[0][1]);
        H[1][0] = invConst * (-Ht[1][0]);
        H[1][1] = invConst * Ht[0][0];
//        cout << "dudu: " << H[0][0] << " dudv: " << H[0][1] << " dvdv: " << H[1][1] << endl;

//        cout << "Delta u: " << -sigma*(H[0][0] * grad[0] + H[0][1] * grad[1]) << endl;
//        cout << "Delta v: " << -sigma*(H[1][0] * grad[0] + H[1][1] * grad[1]) << endl;

        u -= sigma*(H[0][0] * grad[0] + H[0][1] * grad[1]);
        v -= sigma*(H[1][0] * grad[0] + H[1][1] * grad[1]);
        if( u < 0 ){
            u = 0;
        }
        if( v < 0 ){
            v = 0;
        }
        if( u > 1 ){
            u = 1;
        }
        if( v > 1 ){
            v = 1;
        }
        dist_old = dist;
        dist = model.distanceToTopPoint(u,v,A) + lambda*constraintFunc(u,v);
        us.push_back(u);
        vs.push_back(v);
        HModifier /= 2;
    }while(abs(dist_old - dist) > eps);
//    plt::plot(us,vs);
//    plt::show();
//    cout << "Newton's method iterations: " << iter << " u: " << u << " v: " << v << endl;
    return {u,v,dist,u0,v0,static_cast<double>(iter)};
}

double NewtonForLinearCG(Model model, vector<double> A, double u, double v, double dU, double dV, double eps){
    double V = model.distanceToTopPoint(u,v,A);
    double V_old;
    double t = 0;


    do{
        t -= (dU * model.squaredTopDistDerU(u + dU*t,v + dV*t,A) + dV * model.squaredTopDistDerV(u + dU*t,v + dV*t,A))
                /(dU*dU*model.squaredTopDistDerUU(u + dU*t,v + dV*t,A)
                + 2*dU*dV*model.squaredTopDistDerUV(u + dU*t,v + dV*t,A)
                + dV*dV*model.squaredTopDistDerVV(u + dU*t,v + dV*t,A));

//        cout << "t: " << t << endl;
        V_old = V;
        V = model.distanceToTopPoint(u + dU*t, v + dV*t, A);
    }while(abs(V-V_old) > eps && 0 <= u + dU*t && u + dU*t <= 1 && v + dV*t <= 1 && 0 <= v + dV*t);

    vector<double> xs,ys, ys2, ys3;

    /*cout << "Fletcher Reeves Newton" << endl;
    cout << u << " " << dU << " " << -u/dU << " " << (1-u)/dU << (1-2*u)/dU * 0.01 << endl;
    for(double i = 0; i < 1; i+= 0.01 ){
        xs.push_back(i);
        ys.push_back( model.squaredTopDist(u + dU*((1-2*u)/dU * i), v + dV*((1-2*v)/dV * i), A));
        ys2.push_back((dU * model.squaredTopDistDerU(u + dU*((1-2*u)/dU * i),v + dV*((1-2*v)/dV * i),A) + dV * model.squaredTopDistDerV(u + dU*((1-2*u)/dU * i),v + dV*((1-2*v)/dV * i),A))/100);
        ys3.push_back((dU*dU*model.squaredTopDistDerUU(u + dU*((1-2*u)/dU * i),v + dV*((1-2*v)/dV * i),A)
                       + 2*dU*dV*model.squaredTopDistDerUV(u + dU*((1-2*u)/dU * i),v + dV*((1-2*v)/dV * i),A)
                       + dV*dV*model.squaredTopDistDerVV(u + dU*((1-2*u)/dU * i),v + dV*((1-2*v)/dV * i),A))/10000);
        cout << model.squaredTopDist(u + dU*((1-2*u)/dU * i), v + dV*((1-2*u)/dU * i), A) << endl;
    }
    plt::plot(xs,ys,"g");
    plt::plot(xs,ys2,"b");
    plt::plot(xs,ys3,"r");

    plt::show();*/

    return t;

}

double SecantForLinearCG(Model model, vector<double> A, double u0, double v0, double t2, double dU, double dV, double eps){

    double t1 = 0;
    double f1 = dU * model.completeDistanceTopDerU(u0 + dU*t2,v0 + dV*t2,A) + dV * model.completeDistanceTopDerV(u0 + dU*t2,v0 + dV*t2,A);
    double f2;

    do{
        f2 = f1;
        f1 = dU * model.completeDistanceTopDerU(u0 + dU*t1,v0 + dV*t1,A) + dV * model.completeDistanceTopDerV(u0 + dU*t1,v0 + dV*t1,A);
        double temp = t1;
        t1 = t1 - f1 * (t1 - t2)/(f1-f2);
        t2 = temp;

    }while(abs(t1-t2) > eps && 0 <= u0 + dU*t1 && u0 + dU*t1 <= 1 && v0 + dV*t1 <= 1 && 0 <= v0 + dV*t1);


    return t1;
}

void DrawDistWithDirection(Model &model, const vector<double> &P) {
    vector<vector<double>> u, v, dist;

    double min_dist = 100000;
    double max_dist = -100000;
    for (double i = 0; i <= 1.05; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j <= 1.05; j += 0.05) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(model.squaredTopDist(j, i, P));
            if (min_dist > model.squaredTopDist(j, i, P)) {
                min_dist = model.squaredTopDist(j, i, P);
            }
            if (max_dist < model.squaredTopDist(j, i, P)) {
                max_dist = model.squaredTopDist(j, i, P);
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    vector<double> x2, y2, z2, u2, v2, w2;
    for (double i = 0; i <= 1.05; i += 0.05) {
        for (double j = 0; j <= 1.05; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.squaredTopDistDerU(j, i, P);
            tV = model.squaredTopDistDerV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);

    plt::show();

    u.clear(), v.clear(), dist.clear();

    min_dist = 100000;
    max_dist = -100000;
    for (double i = 0; i <= 1.05; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j <= 1.05; j += 0.05) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(model.squaredTopDistDerU(j, i, P));
            if (min_dist > model.squaredTopDistDerU(j, i, P)) {
                min_dist = model.squaredTopDistDerU(j, i, P);
            }
            if (max_dist < model.squaredTopDistDerU(j, i, P)) {
                max_dist = model.squaredTopDistDerU(j, i, P);
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    x2.clear(), y2.clear(), z2.clear(), u2.clear(), v2.clear(), w2.clear();
    for (double i = 0; i <= 1.05; i += 0.05) {
        for (double j = 0; j <= 1.05; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.squaredTopDistDerUU(j, i, P);
            tV = model.squaredTopDistDerUV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);

    plt::show();

    u.clear(), v.clear(), dist.clear();

    min_dist = 100000;
    max_dist = -100000;
    for (double i = 0; i <= 1.05; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j <= 1.05; j += 0.05) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(model.squaredTopDistDerV(j, i, P));
            if (min_dist > model.squaredTopDistDerV(j, i, P)) {
                min_dist = model.squaredTopDistDerV(j, i, P);
            }
            if (max_dist < model.squaredTopDistDerV(j, i, P)) {
                max_dist = model.squaredTopDistDerV(j, i, P);
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    x2.clear(), y2.clear(), z2.clear(), u2.clear(), v2.clear(), w2.clear();
    for (double i = 0; i <= 1.05; i += 0.05) {
        for (double j = 0; j <= 1.05; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.squaredTopDistDerUV(j, i, P);
            tV = model.squaredTopDistDerVV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);

    plt::show();

}



vector<double> BroydenMethod(Model model, vector<double> A, double u1, double v1, double u2, double v2, double eps=0.0001, int plot = 0){
    double dist = model.squaredTopDist(u1,v1,A);
    double old_dist = model.squaredTopDist(u2,v2,A);
    int iter = 0;
    while(abs(dist-old_dist) > eps && iter < 0){



    }
}





vector<double> FletcherReevesCG(Model model, vector<double> A, double u, double v, double eps=0.0001, int plot = 0){
    double dU, dV, alpha, beta, dU2, dV2, u0 = u, v0 = v;
    int iter = 0;
    //cout << "u: " << u << " v: " << v << endl;

    // p 0 <- - Grad f 0 = Grad f(u,v)
    dU2 = -model.squaredTopDistDerU(u,v,A);
    dV2 = -model.squaredTopDistDerV(u,v,A);
    vector<double> xs, ys, is, js,zs,ks;
    xs.push_back(u);
    ys.push_back(v);
    zs.push_back(0);
    is.push_back(dU2);
    js.push_back(dV2);
    ks.push_back(0);
    do{
        dU = dU2;
        dV = dV2;

        // Compute alpha k and set xk1 = xk + alpha k * p k
        alpha = NewtonForLinearCG(model,A,u,v,dU,dV,eps);
//        cout << "alpha: " << alpha << endl;
//        cout << "dU: " << dU << " dV: " << dV << endl;
        if(plot){
            cout << "u: " << u << " v: " << v << " dU: " << dU << " dV: " << dV << " alpha: " << alpha << endl;
        }
        u += alpha * dU;
        v += alpha * dV;
        if (u < 0) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,0,v,0,1,eps);
            cout << "u: " << u << " v: " << v << endl;
            DrawDistWithDirection(model,A);
            return {0,v + alpha, model.distanceToTopPoint(0,v + alpha,A), u0, v0, static_cast<double>(iter)};
        }
        if (u > 1) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,1,v,0,1,eps);
            cout << "u: " << u << " v: " << v << endl;
            DrawDistWithDirection(model,A);
            return {1,v + alpha, model.distanceToTopPoint(1,v + alpha,A), u0, v0, static_cast<double>(iter)};
        }
        if (v < 0) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,u,0,1,0,eps);
            cout << "u: " << u << " v: " << v << endl;
            DrawDistWithDirection(model,A);
            return {u + alpha,0, model.distanceToTopPoint(u + alpha,0,A), u0, v0, static_cast<double>(iter)};
        }
        if (v > 1) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,u,1,1,0,eps);
            cout << "u: " << u << " v: " << v << endl;
            DrawDistWithDirection(model,A);
            return {u + alpha,1, model.distanceToTopPoint(u + alpha,1,A), u0, v0, static_cast<double>(iter)};
        }

        //beta k1 <- (Grad fk1 T * Grad fk1)/(Grad fk T * Grad fk)
        dU2 = -model.squaredTopDistDerU(u,v,A);
        dV2 = -model.squaredTopDistDerV(u,v,A);
        beta = (dU2*dU2 + dV2*dV2)/(dU*dU + dV*dV);

        //p k1 <- - Grad fk1 + beta k1 * p k
        dU2 += beta * dU;
        dV2 += beta * dV;

        //k <- k + 1
        iter++;
        xs.push_back(u);
        ys.push_back(v);
        zs.push_back(50);
        is.push_back(dU2);
        js.push_back(dV2);
        ks.push_back(0);
//        cout << "u: " << u << " v: " << v << endl;
    }while(abs(dU + dV) > eps && iter < 10); // while Grad
    cout << "u: " << u << " v: " << v << endl;
    DrawDistWithDirection(model,A);

    if(plot){
        plt::plot(xs,ys, "g");
        plt::show();
    }
    return {u,v,model.distanceToTopPoint(u,v,A), u0, v0, static_cast<double>(iter)};
}

vector<double> FletcherReevesCGSecant(Model model, vector<double> A, double u, double v, double alpha, double eps=0.0001, int plot = 0, double initDenom = 3){
    double dU, dV, beta, dU2, dV2, u0 = u, v0 = v;
    int iter = 0;
    //cout << "u: " << u << " v: " << v << endl;

    // p 0 <- - Grad f 0 = Grad f(u,v)
    dU2 = -model.completeDistanceTopDerU(u,v,A);
    dV2 = -model.completeDistanceTopDerV(u,v,A);
    vector<double> xs, ys;
    xs.push_back(u);
    ys.push_back(v);
    do{
        dU = dU2;
        dV = dV2;
        // Compute alpha k and set xk1 = xk + alpha k * p k
        alpha = SecantForLinearCG(model,A,u,v,alpha/initDenom, dU,dV,eps);
//        cout << "alpha: " << alpha << endl;
//        cout << "dU: " << dU << " dV: " << dV << endl;

        if(plot){
            cout << "u: " << u << " v: " << v << " dU: " << dU << " dV: " << dV << " alpha: " << alpha << endl;
        }
        u += alpha * dU;
        v += alpha * dV;
        if (u < 0) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = SecantForLinearCG(model,A,0,v,alpha/initDenom, 0,1,eps);
            return {0,v + alpha, model.distanceToTopPoint(0,v + alpha,A), u0, v0, static_cast<double>(iter)};
        }
        if (u > 1) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = SecantForLinearCG(model,A,1,v,alpha/initDenom, 0,1,eps);
            return {1,v + alpha, model.distanceToTopPoint(1,v + alpha,A), u0, v0, static_cast<double>(iter)};
        }
        if (v < 0) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = SecantForLinearCG(model,A,u,0,alpha/initDenom, 1,0,eps);
            return {u + alpha,0, model.distanceToTopPoint(u + alpha,0,A), u0, v0, static_cast<double>(iter)};
        }
        if (v > 1) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = SecantForLinearCG(model,A,u,1,alpha/initDenom, 1,0,eps);
            return {u + alpha,1, model.distanceToTopPoint(u + alpha,1,A), u0, v0, static_cast<double>(iter)};
        }

        //beta k1 <- (Grad fk1 T * Grad fk1)/(Grad fk T * Grad fk)
        dU2 = -model.completeDistanceTopDerU(u,v,A);
        dV2 = -model.completeDistanceTopDerV(u,v,A);
        beta = (dU2*dU2 + dV2*dV2)/(dU*dU + dV*dV);

        //p k1 <- - Grad fk1 + beta k1 * p k
        dU2 += beta * dU;
        dV2 += beta * dV;

        //k <- k + 1
        iter++;
        xs.push_back(u);
        ys.push_back(v);
//        cout << "u: " << u << " v: " << v << endl;
    }while(abs(dU + dV) > eps && iter < 10); // while Grad
//    cout << "u: " << u << " v: " << v << endl;
    if(plot){
        plt::plot(xs,ys, "g");
        plt::show();
    }
    return {u,v,model.distanceToTopPoint(u,v,A), u0, v0, static_cast<double>(iter)};
}

vector<double> getBSPatch(double um, double vm, double ul, double vl, double ur, double vr){
    double dU, dV;
    if(um == ul) {
        dU = ur - um;
    }else if(um == ur){
        dU = um - ul;
    }else if(um - ul < ur - um){
        dU = um - ul;
    }else{
        dU = ur - um;
    }
    if(vm == vl) {
        dV = vr - vm;
    }else if(vm == vr){
        dV = vm - vl;
    }else if(vm - vl < vr - vm){
        dV = vm - vl;
    }else{
        dV = vr - vm;
    }

    return {um, vm, dU, dV};
}

vector<double> findStartingPoint(Model model, vector<double> A, int iterations = 12){
    double ul = 0, ur = 1, vl = 0, vr = 1;
    double dist = model.distanceToTopPoint(0.5,0.5,A);
    double um, vm;
    for (int i = 0; i < iterations; ++i) {
        if( i % 2 == 0 ) {
            um = (ul+ur)/2;
            vm = (vl+vr)/2;
            double u1 = (ul+um)/2;
            double u2 = (ur+um)/2;
            double dist1 = model.distanceToTopPoint(u1,vm,A);
            double dist2 = model.distanceToTopPoint(u2,vm,A);
            if(dist1 < dist2){
                ur = um;
            }else{
                ul = um;
            }
        }else{
            um = (ul+ur)/2;
            vm = (vl+vr)/2;
            double v1 = (vl + vm)/2;
            double v2 = (vr + vm)/2;
            double dist1 = model.distanceToTopPoint(um,v1,A);
            double dist2 = model.distanceToTopPoint(um,v2,A);
            if(dist1 < dist2){
                vr = vm;
            }else{
                vl = vm;
            }
        }
    }
    return getBSPatch(um,vm,ul,vl,ur,vr);
}


vector<double> quadInterpolation1(Model m, vector<double> P, int iter){
    double u1 = 0.5, u2 = 0, v1 = 0.5, v2 = 0;
    double h1,h2,h1d;
    for (int i = 0; i < iter; ++i) {

        if(u1 != 0 && u1 != 1) {
        h1 = m.squaredTopDist(u1, v1, P);
        h2 = m.squaredTopDist(u2, v1, P);
        h1d = m.squaredTopDistDerU(u1, v1, P);
        cout << "u h1: " << h1 << " h2: " << h2 << " h1d: " << h1d << endl;

        double ut = u1;
        cout << "u1: " << u1 << " u2: " << u2 << " h1: " << h1 << " h2: " << h2 << " h1d: " << h1d << endl;

            u1 = u1 - 1.0 / 2.0 * ((u1 - u2) * h1d) / (h1d - (h1 - h2) / (u1 - u2));
            if (u1 < 0) {
                u1 = 0;
            } else if (u1 > 1) {
                u1 = 1;
            }
            u2 = ut;
        }

        cout << "new u: " << u1 << " new v: " << v1 << endl;

        if(v1 != 0 && v1 != 1) {

            h1 = m.squaredTopDist(u1, v1, P);
            h2 = m.squaredTopDist(u1, v2, P);
            h1d = m.squaredTopDistDerV(u1, v1, P);

            cout << "v h1: " << h1 << " h2: " << h2 << " h1d: " << h1d << endl;

            double vt = v1;

            v1 = v1 - 1.0 / 2.0 * ((v1 - v2) * h1d) / (h1d - (h1 - h2) / (v1 - v2));
            if (v1 < 0) {
                v1 = 0;
            } else if (v1 > 1) {
                v1 = 1;
            }
            v2 = vt;

        }

        cout << "new u: " << u1 << " new v: " << v1 << endl;
    }
    return {u1,v1};
}

void drawInterpolation(double x1, double x2, double x3, double h1, double h2, double h3, double v, Model m, vector<double> point){
    double y, y_min = 1000000, x_min;
    vector<double> xs, ys;
    for (double i = 0; i <= 1.0000001; i += 0.01) {
        xs.push_back(i);

        y =     ((i-x2)*(i-x3))/((x1-x2)*(x1-x3))*h1+
                ((i-x1)*(i-x3))/((x2-x1)*(x2-x3))*h2+
                ((i-x1)*(i-x2))/((x3-x1)*(x3-x2))*h3;

        ys.push_back(y);
        if(y < y_min) {
            y_min = y;
            x_min = i;
        }
    }
    cout << x_min << endl;
    plt::plot(xs,ys);
    plt::show();
    y_min = 100000;
    xs.clear(), ys.clear();
    for (double i = 0; i <= 1.0000001; i += 0.01) {
        xs.push_back(i);

        y = m.squaredTopDist(i,v,point);

        ys.push_back(y);
        if(y < y_min) {
            y_min = y;
            x_min = i;
        }
    }
    cout << x_min << endl;
    plt::plot(xs,ys);
    plt::show();


}

vector<double> quadInterpolation2(Model m, vector<double> P, int iter){
    double u1 = 0.3, u2 = 0.6, u3 = 0.9, v1 = 0, v2 = 0.5, v3 = 1;
    double h1,h2,h3,h4;
    double unew = 0.5, vnew = 0.5;
    for (int i = 0; i < iter; ++i) {

        if(unew != 0 && unew != 1) {

            h1 = m.squaredTopDist(u1, vnew, P);
            h2 = m.squaredTopDist(u2, vnew, P);
            h3 = m.squaredTopDist(u3, vnew, P);
            cout << u1 << " " << u2 << " " << u3 << " " << h1 << " " << h2 << " " << h3 << " " << endl;

            drawInterpolation(u1,u2,u3,h1,h2,h3,vnew,m,P);
            cout << u1 << " " << u2 << " " << u3 << " " << h1 << " " << h2 << " " << h3 << " " << endl;
            unew = (u1*u1*(h2-h3) + u2*u2*(h3-h1) + u3*u3*(h1-h2))/(2.0*(u1*(h2-h3) + u2*(h3-h1) + u3*(h1-h2)));
//            unew = 0.5 * ((u2 * u2 - u3 * u3) * h1 + (u3 * u3 - u1 * u1) * h2 + (u1 * u1 - u2 * u2) * h3) /
//                          ((u2 - u3) * h1 + (u3 - u1) * h2 + (u1 - u2) * h3);
            cout << unew << endl;
            if (unew < 0) {
                unew = 0;
            } else if (unew > 1) {
                unew = 1;
            }

            h4 = m.squaredTopDist(unew, vnew, P);

            double hmax = h1;
            if(h2 > hmax)
                hmax = h2;
            if(h3 > hmax)
                hmax = h3;
            if(h4 > hmax)
                hmax = h4;
            if(hmax == h1)
                u1 = unew;
            if(hmax == h2)
                u2 = unew;
            if(hmax == h3)
                u3 = unew;
            if(hmax == h4)
                unew = unew;


        }

        cout << "new u: " << unew << " new v: " << unew << endl;

        if(vnew != 0 && vnew != 1) {

            h1 = m.squaredTopDist(unew, v1, P);
            h2 = m.squaredTopDist(unew, v2, P);
            h3 = m.squaredTopDist(unew, v3, P);

//            drawInterpolation(v1,v2,v3,h1,h2,h3);
//
            vnew = 0.5 * ((v2 * v2 - v3 * v3) * h1 + (v3 * v3 - v1 * v1) * h2 + (v1 * v1 - v2 * v2) * h3) /
                   ((v2 - v3) * h1 + (v3 - v1) * h2 + (v1 - v2) * h3);

            if (vnew < 0) {
                vnew = 0;
            } else if (vnew > 1) {
                vnew = 1;
            }

            h4 = m.squaredTopDist(unew, vnew, P);

            double hmax = h1;
            if(h2 > hmax)
                hmax = h2;
            if(h3 > hmax)
                hmax = h3;
            if(h4 > hmax)
                hmax = h4;
            if(hmax == h1)
                v1 = vnew;
            if(hmax == h2)
                v2 = vnew;
            if(hmax == h3)
                v3 = vnew;
            if(hmax == h4)
                vnew = vnew;

        }

        cout << "new u: " << unew << " new v: " << unew << endl;
    }
    return {unew,unew};
}


double B(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, int dim){
    return (1-t)*(1-t)*Q0[dim] + 2*(1-t)*t*Q1[dim] + t*t*Q2[dim];
}

double B1(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, int dim){
    return -2*(1-t)*Q0[dim] + 2*(1-t)*Q1[dim] - 2*t*Q1[dim] + 2*t*Q2[dim];
}

double B2(vector<double> Q0, vector<double> Q1, vector<double> Q2, int dim){
    return 2*Q0[dim] - 4*Q1[dim] + 2*Q2[dim];
}


double C(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, vector<double> Q3, int dim){
    return (1-t)*B(t,Q0,Q1,Q2,dim) + t*B(t,Q1,Q2,Q3,dim);
}

double C1(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, vector<double> Q3, int dim){
    return 3*(1-t)*(1-t)*(Q1[dim] - Q0[dim]) + 6*(1-t)*t*(Q2[dim] - Q1[dim]) - 3*t*t*(Q3[dim] - Q2[dim]) ;
}

double C2(double t, vector<double> Q0, vector<double> Q1, vector<double> Q2, vector<double> Q3, int dim){
    return 6*(1-t)*(Q2[dim] - 2*Q1[dim] + Q0[dim]) + 6*t*(Q3[dim] - 2*Q2[dim] + Q1[dim]);
}

void drawB(vector<double> Q0, vector<double> Q1, vector<double> Q2){
    vector<double> xs, ys;
    for(double i = 0; i <= 1.0; i+= 0.01){
        xs.push_back(B(i,Q0,Q1,Q2,0));
        ys.push_back(B(i,Q0,Q1,Q2,1));
    }
    plt::plot(xs,ys);
}
void drawC(vector<double> Q0, vector<double> Q1, vector<double> Q2, vector<double> Q3){
    vector<double> xs, ys;
    for(double i = 0; i <= 1.0; i+= 0.01){
        xs.push_back(C(i,Q0,Q1,Q2,Q3,0));
        ys.push_back(C(i,Q0,Q1,Q2,Q3,1));
    }
    plt::plot(xs,ys);
}
void drawLine(vector<double> P0, vector<double> P1){
    vector<double> xs, ys;
    for(double i = 0; i <= 1.0; i+= 0.01){
        xs.push_back(P0[0] + (P1[0] - P0[0])*i);
        ys.push_back(P0[1] + (P1[1] - P0[1])*i);
    }
    plt::plot(xs,ys);
}
vector<double> add(vector<double> P, vector<double> Q){
    assert(P.size() == Q.size());
    for (int i = 0; i < P.size(); ++i) {
        P[i] += Q[i];
    }
    return P;
}
vector<double> sub(vector<double> P, vector<double> Q){
    assert(P.size() == Q.size());
    for (int i = 0; i < P.size(); ++i) {
        P[i] -= Q[i];
    }
    return P;
}
vector<double> mul(vector<double> P, double q){
    for (int i = 0; i < P.size(); ++i) {
        P[i] *= q;
    }
    return P;
}

vector<double> div(vector<double> P, double q){
    for (int i = 0; i < P.size(); ++i) {
        P[i] /= q;
    }
    return P;
}

double det(vector<double> P, vector<double> Q){
    return P[0]*Q[1] - P[1]*Q[0];
}

double sign(vector<double> P, vector<double> Q) {
    if(P[0]*Q[0] + P[1]*Q[1] > 0){
        return 1.0;
    }else{
        return -1.0;
    }
}

double magnitude(vector<double> v){
    double sum = 0;
    for (int i = 0; i < v.size(); ++i) {
        sum += v[i]*v[i];
    }
    return sqrt(sum);
}

double curvature(vector<double> c1, vector<double> c2){
    double mag = magnitude(c1);
    return det(c1,c2)/(mag*mag*mag);
}

vector<double> circleCenterDir(vector<double> c1, vector<double> c2){
    double sign = 0;
    if(det(c1,c2) > 0){
        sign = 1;
    }else{
        sign = -1;
    }
    return div({c1[1] * (-sign),c1[0] * sign},magnitude(c1));
}

vector<double> orthogonalProjectionTest(){

    vector<double> Q0 = {0,0}, Q1 = {1,0}, Q2 = {1,1}, Q3 = {2,2};
    vector<double> P = {0.5,0.1}, Q = {0,0}, F = {0,0}, F1 = {0,0}, F2 = {0,0};
    double t = 0.1, dt = 1;
    int iter = 10;
    int i = 0;
    double eps = 0.0000001;

    plt::xlim(-1.0,2.0);
    plt::ylim(-1.0,2.0);


    while(abs(dt) > eps && i < iter){
        F = {C(t,Q0,Q1,Q2,Q3,0), C(t,Q0,Q1,Q2,Q3,1)};
        F1 = {C1(t,Q0,Q1,Q2,Q3,0), C1(t,Q0,Q1,Q2,Q3,1)};
        F2 = {C2(t,Q0,Q1,Q2,Q3,0), C2(t,Q0,Q1,Q2,Q3,1)};

        double k = curvature(F1,F2);
        vector<double> c = circleCenterDir(F1,F2);
        vector<double> s = add(F,div(c,k));
        vector<double> u = sub(P,s);
        Q = add(s, div(u,k*magnitude(u)));


        drawC(Q0,Q1,Q2,Q3);
        drawLine(F,P);
        drawLine(F,s);
        drawLine(F,add(F,F2));
        plt::show();

        dt = sign(F1,sub(Q,F))*sqrt(abs( 2*det(F1,sub(Q,F))/det(F1,F2)));

        t += dt;
        i++;
        cout << "dt: " << dt << " iter: " << i << endl;
    }

    drawC(Q0,Q1,Q2,Q3);
    drawLine(F,P);
    plt::show();
    return {B(t,Q0,Q1,Q2,0), B(t,Q0,Q1,Q2,1)};
}

int factorial(int n){
    int f = 1;
    for (int i = 1; i <= n; ++i) {
        f *= i;
    }
    return f;
}
int binom(int n, int i){
    if(n < 2){
        return 1;
    }
    return factorial(n)/(factorial(i)*factorial(n-i));
}

vector<double> BSur(double u, double v, vector<vector<double>> k, int n, int m){
    vector<double> sum = {0,0,0};

    for (int i = 0; i < n+1; ++i) {
        for (int j = 0; j < m+1; ++j) {
            double param = binom(n,i)*pow(u,i)*pow(1-u,n-i)*binom(m,j)*pow(v,j)*pow(1-v,m-j);
            sum[0] += param*k[i*(n+1) + j][0];
            sum[1] += param*k[i*(n+1) + j][1];
            sum[2] += param*k[i*(n+1) + j][2];
        }
    }
    return sum;
}

vector<double> BSurDu(double u, double v, vector<vector<double>> k, int n, int m){
    vector<double> sum = {0,0,0};

    for (int i = 0; i <= n-1; ++i) {
        for (int j = 0; j <= m; ++j) {
            double param = binom(m, j) * pow(v, j) * pow(1 - v, m - j) *
                           binom(n - 1, i) * pow(u, i) * pow(1 - u, n - 1 - i);
            sum[0] += param * n * (k[(i + 1) * (n+1) + j][0] - k[i * (n+1) + j][0]);
            sum[1] += param * n * (k[(i + 1) * (n+1) + j][1] - k[i * (n+1) + j][1]);
            sum[2] += param * n * (k[(i + 1) * (n+1) + j][2] - k[i * (n+1) + j][2]);
        }
    }
    return sum;
}

vector<double> BSurDv(double u, double v, vector<vector<double>> k, int n, int m){
    vector<double> sum = {0,0,0};

    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m-1; ++j) {
            double param = binom(m-1,j)*pow(v,j)*pow(1-v,m-1-j)*
                           binom(n,i)*pow(u,i)*pow(1-u,n-i);
            sum[0] += param*m*(k[i*(n+1) + j+1][0] - k[i*(n+1) + j][0]);
            sum[1] += param*m*(k[i*(n+1) + j+1][1] - k[i*(n+1) + j][1]);
            sum[2] += param*m*(k[i*(n+1) + j+1][2] - k[i*(n+1) + j][2]);
        }
    }
    return sum;
}

vector<double> BSurDuu(double u, double v, vector<vector<double>> k, int n, int m){
    vector<double> sum = {0,0,0};

    for (int i = 0; i <= n-2; ++i) {
        for (int j = 0; j <= m; ++j) {
            double param = binom(m, j) * pow(v, j) * pow(1 - v, m - j) *
                           binom(n - 2, i) * pow(u, i) * pow(1 - u, n - 2 - i);
            sum[0] += param * n * (n-1) * (k[(i + 2) * (n+1) + j][0] - 2*k[(i + 1) * (n+1) + j][0] + k[i * (n+1) + j][0]);
            sum[1] += param * n * (n-1) * (k[(i + 2) * (n+1) + j][1] - 2*k[(i + 1) * (n+1) + j][1] + k[i * (n+1) + j][1]);
            sum[2] += param * n * (n-1) * (k[(i + 2) * (n+1) + j][2] - 2*k[(i + 1) * (n+1) + j][2] + k[i * (n+1) + j][2]);
        }
    }
    return sum;
}

vector<double> BSurDuv(double u, double v, vector<vector<double>> k, int n, int m) {
    vector<double> sum = {0, 0, 0};

    for (int i = 0; i <= n - 1; ++i) {
        for (int j = 0; j <= m - 1; ++j) {
            double param = binom(m - 1, j) * pow(v, j) * pow(1 - v, m - 1 - j) *
                           binom(n - 1, i) * pow(u, i) * pow(1 - u, n - 1 - i);
            sum[0] += param * n * m * (k[(i + 1) * (n+1) + j + 1][0] - k[(i + 1) * (n+1) + j][0] - k[i * (n+1) + j + 1][0] + k[i * (n+1) + j][0]);
            sum[1] += param * n * m * (k[(i + 1) * (n+1) + j + 1][1] - k[(i + 1) * (n+1) + j][1] - k[i * (n+1) + j + 1][1] + k[i * (n+1) + j][1]);
            sum[2] += param * n * m * (k[(i + 1) * (n+1) + j + 1][2] - k[(i + 1) * (n+1) + j][2] - k[i * (n+1) + j + 1][2] + k[i * (n+1) + j][2]);
        }
    }
    return sum;
}

vector<double> BSurDvv(double u, double v, vector<vector<double>> k, int n, int m){
    vector<double> sum = {0,0,0};

    for (int i = 0; i <= n-2; ++i) {
        for (int j = 0; j <= m; ++j) {
            double param = binom(m - 2, j) * pow(v, j) * pow(1 - v, m - 2 - j) *
                           binom(n, i) * pow(u, i) * pow(1 - u, n - i);
            sum[0] += param * m * (m-1) * (k[i * (n+1) + j + 2][0] - 2*k[i * (n+1) + j + 1][0] + k[i * (n+1) + j][0]);
            sum[1] += param * m * (m-1) * (k[i * (n+1) + j + 2][1] - 2*k[i * (n+1) + j + 1][1] + k[i * (n+1) + j][1]);
            sum[2] += param * m * (m-1) * (k[i * (n+1) + j + 2][2] - 2*k[i * (n+1) + j + 1][2] + k[i * (n+1) + j][2]);
        }
    }
    return sum;
}

vector<double> BSurParam(double u0, double v0, double t, vector<vector<double>> k, int n, int m, double lamu, double lamv){
    return BSur(u0+lamu*t,v0+lamv*t,k,n,m);
}
vector<double> BSurParamDt(double u0, double v0, double t, vector<vector<double>> k, int n, int m, double lamu, double lamv){
    return add(mul(BSurDu(u0+lamu*t,v0+lamv*t,k,n,m),lamu),mul(BSurDv(u0+lamu*t,v0+lamv*t,k,n,m),lamv));
}
vector<double> BSurParamDtt(double u0, double v0, double t, vector<vector<double>> k, int n, int m, double lamu, double lamv){
    return add(add(mul(BSurDuu(u0+lamu*t,v0+lamv*t,k,n,m),lamu*lamu),mul(BSurDuv(u0+lamu*t,v0+lamv*t,k,n,m),2*lamu*lamv)),mul(BSurDvv(u0+lamu*t,v0+lamv*t,k,n,m),lamv*lamv));
}

void drawSur(vector<vector<double>> C, int n, int m){
    vector<vector<double>> xs, ys, zs;
    for (double u = 0; u <= 1.0; u+=0.1) {
        vector<double> xr, yr, zr;
        for (double v = 0; v <= 1.0; v+=0.1) {
            vector<double> p = BSur(u,v,C,n,m);
            xr.push_back(p[0]);
            yr.push_back(p[1]);
            zr.push_back(p[2]);
        }
        xs.push_back(xr);
        ys.push_back(yr);
        zs.push_back(zr);
    }
    plt::plot_surface(xs,ys,zs);
}

void drawSurWithLine(vector<vector<double>> C, int n, int m, vector<double> line){
    assert(line.size() == 6);
    vector<vector<double>> xs, ys, zs;
    for (double u = 0; u <= 1.0; u+=0.1) {
        vector<double> xr, yr, zr;
        for (double v = 0; v <= 1.0; v+=0.1) {
            vector<double> p = BSur(u,v,C,n,m);
            xr.push_back(p[0]);
            yr.push_back(p[1]);
            zr.push_back(p[2]);
        }
        xs.push_back(xr);
        ys.push_back(yr);
        zs.push_back(zr);
    }
    double len = sqrt((line[3])*(line[3]) + (line[4])*(line[4]) + (line[5])*(line[5]));
    plt::plot_surface_with_line_and_axes_set(xs,ys,zs,{-250,50,-200,150,-30,40},{line[0]},{line[1]},{line[2]},{line[3]},{line[4]},{line[5]},len);
}

void drawFillTopWithLine(Model m, vector<double> line, vector<double> axes){
    assert(line.size() == 6);
    vector<vector<double>> xs, ys, zs;
    for (double u = 0; u <= 1.0; u+=0.1) {
        vector<double> xr, yr, zr;
        for (double v = 0; v <= 1.0; v+=0.1) {
            vector<double> p = m.fillTop(u,v);
            xr.push_back(p[0]);
            yr.push_back(p[1]);
            zr.push_back(p[2]);
        }
        xs.push_back(xr);
        ys.push_back(yr);
        zs.push_back(zr);
    }
    double len = sqrt((line[3])*(line[3]) + (line[4])*(line[4]) + (line[5])*(line[5]));
    plt::plot_surface_with_line_and_axes_set(xs,ys,zs,axes,{line[0]},{line[1]},{line[2]},{line[3]},{line[4]},{line[5]},len);
}

double dot(vector<double> a, vector<double> b){
    double sum = 0;
    assert(a.size() == b.size());
    for (int i = 0; i < a.size(); ++i) {
        sum += a[i]*b[i];
    }
    return sum;
}

vector<double> cross(vector<double> a, vector<double> b){
    assert(a.size() == 3);
    assert(b.size() == 3);
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}

vector<double> solve3Dlinear(vector<double> x, vector<double> y, vector<double> z){
    assert(x.size() == 4);
    assert(y.size() == 4);
    assert(z.size() == 4);


//    cout << " Initial linear system: " << endl;
//    cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << endl;
//    cout << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << endl;
//    cout << z[0] << " " << z[1] << " " << z[2] << " " << z[3] << " " << endl;
//    cout << " ---------------------- " << endl;


    if(x[2] == 0){
        if(z[2] != 0){
            vector<double> temp = x;
            x = z;
            z = temp;
        }else{
            vector<double> temp = x;
            x = y;
            y = temp;
        }
    }
    if(y[1] == 0){
        if(z[1] != 0){
            vector<double> temp = y;
            y = z;
            z = temp;
        }else{
            vector<double> temp = y;
            y = x;
            x = temp;
        }
    }
    if(z[0] == 0){
        if(y[0] != 0){
            vector<double> temp = y;
            y = z;
            z = temp;
        }else{
            vector<double> temp = z;
            z = x;
            x = temp;
        }
    }
    z = sub(z,mul(x,z[2]/x[2]));
    y = sub(y,mul(x,y[2]/x[2]));
    z = sub(z,mul(y,z[1]/y[1]));
    y = sub(y,mul(z,y[0]/z[0]));
    x = sub(x,mul(z,x[0]/z[0]));
    x = sub(x,mul(y,x[1]/y[1]));

//    cout << " Solved linear system: " << endl;
//    cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << endl;
//    cout << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << endl;
//    cout << z[0] << " " << z[1] << " " << z[2] << " " << z[3] << " " << endl;
//    cout << " ---------------------- " << endl;

    return {z[3]/z[0],y[3]/y[1],x[3]/x[2]};

}
double sign3D(vector<double> P, vector<double> Q) {
    if(dot(P,Q) > 0){
        return 1.0;
    }else{
        return -1.0;
    }
}

vector<double> orthogonalProjectionTestSur3D(){

    vector<double> axes = {-250,50,-200,150,-30,40};

    vector<vector<double>> C = {{-236,-197,-22},{-206,-117,-22},{-216,-27,8},{-246,62,-22},
                                {-156,-177,8},{-176,-97,38},{-157,20,126},{-186,142,8},
                                {-86,-157,8},{-138,-113,-146},{-104,14,-60},{-96,-102,8},
                                {-6,-197,-22},{-47,-96,-33},{25,32,95},{-6,102,8}};
    vector<vector<double>> S = {{0,0,0},{1,0,0},{2,0,0},
                                {0,1,0},{1,1,10},{2,1,0},
                                {0,2,0},{1,2,0},{2,2,0}};
    vector<double> P = {-175,0,20}, Q = {0,0}, F = {0,0}, F1 = {0,0}, F2 = {0,0};
    double u = 0.2, v = 0.6;
    double lam1 = 1, lam2 = 2;
    double dt = 1;
    int iter = 1000;
    int i = 0;
    double eps = 0.0000001;

    cout << C[4][2] << endl;
    drawSur(C,3,3);
    plt::show();


    while(abs(dt*lam1 + dt*lam2) > eps && i < iter){
        vector<double> p0 = BSur(u,v,C,3,3),p0P;
        vector<double> x,y,z,du,duu,duv,dvv,dv,n;

        p0P = sub(P,p0);

//        cout << "u: " << u << " v: " << v << endl;
        drawSurWithLine(C,3,3, {p0[0],p0[1],p0[2],p0P[0],p0P[1],p0P[2]});
        plt::show();
        du = BSurDu(u,v,C,3,3);
        dv = BSurDv(u,v,C,3,3);
        duu = BSurDuu(u,v,C,3,3);
        duv = BSurDuv(u,v,C,3,3);
        dvv = BSurDvv(u,v,C,3,3);

        n = cross(du,dv);
        n = div(n,magnitude(n));



        vector<double> coefs = solve3Dlinear({du[0],dv[0],n[0],p0P[0]},{du[1],dv[1],n[1],p0P[1]},{du[2],dv[2],n[2],p0P[2]});
        lam1 = coefs[0];
        lam2 = coefs[1];
        double nu = coefs[2];
        cout << "lam1: " << lam1 << " lam2: " << lam2 << " nu: " << nu << endl;




        double g11, g12, g22, h11, h12, h22, g21, h21;

        g11 = dot(du,du);
        g12 = dot(du,dv);
        g21 = dot(dv,du);
        g22 = dot(dv,dv);
        h11 = dot(duu,n);
        h12 = dot(duv,n);
        h21 = dot(duv,n);
        h22 = dot(dvv,n);

        double k = (h11*lam1*lam1 + h12*lam1*lam2 + h22*lam2*lam2 + h21*lam2*lam1)/(g11*lam1*lam1 + g12*lam1*lam2 + g22*lam2*lam2 + g21*lam2*lam1);

//        cout << "curvature: " << k << " radius: " << 1.0/k << endl;

        drawSurWithLine(C,3,3, {p0[0],p0[1],p0[2],n[0]/k,n[1]/k,n[2]/k});
        plt::show();

        vector<double> s = add(p0,mul(n,1.0/k));

//        drawSurWithLine(C,2,2, {0,0,0,s[0],s[1],s[2]});
//        plt::show();

        vector<double> Ps = sub(P,s);
        vector<double> q = add(s,div(Ps,magnitude(Ps)*abs(k)));
        double drawMag = magnitude(Ps)*abs(k);
        drawSurWithLine(C,3,3, {s[0],s[1],s[2],Ps[0]/drawMag,Ps[1]/drawMag,Ps[2]/drawMag});
        plt::show();

        vector<double> c = BSurParam(u,v,0,C,3,3,lam1,lam2);
        vector<double> c1 = BSurParamDt(u,v,0,C,3,3,lam1,lam2);
        vector<double> c2 = BSurParamDtt(u,v,0,C,3,3,lam1,lam2);


        dt = sign3D(c1,sub(q,c))*sqrt(abs(magnitude(cross(c1,sub(q,c)))/magnitude(cross(c1,c2))));

        double dist = magnitude(p0P);

        u += dt*lam1;
        v += dt*lam2;
        i++;

        cout << "dt: " << dt*lam1 + dt*lam2 << " ddist: " << dist - magnitude(sub(P,BSur(u,v,C,3,3))) << endl;
    }

}


vector<double> fillTopParam(Model m, double u0, double v0, double t, double lamu, double lamv){
    return m.fillTop(u0+lamu*t,v0+lamv*t);
}

vector<double> fillTopParamDt(Model m, double u0, double v0, double t, double lamu, double lamv){
    return add(mul(m.fillTopDerU(u0+lamu*t,v0+lamv*t),lamu),mul(m.fillTopDerV(u0+lamu*t,v0+lamv*t),lamv));
}

vector<double> fillTopParamDtt(Model m, double u0, double v0, double t, double lamu, double lamv){
    return add(add(mul(m.fillTopDerUU(u0+lamu*t,v0+lamv*t),lamu*lamu),mul(m.fillTopDerUV(u0+lamu*t,v0+lamv*t),2*lamu*lamv)),mul(m.fillTopDerVV(u0+lamu*t,v0+lamv*t),lamv*lamv));
}


vector<double> orthogonalProjection1DU(Model m, vector<double> P, double u, double v, double eps, int iter) {
    double u0 = u, v0 = v;
    double dist = 100000, distOld;
    double lam1 = 1, lam2 = 2;
    double dt = 1;
    int i = 0;

    vector<double> p0,p0P,du,duu,duv,dvv,dv,n,coefs,Ps,q,c1,c2,cq,s;
    double g11, g12, g22, h11, h12, h22, g21, h21,k;

    p0 = m.fillTop(u, v);

    p0P = sub(P, p0);
    distOld = dist;
    dist = magnitude(p0P);
    do {
//        cout << "u: " << u << " v: " << v << endl;
//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],p0P[0],p0P[1],p0P[2]}, {0,15,0,15,5,20});
//        plt::show();
        du = m.fillTopDerU(u, v);
        dv = m.fillTopDerV(u, v);
        duu = m.fillTopDerUU(u, v);
        duv = m.fillTopDerVV(u, v);
        dvv = m.fillTopDerUV(u, v);

        n = cross(du, dv);
        n = div(n, magnitude(n));

//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],du[0],du[1],du[2]}, {0,15,0,15,5,20});
//        plt::show();
//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],dv[0],dv[1],dv[2]}, {0,15,0,15,5,20});
//        plt::show();
//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],n[0],n[1],n[2]}, {0,15,0,15,5,20});
//        plt::show();
        coefs = solve3Dlinear({du[0], dv[0], n[0], p0P[0]}, {du[1], dv[1], n[1], p0P[1]}, {du[2], dv[2], n[2], p0P[2]});
        lam1 = coefs[0];
        lam2 = coefs[1];
        //double nu = coefs[2];

//        cout << "lam1: " << lam1 << " lam2: " << lam2 << " nu: " << nu << endl;


        g11 = dot(du, du);
        g12 = dot(du, dv);
        g21 = dot(dv, du);
        g22 = dot(dv, dv);
        h11 = dot(duu, n);
        h12 = dot(duv, n);
        h21 = h12;
        h22 = dot(dvv, n);

        k = (h11 * lam1 * lam1 + h12 * lam1 * lam2 + h21 * lam2 * lam1 + h22 * lam2 * lam2) /
            (g11 * lam1 * lam1 + g12 * lam1 * lam2 + g21 * lam2 * lam1 + g22 * lam2 * lam2);

//        cout << "curvature: " << k << " radius: " << 1.0/k << endl;

//        drawSurWithLine(C,2,2, {p0[0],p0[1],p0[2],n[0]/k,n[1]/k,n[2]/k});
//        plt::show();

        s = add(p0, mul(n, 1.0 / k));
//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],n[0]/k,n[1]/k,n[2]/k}, {0,15,0,15,5,20});
//        plt::show();
//        drawSurWithLine(C,2,2, {0,0,0,s[0],s[1],s[2]});
//        plt::show();

        Ps = sub(P, s);
        q = add(s, div(Ps, magnitude(Ps) * abs(k)));
        //double drawMag = magnitude(Ps)*abs(k);

//        drawFillTopWithLine(m, {s[0],s[1],s[2],Ps[0]/drawMag,Ps[1]/drawMag,Ps[2]/drawMag}, {0,15,0,15,5,20});
//        plt::show();

        c1 = add(mul(du, lam1), mul(dv, lam2));
        c2 = add(add(mul(duu, lam1 * lam1), mul(duv, 2 * lam1 * lam2)), mul(dvv, lam2 * lam2));

//        vector<double> c = fillTopParam(m,u,v,0,lam1,lam2);
//        vector<double> c1 = fillTopParamDt(m,u,v,0,lam1,lam2);
//        vector<double> c2 = fillTopParamDtt(m,u,v,0,lam1,lam2);

        cq = sub(q, p0);
        dt = sign3D(c1, cq) * sqrt(abs(magnitude(cross(c1, cq)) / magnitude(cross(c1, c2))));

        u += dt * lam1;
        v += dt * lam2;
        i++;
        if (u < 0) {

        }
        if (u > 1) {

        }
        if (v < 0) {

        }
        if (v > 1) {

        }

        p0 = m.fillTop(u, v);

        p0P = sub(P, p0);
        distOld = dist;
        dist = magnitude(p0P);
        //cout << "dt: " << dt*lam1 + dt*lam2 << " ddist: " << abs(dist - magnitude(sub(P,m.fillTop(u,v)))) << endl;
    }while(abs(dist - distOld) > eps && i < iter);

    return {u,v,magnitude(sub(P,m.fillTop(u,v))),u0,v0,static_cast<double>(i)};
}

vector<double> orthogonalProjection1DV(Model m, vector<double> P, double u, double v, double eps, int iter) {

}

vector<double> orthogonalProjection(Model m, vector<double> P, double u, double v, double eps = 0.0000001, int iter = 100){

    double u0 = u, v0 = v;
    double dist = 100000, distOld;
    double lam1 = 1, lam2 = 2;
    double dt = 1;
    int i = 0;

    vector<double> p0,p0P,du,duu,duv,dvv,dv,n,coefs,Ps,q,c1,c2,cq,s;
    double g11, g12, g22, h11, h12, h22, g21, h21,k;

    p0 = m.fillTop(u, v);

    p0P = sub(P, p0);
    distOld = dist;
    dist = magnitude(p0P);
    do {
//        cout << "u: " << u << " v: " << v << endl;
//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],p0P[0],p0P[1],p0P[2]}, {0,15,0,15,5,20});
//        plt::show();
        du = m.fillTopDerU(u, v);
        dv = m.fillTopDerV(u, v);
        duu = m.fillTopDerUU(u, v);
        duv = m.fillTopDerVV(u, v);
        dvv = m.fillTopDerUV(u, v);

        n = cross(du, dv);
        n = div(n, magnitude(n));

//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],du[0],du[1],du[2]}, {0,15,0,15,5,20});
//        plt::show();
//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],dv[0],dv[1],dv[2]}, {0,15,0,15,5,20});
//        plt::show();
//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],n[0],n[1],n[2]}, {0,15,0,15,5,20});
//        plt::show();
        coefs = solve3Dlinear({du[0], dv[0], n[0], p0P[0]}, {du[1], dv[1], n[1], p0P[1]}, {du[2], dv[2], n[2], p0P[2]});
        lam1 = coefs[0];
        lam2 = coefs[1];
        //double nu = coefs[2];

//        cout << "lam1: " << lam1 << " lam2: " << lam2 << " nu: " << nu << endl;


        g11 = dot(du, du);
        g12 = dot(du, dv);
        g21 = dot(dv, du);
        g22 = dot(dv, dv);
        h11 = dot(duu, n);
        h12 = dot(duv, n);
        h21 = h12;
        h22 = dot(dvv, n);

        k = (h11 * lam1 * lam1 + h12 * lam1 * lam2 + h21 * lam2 * lam1 + h22 * lam2 * lam2) /
            (g11 * lam1 * lam1 + g12 * lam1 * lam2 + g21 * lam2 * lam1 + g22 * lam2 * lam2);

//        cout << "curvature: " << k << " radius: " << 1.0/k << endl;

//        drawSurWithLine(C,2,2, {p0[0],p0[1],p0[2],n[0]/k,n[1]/k,n[2]/k});
//        plt::show();

        s = add(p0, mul(n, 1.0 / k));
//        drawFillTopWithLine(m, {p0[0],p0[1],p0[2],n[0]/k,n[1]/k,n[2]/k}, {0,15,0,15,5,20});
//        plt::show();
//        drawSurWithLine(C,2,2, {0,0,0,s[0],s[1],s[2]});
//        plt::show();

        Ps = sub(P, s);
        q = add(s, div(Ps, magnitude(Ps) * abs(k)));
        //double drawMag = magnitude(Ps)*abs(k);

//        drawFillTopWithLine(m, {s[0],s[1],s[2],Ps[0]/drawMag,Ps[1]/drawMag,Ps[2]/drawMag}, {0,15,0,15,5,20});
//        plt::show();

        c1 = add(mul(du, lam1), mul(dv, lam2));
        c2 = add(add(mul(duu, lam1 * lam1), mul(duv, 2 * lam1 * lam2)), mul(dvv, lam2 * lam2));

//        vector<double> c = fillTopParam(m,u,v,0,lam1,lam2);
//        vector<double> c1 = fillTopParamDt(m,u,v,0,lam1,lam2);
//        vector<double> c2 = fillTopParamDtt(m,u,v,0,lam1,lam2);

        cq = sub(q, p0);
        dt = sign3D(c1, cq) * sqrt(abs(magnitude(cross(c1, cq)) / magnitude(cross(c1, c2))));

        u += dt * lam1;
        v += dt * lam2;
        i++;
        if (u < 0) {

        }
        if (u > 1) {

        }
        if (v < 0) {

        }
        if (v > 1) {

        }

        p0 = m.fillTop(u, v);

        p0P = sub(P, p0);
        distOld = dist;
        dist = magnitude(p0P);
        //cout << "dt: " << dt*lam1 + dt*lam2 << " ddist: " << abs(dist - magnitude(sub(P,m.fillTop(u,v)))) << endl;
    }while(abs(dist - distOld) > eps && i < iter);

    return {u,v,magnitude(sub(P,m.fillTop(u,v))),u0,v0,static_cast<double>(i)};
}

void plotSurface(vector<double> P, Model m, double ul = -.2, double vl = -.2, double ur = 1.2, double vr = 1.2, double dU = 0.1, double dV = 0.1){
    vector<vector<double>> u,v,dist;
    double min_dist = 100000;
    for (double i = ul-dU; i < ur+dU; i += (ur+dU - (ul-dU))/10) {
        vector<double> uRow, vRow, distRow;
        for (double j = vl-dV; j < vr+dV; j += (vr+dV - (vl-dV))/10) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(m.distanceToTopPoint(j,i,P) + 100*constraintFunc(j,i));
            if(min_dist > m.distanceToTopPoint(j,i,P) + 100*constraintFunc(j,i))
                min_dist = m.distanceToTopPoint(j,i,P) + 100*constraintFunc(j,i);
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    plt::plot_surface(u,v,dist);
}

double gridSearch(Model model, vector<double> A, double inc = 0.05){
    double min_dist = 10000;
    for (double i = 0; i < 1; i += inc) {
        for (double j = 0; j < 1; j += inc) {
            if(min_dist > model.distanceToTopPoint(j,i,A))
                min_dist = model.distanceToTopPoint(j,i,A);
        }
    }
    return min_dist;
}

vector<double> BinarySearch(Model model, vector<double> A, double eps = 0.0001, int iterations = 100){
    double ul = 0, ur = 1, vl = 0, vr = 1;
    double dist = model.distanceToTopPoint(0.5,0.5,A);
    double dist_old = DBL_MAX;
    double um, vm;
    int i = 0;
    while(abs(dist_old - dist) > eps && i < iterations){
        dist_old = dist;
        if( i % 2 == 0 ) {
            um = (ul+ur)/2;
            vm = (vl+vr)/2;
            double u1 = (ul+um)/2;
            double u2 = (ur+um)/2;
            double dist1 = model.distanceToTopPoint(u1,vm,A);
            double dist2 = model.distanceToTopPoint(u2,vm,A);
            if(dist1 < dist2){
                ur = um;
                dist = dist1;
            }else{
                ul = um;
                dist = dist2;
            }
        }else{
            um = (ul+ur)/2;
            vm = (vl+vr)/2;
            double v1 = (vl + vm)/2;
            double v2 = (vr + vm)/2;
            double dist1 = model.distanceToTopPoint(um,v1,A);
            double dist2 = model.distanceToTopPoint(um,v2,A);
            if(dist1 < dist2){
                vr = vm;
                dist = dist1;
            }else{
                vl = vm;
                dist = dist2;
            }
        }
        i++;
    }

    return {(ul + ur)/2,(vl + vr)/2,dist,0.5,0.5, static_cast<double>(i)};
}

void DrawDist(Model &model, const vector<double> &P) {
    vector<vector<double>> u, v, dist;
    double min_dist = 100000;
    for (double i = -0.5; i < 1.5; i += 0.15) {
        vector<double> uRow, vRow, distRow;
        for (double j = -0.5; j < 1.5; j += 0.15) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(model.distanceToTopPoint(j, i, P));
            if (min_dist > model.distanceToTopPoint(j, i, P)) {
                min_dist = model.distanceToTopPoint(j, i, P);
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    vector<double> x2, y2, z2, u2, v2, w2;
    for (double i = -0.5; i < 1.5; i += 0.15) {
        for (double j = -0.5; j < 1.5; j += 0.15) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.completeDistanceTopDerU(j, i, P);
            tV = model.completeDistanceTopDerV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);
    plt::show();
}

void TestAlgorithmOnRandomUniform(Model model, int algorithm, int preprocess, int seed, int plot = 0, vector<double> range = {4, 4, 4, 10, 10, 10}, int iterations = 1000, double GSAcc = 0.05, double eps = 0.000001){
    int reliability = 0;
    double avgError = 0;
    int sumIters = 0;
    srand(seed);
    for (int i = 0; i < iterations; i++) {
        vector<double> P = {
                (double(rand())/RAND_MAX)*(range[3]-range[0]) + range[0],
                (double(rand())/RAND_MAX)*(range[4]-range[1]) + range[1],
                (double(rand())/RAND_MAX)*(range[5]-range[2]) + range[2]
        };
        vector<double> start;
        double gridV = gridSearch(model,P,GSAcc);
        vector<double> algV;
        if(algorithm == 5)
            preprocess = 0;

        if(preprocess)
            start = findStartingPoint(model,P,preprocess);
        else
            start = {0.5,0.5};

        switch (algorithm) {
            case 0:
                algV = gradientDescent(model,P,start[0],start[1],eps);
                if(plot)
                    cout << "Gradient Descent ";
                break;
            case 1:
                algV = FletcherReevesCG(model,P,start[0],start[1],eps);
                if(plot)
                    cout << "Fletcher-Reeves CG ";
                break;
            case 2:
                algV = NewtonMethod(model,P, start[0], start[1], eps);
                if(plot)
                    cout << "Newton's Method ";
                break;
            case 3:
                algV = NewtonMethodProjected(model,P,start[0],start[1],eps,1,1600);
                if(plot)
                    cout << "Newton's Method constrained ";
                break;
            case 4:
                // TODO: UNFINISHED!!
                algV = FletcherReevesCGSecant(model,P,start[0],start[1],1.0/preprocess,eps);
                if(plot)
                    cout << "Fletcher-Reeves CG Secant constrained ";
                break;
            case 5:
                algV = BinarySearch(model,P);
                if(plot)
                    cout << "Binary Search ";
                break;
            default:
                cout << "Algorithm number out of bounds" << endl;
                i = iterations;
                break;
        }
        if(plot)
            cout << "u: " << algV[0] << " v: " << algV[1] << " V: " << algV[2] << " u0: " << algV[3] << " v0: " << algV[4] << " iterations: " << algV[5] << endl;

        if(algV[2] <= gridV) {
            reliability++;
        }else{
            avgError += (algV[2] - gridV)/gridV;
            if(plot) {
                cout << "Grid Search: " << gridV << endl;
                DrawDist(model, P);
            }
        }

        sumIters += int(algV[5]);

    }
    switch (algorithm) {
        case 0:
            cout << "Gradient Descent ";
            break;
        case 1:
            cout << "Fletcher-Reeves CG ";
            break;
        case 2:
            cout << "Newton's Method ";
            break;
        case 3:
            cout << "Newton's Method constrained ";
            break;
        case 4:
            cout << "Fletcher-Reeves CG Secant ";
            break;
        case 5:
            cout << "Binary Search ";
            break;
        default:
            break;
    }
    cout << "reliability: " << double(reliability)/iterations << " average error: " << avgError/iterations << " average iterations: " << double(sumIters)/iterations << endl;
}

void TestAlgorithmOnRandomCloseToSurface(Model model, int algorithm, int preprocess, int seed, int plot = 0, int iterations = 1000, double GSAcc = 0.05, double eps = 0.0000001){
    int reliability = 0;
    double avgError = 0;
    int sumIters = 0;
    srand(seed);

    time_t avgTime = 0;

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration1;
    chrono::microseconds t_duration2;


    for (int i = 0; i < iterations; i++) {
        vector<double> P = model.fillTop(double(rand())/RAND_MAX,double(rand())/RAND_MAX);

        P[0] += (double(rand())/RAND_MAX)*0.5;
        P[1] += (double(rand())/RAND_MAX)*0.5;
        P[2] += (double(rand())/RAND_MAX)*0.5;

        vector<double> start;
        double gridV = gridSearch(model,P,GSAcc);
        vector<double> algV;
        if(algorithm == 5)
            preprocess = 0;

        if(preprocess) {
            t_start = chrono::high_resolution_clock::now();
            start = findStartingPoint(model, P, preprocess);
            t_stop = chrono::high_resolution_clock::now();
            t_duration1 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);

        }else
            start = {0.5,0.5};

        switch (algorithm) {
            case 0:
                algV = gradientDescent(model,P,start[0],start[1],eps);
                if(plot)
                    cout << "Gradient Descent ";
                break;
            case 1:
                algV = FletcherReevesCG(model,P,start[0],start[1],eps);
                if(plot)
                    cout << "Fletcher-Reeves CG ";
                break;
            case 2:
                algV = NewtonMethod(model,P, start[0], start[1], eps);
                if(plot)
                    cout << "Newton's Method ";
                break;
            case 3:
                t_start = chrono::high_resolution_clock::now();
                algV = NewtonMethodProjected(model,P,start[0],start[1],eps,0.01,1600);
                t_stop = chrono::high_resolution_clock::now();
                t_duration2 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
                if(preprocess) {
                    avgTime += t_duration1.count() + t_duration2.count();
                    cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
                    cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
                }
                if(plot)
                    cout << "Newton's Method constrained ";
                break;
            case 4:
                // TODO: UNFINISHED!!
                algV = FletcherReevesCGSecant(model,P,start[0],start[1],1.0/preprocess,eps);
                if(plot)
                    cout << "Fletcher-Reeves CG Secant constrained ";
                break;
            case 5:
                algV = BinarySearch(model,P);
                if(plot)
                    cout << "Binary Search ";
                break;
            case 6:
                t_start = chrono::high_resolution_clock::now();
                algV = orthogonalProjection(model,P,start[0], start[1],eps);
                t_stop = chrono::high_resolution_clock::now();
                t_duration2 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
                if(preprocess) {
                    avgTime += t_duration1.count() + t_duration2.count();
                    cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
                    cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
                }
                if(plot){
                    cout << "Geometric search: ";
                }
                break;
            default:
                cout << "Algorithm number out of bounds" << endl;
                i = iterations;
                break;
        }
        if(plot)
            cout << "u: " << algV[0] << " v: " << algV[1] << " V: " << algV[2] << " u0: " << algV[3] << " v0: " << algV[4] << " iterations: " << algV[5] << endl;

        if(algV[2] <= gridV) {
            reliability++;
        }else{
            avgError += (algV[2] - gridV)/gridV;
            if(plot) {
                cout << "Grid Search: " << gridV << endl;
                vector<vector<double>> u, v, dist;
                double min_dist = 100000;
                for (double i = -0.5; i < 1.5; i += 0.15) {
                    vector<double> uRow, vRow, distRow;
                    for (double j = -0.5; j < 1.5; j += 0.15) {
                        uRow.push_back(j);
                        vRow.push_back(i);
                        distRow.push_back(model.squaredTopDist(j, i, P));
                        if (min_dist > model.squaredTopDist(j, i, P)) {
                            min_dist = model.squaredTopDist(j, i, P);
                        }
                    }
                    u.push_back(uRow);
                    v.push_back(vRow);
                    dist.push_back(distRow);
                }
                vector<double> x2, y2, z2, u2, v2, w2;
                for (double i = -0.5; i < 1.5; i += 0.15) {
                    for (double j = -0.5; j < 1.5; j += 0.15) {
                        x2.push_back(j);
                        y2.push_back(i);
                        z2.push_back(min_dist);
                        double tU, tV;
                        tU = model.squaredTopDistDerU(j, i, P);
                        tV = model.squaredTopDistDerV(j, i, P);
                        u2.push_back(tU);
                        v2.push_back(tV);
                        w2.push_back(0);
                    }
                }
                plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);
                plt::show();
            }
        }

        sumIters += int(algV[5]);

    }
    switch (algorithm) {
        case 0:
            cout << "Gradient Descent ";
            break;
        case 1:
            cout << "Fletcher-Reeves CG ";
            break;
        case 2:
            cout << "Newton's Method ";
            break;
        case 3:
            cout << "Newton's Method constrained ";
            break;
        case 4:
            cout << "Fletcher-Reeves CG Secant ";
            break;
        case 5:
            cout << "Binary Search ";
            break;
        case 6:
            cout << "Geometric Search ";
        default:
            break;
    }
    cout << "reliability: " << double(reliability)/iterations << " average error: " << avgError/iterations << " average iterations: " << double(sumIters)/iterations << " average time: " << avgTime/iterations << " microseconds " << endl;
}

int main() {


//    orthogonalProjectionTestSur3D();
//
//    vector<double> c = orthogonalProjectionTest();
//
//    cout << c[0] << " " << c[1] << endl;


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
    vector<double> point = p1.fillTop(0.7,0.6);
    point[0] += 0.1;
    point[1] += 0.1;
    point[2] += 0.4;
    //orthogonalProjection(p1,point,0.5,0.5);
    //TestAlgorithmOnRandomCloseToSurface(p1,6,0,0,1);

    vector<vector<double>> *Ps;
    Ps = static_cast<vector<vector<double>> *>(malloc(1000 * sizeof(vector<double>(3))));

    for (int i = 0; i < 1000; ++i) {
        Ps[i].push_back({1.0*i,2.0*i,3.0*i});
    }
    double a = Ps[999][2];
    cout << a << endl;
    free(Ps);
    /*double dist = p1.squaredTopDist(0.3,0.5,point);

    double dU = p1.squaredTopDistDerU(0.3,0.5, point);
    double dV = p1.squaredTopDistDerV(0.3,0.5, point);

    double dUU = p1.squaredTopDistDerUU(0.3,0.5, point);
    double dUV = p1.squaredTopDistDerUV(0.3,0.5, point);
    double dVV = p1.squaredTopDistDerVV(0.3,0.5, point);


    cout << dist << " " << dU << " " << dV << " " << dUU << " " << dUV << " " << dVV << " " << endl;

    quadInterpolation2(p1, point, 8);
    vector<double> res = BinarySearch(p1,point,0.0001,8);

    cout << "binary " << res[0] << " " << res[1] << endl;

    DrawDistWithDirection(p1,point);*/
    //TestAlgorithmOnRandomCloseToSurface(p1,1,12,0,1);

//    drawPart(p2, x, y, z);
//
//    drawPart(p3, x, y, z);
//    drawPart(p4, x, y, z);
//
//    drawPart(p5, x, y, z);
//    drawPart(p6, x, y, z);

    /*plt::plot_surface(x,y,z);
    plt::xlim(-20,20);
    plt::ylim(-20,20);

    plt::show();


    vector<double> A = {4,4,4};
    vector<vector<double>> u,v,dist;
    double min_dist = 100000;
    for (double i = -0.2; i < 1.2; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = -0.2; j < 1.2; j += 0.05) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(p1.distanceToTopPoint(j,i,A) + 100*constraintFunc(j,i));
            if(min_dist > p1.distanceToTopPoint(j,i,A) + 100*constraintFunc(j,i))
                min_dist = p1.distanceToTopPoint(j,i,A) + 100*constraintFunc(j,i);
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }

//    plt::plot_surface(u,v,dist);
//    plt::show();
    double scale = 0.10;
    vector<double> x2,y2,z2,u2,v2,w2;
    for (double i = -0.2; i < 1.2; i += 0.05) {
        for (double j = -0.2; j < 1.2; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = p1.completeDistanceTopDerU(j,i,A) + 100*constraintFuncdU(j,i);
            tV = p1.completeDistanceTopDerV(j,i,A) + 100*constraintFuncdV(j,i);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    cout << "Rough expected distance: " << min_dist << endl;
    vector<vector<double>> CGus, CGvs, CGVs;
    vector<vector<double>> Nus, Nvs, NVs;

    for (double i = 0; i < 1; i += 0.05) {
        vector<double> CGusR, CGvsR, CGVsR;
        vector<double> NusR, NvsR, NVsR;
        for (double j = 0; j < 1; j += 0.05) {
            vector<double> res;
            res = NewtonMethod(p1,A,j,i,0.00001,.5);
            NusR.push_back(res[3]);
            NvsR.push_back(res[4]);
            if(res[0] < 1 && res[0] > 0 && res[1] < 1 && res[1] > 0){
                NVsR.push_back(0);
                cout << "Newton u: " << res[0] << " v: " << res[1] << " V: " << res[2] << " u0: " << res[3] << " v0: " << res[4] << " iterations: " << res[5] << endl;

            }else {
                NVsR.push_back(1);
                cout << "Newton u: " << res[0] << " v: " << res[1] << " V: " << res[2] << " u0: " << res[3] << " v0: " << res[4] << " iterations: " << res[5] << endl;
            }

            res = FletcherReevesCG(p1,A,j,i,0.00001);
            CGusR.push_back(res[3]);
            CGvsR.push_back(res[4]);
            if(res[0] < 1 && res[0] > 0 && res[1] < 1 && res[1] > 0){
                CGVsR.push_back(0);
            }else {
                CGVsR.push_back(1);
            }
            if( res[0] < 1 && res[0] > 0 && res[1] < 1 && res[1] > 0 && res[2] > min_dist) {
                cout << "CG     u: " << res[0] << " v: " << res[1] << " V: " << res[2] << " u0: " << res[3] << " v0: "
                     << res[4] << endl;
                FletcherReevesCG(p1, A, j, i, 0.00001, 0);
            }
        }
        CGus.push_back(CGusR);
        CGvs.push_back(CGvsR);
        CGVs.push_back(CGVsR);
        Nus.push_back(NusR);
        Nvs.push_back(NvsR);
        NVs.push_back(NVsR);
    }
    plt::plot_surface(Nus,Nvs,NVs);
    plt::show();
    plt::plot_surface(CGus,CGvs,CGVs);
    plt::show();
    cout << "Distance via Newton's Method for 2D: " << NewtonMethod(p1,A,0.8,0.2,0.00001,.5)[2] << endl;
    plt::plot_surface_with_vector_field(u,v,dist,x2,y2,z2,u2,v2,w2,0.03);
    plt::show();
    cout << "Distance via grid search: " << min_dist << endl;
    cout << "Distance via Fletcher-Reeves CG: " << FletcherReevesCG(p1,A,0.8,0.2,0.00001)[2] << endl;


    plt::xlim(0,1);
    plt::ylim(0,1);
    plt::show();*/

    //TestAlgorithmOnRandomUniform(p1,3,12,0);
    /*vector<double> P = {0,0,0};

    for(double i = 4; i < 10; i += 0.2){
        for(double j = 4; j < 10; j += 0.2){
            for(double k = 4; k < 10; k += 0.2){
                P = {i,j,k};
                cout << "Point: " << i << ", " << j << ", " << k << endl;
                vector<double> start = findStartingPoint(p1,P);
                cout << "Start: " << start[0] << ", " << start[1] << ", " << start[2] << ", " << start[3] << ", " << start[4] << ", " << start[5] << endl;
                plotSurface(P,p1,start[0] - start[2], start[1] - start[3], start[0] + start[2], start[1] + start[3]);
                plt::show();
                vector<double> res;
                res = NewtonMethod(p1,P,start[0],start[1],0.00001,.5);
                cout << "Newton u: " << res[0] << " v: " << res[1] << " V: " << res[2] << " u0: " << res[3] << " v0: " << res[4] << " iterations: " << res[5] << endl;
                res = gradientDescent(p1,P,start[0],start[1],0.00001);
                cout << "GD u: " << res[0] << " v: " << res[1] << " V: " << res[2] << " iterations: " << res[3] << endl;
                res = NewtonMethodProjected(p1,P,start[0],start[1],0.00001,.5);
                cout << "NewtonProjected u: " << res[0] << " v: " << res[1] << " V: " << res[2] << " u0: " << res[3] << " v0: " << res[4] << " iterations: " << res[5] << endl;
                cout << "GridSearch: " << gridSearch(p1,P) << endl;
            }
        }
    }*/

//    cout << p1.fillTop(0.1,0.1)[0] << " " << p1.fillTop(0.1,0.1)[1] << " " <<p1.fillTop(0.1,0.1)[2] << " " <<endl;
//    cout << p1.fillTopDerU(0.1,0.1)[0]<<" " << p1.fillTopDerU(0.1,0.1)[1]<< " " <<p1.fillTopDerU(0.1,0.1)[2] <<" " << endl;
//    cout << p1.fillTopDerUU(0.1,0.1)[0]<< " " <<p1.fillTopDerUU(0.1,0.1)[1]<<" " << p1.fillTopDerUU(0.1,0.1)[2] << " " <<endl;
//    cout << p1.fillTopDerUV(0.1,0.1)[0]<<" " << p1.fillTopDerUV(0.1,0.1)[1]<< " " <<p1.fillTopDerUV(0.1,0.1)[2] <<" " << endl;



    /*vector<double> T, Tu, Tuu, Tuv, axis;
    for(double x = 0; x <= 1; x += 0.01){
        axis.push_back(x);
        T.push_back(p1.distanceToTopPoint(x,0.1,A));
        Tu.push_back(p1.completeDistanceTopDerV(x,0.1,A));
        Tuu.push_back(p1.completeDistanceTopDerUV2(x,0.1,A));
        Tuv.push_back(p1.completeDistanceTopDerVV2(x,0.1,A));
    }

    plt::plot(axis,T,"g");

    plt::plot(axis,Tu, "b");

    plt::plot(axis,Tuu, "r");
    plt::show();
    for(int i = 0; i < 3; i++) {
        vector<double> Tx, Txu, Txuu, Txuv, axisX;
        for (double x = 0; x <= 1; x += 0.01) {
            axisX.push_back(x);
            Tx.push_back(p1.fillTop(x, 0.1)[i]);
            Txu.push_back(p1.fillTopDerU(x, 0.1)[i]);
            Txuu.push_back(p1.fillTopDerUU(x, 0.1)[i]);
            Txuv.push_back(p1.fillTopDerUV(x, 0.1)[i]);
        }

        plt::plot(axis, Tx, "g");

        plt::plot(axis, Txu, "b");

        plt::plot(axis, Txuu, "r");
        plt::show();
    }*/
    //cout << "Distance via grid seach: " << min_dist << endl;
    //cout << "Distance via gradient descent: " << gradientDescent(p1,A,0.8,0.2,0.00001, 0.001) << endl;

    /*cout << "Distance via Newton's method: " << NewtonMethod(p1,A,0.5,0.5, 0.0001, .01) << endl;
    vector<vector<double>> xs,ys,zs;

    for (double y = 0; y < 20; y += 0.5) {
        vector<double> xr,yr,zr;
        for (double x = 5; x < 25; x += 0.5) {
            xr.push_back(x);
            yr.push_back(y);
            zr.push_back(f(x,y));
        }
        xs.push_back(xr);
        ys.push_back(yr);
        zs.push_back(zr);
    }

    vector<double> x2,y2,z2,u2,v2,w2;
    for (double y = 0; y < 20; y += 0.5) {
        for (double x = 5; x < 25; x += 0.5) {
            x2.push_back(x);
            y2.push_back(y);
            z2.push_back(-10);
            double tU, tV;
            tU = fdx(x,y);
            tV = fdy(x,y);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }


    plt::plot_surface_with_vector_field(xs,ys,zs,x2,y2,z2,u2,v2,w2);

    cout << NewtonMethodTest(18,14,0.000001,.1) << endl;
    plt::show();*/

//    cout << NewtonMethodSimpleFunc(1,0.00001);
}



