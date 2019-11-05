#include <iostream>
#include <vector>
#include <ctime>
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
    }while(abs(grad[0]) + abs(grad[1]) > eps && iter < 100);

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

void DrawDistWithDirection(Model &model, const vector<double> &P, vector<double> xs,vector<double> ys,vector<double> zs,vector<double> is,vector<double> js,vector<double> ks) {
    vector<vector<double>> u, v, dist;

    double min_dist = 100000;
    double max_dist = -100000;
    for (double i = 0; i < 1; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j < 1; j += 0.05) {
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
    for (double i = 0; i < 1; i += 0.05) {
        for (double j = 0; j < 1; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.squaredTopDistDerU(j, i, P);
            tV = model.squaredTopDistDerU(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, xs, ys, zs, is, js, ks, 0.1);

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
            DrawDistWithDirection(model,A,xs,ys,zs,is,js,ks);
            return {0,v + alpha, model.distanceToTopPoint(0,v + alpha,A), u0, v0, static_cast<double>(iter)};
        }
        if (u > 1) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,1,v,0,1,eps);
            cout << "u: " << u << " v: " << v << endl;
            DrawDistWithDirection(model,A,xs,ys,zs,is,js,ks);
            return {1,v + alpha, model.distanceToTopPoint(1,v + alpha,A), u0, v0, static_cast<double>(iter)};
        }
        if (v < 0) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,u,0,1,0,eps);
            cout << "u: " << u << " v: " << v << endl;
            DrawDistWithDirection(model,A,xs,ys,zs,is,js,ks);
            return {u + alpha,0, model.distanceToTopPoint(u + alpha,0,A), u0, v0, static_cast<double>(iter)};
        }
        if (v > 1) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,u,1,1,0,eps);
            cout << "u: " << u << " v: " << v << endl;
            DrawDistWithDirection(model,A,xs,ys,zs,is,js,ks);
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
    DrawDistWithDirection(model,A,xs,ys,zs,is,js,ks);

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

void TestAlgorithmOnRandomCloseToSurface(Model model, int algorithm, int preprocess, int seed, int plot = 0, int iterations = 1000, double GSAcc = 0.05, double eps = 0.0001){
    int reliability = 0;
    double avgError = 0;
    int sumIters = 0;
    srand(seed);
    for (int i = 0; i < iterations; i++) {
        vector<double> P = model.fillTop(double(rand())/RAND_MAX,double(rand())/RAND_MAX);

        P[0] += (double(rand())/RAND_MAX)*0.01;
        P[1] += (double(rand())/RAND_MAX)*0.01;
        P[2] += (double(rand())/RAND_MAX)*0.01;

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
                algV = NewtonMethodProjected(model,P,start[0],start[1],eps,0.01,0);
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
        default:
            break;
    }
    cout << "reliability: " << double(reliability)/iterations << " average error: " << avgError/iterations << " average iterations: " << double(sumIters)/iterations << endl;
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
    vector<vector<double>> x, y, z;


    drawPart(p1, x, y, z);

    TestAlgorithmOnRandomCloseToSurface(p1,1,12,0,1);

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



