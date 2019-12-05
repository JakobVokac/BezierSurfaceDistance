#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <float.h>
#include "model.h"
#include "vectorMath.h"
#include "plotting.h"

int onedimcalls = 0;

double inrange(double x, double low, double high){
    return max(min(x, high),low);
}

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

//double constraintFunc(double u, double v){
//    return pow((u < 0 ? -u : 0),3.0)
//           + pow((u > 1 ? u-1 : 0),3.0)
//           + pow((v < 0 ? -v : 0),3.0)
//           + pow((v > 1 ? v-1 : 0),3.0);
//}
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
       /* if( u < 0 ){
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
        }*/
        dist_old = dist;
        dist = model.distanceToTopPoint(u,v,A) + lambda*constraintFunc(u,v);
        us.push_back(u);
        vs.push_back(v);
        HModifier /= 2;
    }while(abs(dist_old - dist) > eps && iter < 20);
//    plt::plot(us,vs);
//    plt::show();
//    cout << "Newton's method iterations: " << iter << " u: " << u << " v: " << v << endl;
    return {u,v,dist,u0,v0,static_cast<double>(iter)};
}

double NewtonForLinearCG(Model model, vector<double> A, double u, double v, double dU, double dV, double eps){
    double V = model.distanceToTopPoint(u,v,A);
    double V_old;
    double t = 0;
    int iter = 0;

    do{
        t -= (dU * model.squaredTopDistDerU(u + dU*t,v + dV*t,A) + dV * model.squaredTopDistDerV(u + dU*t,v + dV*t,A))
                /(dU*dU*model.squaredTopDistDerUU(u + dU*t,v + dV*t,A)
                + 2*dU*dV*model.squaredTopDistDerUV(u + dU*t,v + dV*t,A)
                + dV*dV*model.squaredTopDistDerVV(u + dU*t,v + dV*t,A));

//        cout << "t: " << t << endl;
        V_old = V;
        V = model.distanceToTopPoint(u + dU*t, v + dV*t, A);
        iter++;
    }while(abs(V-V_old) > eps && 0 <= u + dU*t && u + dU*t <= 1 && v + dV*t <= 1 && 0 <= v + dV*t && iter < 10);

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



vector<double> FletcherReevesCG(Model model, vector<double> A, double u, double v, double eps=0.0001, int plot = 0){
    double dU, dV, alpha, beta, dU2, dV2, u0 = u, v0 = v;
    int iter = 0;
    //cout << "u: " << u << " v: " << v << endl;

    // p 0 <- - Grad f 0 = Grad f(u,v)
    dU2 = -model.squaredTopDistDerU(u,v,A);
    dV2 = -model.squaredTopDistDerV(u,v,A);
//    vector<double> xs, ys, is, js,zs,ks;
//    xs.push_back(u);
//    ys.push_back(v);
//    zs.push_back(0);
//    is.push_back(dU2);
//    js.push_back(dV2);
//    ks.push_back(0);
    do{
        dU = dU2;
        dV = dV2;

        // Compute alpha k and set xk1 = xk + alpha k * p k
        alpha = NewtonForLinearCG(model,A,u,v,dU,dV,eps);
//        cout << "alpha: " << alpha << endl;
//        cout << "dU: " << dU << " dV: " << dV << endl;
//        if(plot){
//            cout << "u: " << u << " v: " << v << " dU: " << dU << " dV: " << dV << " alpha: " << alpha << endl;
//        }
        u += alpha * dU;
        v += alpha * dV;

        if (u < 0) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,0,v,0,1,eps);
//            cout << "u: " << u << " v: " << v << endl;
//            DrawDistWithDirection(model,A);
            return {0,v + alpha, model.distanceToTopPoint(0,v + alpha,A), u0, v0, static_cast<double>(iter)};
        }
        if (u > 1) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,1,v,0,1,eps);
//            cout << "u: " << u << " v: " << v << endl;
//            DrawDistWithDirection(model,A);
            return {1,v + alpha, model.distanceToTopPoint(1,v + alpha,A), u0, v0, static_cast<double>(iter)};
        }
        if (v < 0) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,u,0,1,0,eps);
//            cout << "u: " << u << " v: " << v << endl;
//            DrawDistWithDirection(model,A);
            return {u + alpha,0, model.distanceToTopPoint(u + alpha,0,A), u0, v0, static_cast<double>(iter)};
        }
        if (v > 1) {
//            plt::plot(us,vs);
//            plt::show();
            alpha = NewtonForLinearCG(model,A,u,1,1,0,eps);
//            cout << "u: " << u << " v: " << v << endl;
//            DrawDistWithDirection(model,A);
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
//        xs.push_back(u);
//        ys.push_back(v);
//        zs.push_back(50);
//        is.push_back(dU2);
//        js.push_back(dV2);
//        ks.push_back(0);
//        cout << "u: " << u << " v: " << v << endl;
    }while(abs(dU + dV) > eps && iter < 10); // while Grad
//    cout << "u: " << u << " v: " << v << endl;
//    DrawDistWithDirection(model,A);
//
//    if(plot){
//        plt::plot(xs,ys, "g");
//        plt::show();
//    }
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

vector<double> orthogonalProjection1D(CBezier b, vector<double> P, double t, double eps, int iter) {
    double dist = 100000, distOld;
    double dt = 0;
    int i = 0;
    vector<double> N = cross(sub(b.C2,b.C0),sub(b.C1,b.C0));
    double D = abs(dot(N,P)-dot(N,b.C0))/magnitude(N);
    vector<double> Pnew = sub(P,mul(N,D));

    vector<double> p0,p0P,n,q,c1,c2,cq,c,s,u;
    double k;

    p0 = cubicBezier(t,b);

    p0P = sub(Pnew, p0);
    dist = magnitude(p0P);
    do {

        c1 = cubicBezierDer1(t,b);
        c2 = cubicBezierDer2(t,b);

        k = curvature3D(c1,c2);
        c = circleCenterDir3D(c1,c2);
        s = add(p0,div(c,k));
        u = sub(Pnew,s);
        q = add(s, div(u,k*magnitude(u)));
        cq = sub(q,p0);
        dt = sign3D(c1, cq) * sqrt(abs( 2*magnitude(cross(c1, cq))/magnitude(cross(c1, c2))));

        t += dt;
        double t2 = inrange(t,0,1);
        if( t != t2 ){
            return {t2,magnitude(sub(P,cubicBezier(t2,b))), static_cast<double>(i)};
        }

        p0 = cubicBezier(t,b);

        p0P = sub(Pnew, p0);
        distOld = dist;
        dist = magnitude(p0P);

        i++;

    }while(abs(dist - distOld) > eps && i < iter);

    return {t,magnitude(sub(P,cubicBezier(t,b))), static_cast<double>(i)};
}
CBezier raiseDeg(QBezier q){
    return {q.C0,add(mul(q.C0,1.0/3.0),mul(q.C1,2.0/3.0)),add(mul(q.C1,2.0/3.0),mul(q.C2,1.0/3.0)),q.C2};
}

double bisectionBezier(const CBezier &c, const vector<double> &P, double t1, double t2){
    double tm, hm, ht;
    for (int i = 0; i < 8; i++) {
        tm = (t1+t2)/2;
        hm = sqDistBezierDer1(c,tm,P);
        ht = sqDistBezierDer1(c,t1,P);
        if(hm == 0)
            return tm;
        if((hm < 0 && ht < 0) || (hm > 0 && ht > 0)){
            t1 = tm;
        }else{
            t2 = tm;
        }
    }
    return tm;
}

double QuadInterBezier(const CBezier &c, const vector<double> &P, double t1, double t2, double t3){
    double h1, h2, h3, hn;
    double tn;
    for (int i = 0; i < 12; i++) {
        h1 = sqDistBezier(c,t1,P);
        h2 = sqDistBezier(c,t2,P);
        h3 = sqDistBezier(c,t3,P);

        if(t1 == t2 || t2 == t3 || t3 == t1){
            break;
        }
        tn = (t1*t1*(h2-h3) + t2*t2*(h3-h1) + t3*t3*(h1-h2))/(2.0*(t1*(h2-h3) + t2*(h3-h1) + t3*(h1-h2)));

        if (tn < 0) {
            tn = 0;
        } else if (tn > 1) {
            tn = 1;
        }


        hn = sqDistBezier(c,tn,P);


        double hmax = h1;
        if(h2 > hmax)
            hmax = h2;
        if(h3 > hmax)
            hmax = h3;
        if(hn > hmax)
            hmax = hn;
        if(hmax == h1)
            t1 = tn;
        if(hmax == h2)
            t2 = tn;
        if(hmax == h3)
            t3 = tn;
        if(hmax == hn)
            tn = tn;
    }
    return inrange(tn,0,1);
}

vector<double> QuadInterToNewtonEdge(CBezier c, vector<double> P, double t, double eps, int iter, int quadIter = 3, int plot = 0){
    double dist = sqDistBezier(c,t,P);
    double dist_old;
    double t1 = 0,t2,t3 = 1;

    if(plot)
        cout << "t before quad: " << t << endl;
    t2 = t;
    if(t < 0){
        if(plot)
            cout << "final t: " << t << endl;

        return {0,sqrt(sqDistBezier(c,0,P)),static_cast<double>(iter)};
    }else if(t > 1){
        if(plot)
            cout << "final t: " << t << endl;

        return {1,sqrt(sqDistBezier(c,1,P)),static_cast<double>(iter)};
    }

    double temp1 = QuadInterBezier(c, P, t1, t2, t3);
    double temp2 = bisectionBezier(c, P, t1, t3);
    if(sqDistBezier(c,temp1,P) < sqDistBezier(c,temp2,P))
        t = temp1;
    else
        t = temp2;
    if(t <= 0){
        if(plot)
            cout << "final t: " << t << endl;

        return {0,sqrt(sqDistBezier(c,0,P)),static_cast<double>(iter)};
    }else if(t >= 1){
        if(plot)
            cout << "final t: " << t << endl;

        return {1,sqrt(sqDistBezier(c,1,P)),static_cast<double>(iter)};
    }
    return {t,sqrt(sqDistBezier(c,t,P)),static_cast<double>(iter)};
    if(plot)
        cout << "t after quad: " << t << endl;
    do{
        iter++;
        t +=  sqDistBezierDer1(c,t,P)/sqDistBezierDer2(c,t,P);
        dist_old = dist;
        dist =  sqDistBezier(c,t,P);
        if(t < 0){
            if(plot)
                cout << "final t: " << t << endl;

            return {0,sqrt(sqDistBezier(c,0,P)),static_cast<double>(iter)};
        }else if(t > 1){
            if(plot)
                cout << "final t: " << t << endl;

            return {1,sqrt(sqDistBezier(c,1,P)),static_cast<double>(iter)};
        }
        if(plot)
            cout << "t: " << t << endl;
    }while(abs(dist - dist_old) > eps && iter < 100);
    if(plot)
        cout << "final t: " << t << endl;

    return {t,sqrt(sqDistBezier(c,t,P)),static_cast<double>(iter)};
};

vector<double> Geometric1D(CBezier c,vector<double> P,double v,double eps,int iter){
    return QuadInterToNewtonEdge(c,P,v,eps,iter);
//    return orthogonalProjection1D(c,P,v,eps,iter);
}

vector<double> orthogonalProjection(Model m, vector<double> P, double u, double v, double eps = 0.0000001, int iter = 15, int plot = 0){

    double u0 = u, v0 = v;
    double dist = 100000, distOld;
    double lam1 = 1, lam2 = 2;
    double dt = 1;
    int i = 0;
    int rangecount = 1;
    vector<double> p0,p0P,du,duu,duv,dvv,dv,n,coefs,Ps,q,c1,c2,cq,s;
    double g11, g12, g22, h11, h12, h22, g21, h21, k;

    p0 = m.fillTop(u, v);

    p0P = sub(P, p0);
    distOld = dist;
    dist = magnitude(p0P);
    do {

        du = m.fillTopDerU(u, v);
        dv = m.fillTopDerV(u, v);
        duu = m.fillTopDerUU(u, v);
        duv = m.fillTopDerVV(u, v);
        dvv = m.fillTopDerUV(u, v);

        n = cross(du, dv);
        n = div(n, magnitude(n));

        coefs = solve3Dlinear({du[0], dv[0], n[0], p0P[0]}, {du[1], dv[1], n[1], p0P[1]}, {du[2], dv[2], n[2], p0P[2]});
        lam1 = coefs[0];
        lam2 = coefs[1];

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


        s = add(p0, mul(n, 1.0 / k));

        Ps = sub(P, s);
        q = add(s, div(Ps, magnitude(Ps) * abs(k)));

        c1 = add(mul(du, lam1), mul(dv, lam2));
        c2 = add(add(mul(duu, lam1 * lam1), mul(duv, 2 * lam1 * lam2)), mul(dvv, lam2 * lam2));


        cq = sub(q, p0);
        dt = sign3D(c1, cq) * sqrt(abs(2 * magnitude(cross(c1, cq)) / magnitude(cross(c1, c2))));

        u += dt * lam1;
        v += dt * lam2;
        i++;

        if(u < 0 && rangecount > 2){
            onedimcalls++;
            u = 0;
            if(v < 0){
                v = 0;
                vector<double> res1, res2;


                res1 = Geometric1D(m.symCurveTop,P,1-v,eps,iter-i);
                res2 = Geometric1D(m.bendCurve,P,u,eps,iter-i);

                if(plot) {
                    plotEdgeDistSq(m.symCurveTop, P);
                    plotEdgeDistSq(m.bendCurve, P);
                }
                if(res1[1] < res2[1]){
                    v = 1-res1[0];
                }else{
                    u = res1[0];
                }
                double distFinal = res1[1];
                i += res1[2];
                return {u,v,distFinal,u0,v0,static_cast<double>(i)};
            }
            if(v > 1){
                v = 1;
                vector<double> res1, res2;

                res1 = Geometric1D(m.symCurveTop,P,1-v,eps,iter-i);
                res2 = Geometric1D(m.leafCurve,P,u,eps,iter-i);
                if(plot) {
                    plotEdgeDistSq(m.symCurveTop, P);
                    plotEdgeDistSq(m.leafCurve, P);
                }
                if(res1[1] < res2[1]){
                    v = 1-res1[0];
                }else{
                    u = res1[0];
                }
                double distFinal = res1[1];
                i += res1[2];
                return {u,v,distFinal,u0,v0,static_cast<double>(i)};
            }

            vector<double> res = QuadInterToNewtonEdge(m.symCurveTop,P,1-v,eps,iter-i);
            if(plot)
            plotEdgeDistSq(m.symCurveTop, P);
            v = res[0];
            double distFinal = res[1];
            i += res[2];
            return {u,v,distFinal,u0,v0,static_cast<double>(i)};
        } else
            rangecount++;
        if(u > 1 && rangecount > 2){
            onedimcalls++;

            u = 1;
            if(v < 0){
                v = 0;
                vector<double> res1, res2;

                res1 = Geometric1D(raiseDeg(m.sinCurveTop),P,1-v,eps,iter-i);
                res2 = Geometric1D(m.bendCurve,P,u,eps,iter-i);
                if(plot) {
                    plotEdgeDistSq(raiseDeg(m.sinCurveTop), P);
                    plotEdgeDistSq(m.bendCurve, P);
                }
                if(res1[1] < res2[1]){
                    v = 1-res1[0];
                }else{
                    u = res1[0];
                }
                double distFinal = res1[1];
                i += res1[2];
                return {u,v,distFinal,u0,v0,static_cast<double>(i)};
            }
            if(v > 1){
                v = 1;
                vector<double> res1, res2;

                res1 = Geometric1D(raiseDeg(m.sinCurveTop),P,1-v,eps,iter-i);
                res2 = Geometric1D(m.leafCurve,P,u,eps,iter-i);
                if(plot) {
                    plotEdgeDistSq(raiseDeg(m.sinCurveTop), P);
                    plotEdgeDistSq(m.leafCurve, P);
                }
                if(res1[1] < res2[1]){
                    v = 1-res1[0];
                }else{
                    u = res1[0];
                }
                double distFinal = res1[1];
                i += res1[2];
                return {u,v,distFinal,u0,v0,static_cast<double>(i)};
            }
            vector<double> res = Geometric1D(raiseDeg(m.sinCurveTop),P,1-v,eps,iter-i);
            if(plot)
                plotEdgeDistSq(raiseDeg(m.sinCurveTop), P);

            v = 1-res[0];
            double distFinal = res[1];
            i += res[2];
            return {u,v,distFinal,u0,v0,static_cast<double>(i)};
        } else
            rangecount++;
        if(v < 0 && rangecount > 2){
            onedimcalls++;

            v = 0;
            vector<double> res = Geometric1D(m.bendCurve,P,u,eps,iter-i);
            if(plot)
                plotEdgeDistSq(m.bendCurve, P);

            u = res[0];
            double distFinal = res[1];
            i += res[2];
            return {u,v,distFinal,u0,v0,static_cast<double>(i)};
        } else
            rangecount++;
        if(v > 1 && rangecount > 2){
            onedimcalls++;

            v = 1;
            vector<double> res = Geometric1D(m.leafCurve,P,u,eps,iter-i);
            if(plot)
                plotEdgeDistSq(m.leafCurve, P);
            u = res[0];
            double distFinal = res[1];
            i += res[2];
            return {u,v,distFinal,u0,v0,static_cast<double>(i)};
        } else
            rangecount++;

        p0 = m.fillTop(u, v);

        p0P = sub(P, p0);
        distOld = dist;
        dist = magnitude(p0P);

    }while((abs(dist - distOld) > eps && i < iter) || u != inrange(u,0,1) || v != inrange(v,0,1) );

    return {u,v,magnitude(sub(P,m.fillTop(u,v))),u0,v0,static_cast<double>(i)};
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
void TestAlgorithmOnRandomFarFromSurface(Model model, int algorithm, int preprocess, int seed, int plot = 0, int iterations = 1000, double GSAcc = 0.05, double eps = 0.0000001){
    cout << "Far from surface test: Preprocessing: " << preprocess << " ";
    int reliability = 0;
    double avgError = 0;
    int sumIters = 0;
    srand(seed);
    int count1D = 0;
    time_t avgTime = 0;

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration1;
    chrono::microseconds t_duration2;


    for (int i = 0; i < iterations; i++) {
        vector<double> P = model.fillTop(double(rand())/RAND_MAX,double(rand())/RAND_MAX);

        P[0] += (double(rand())/RAND_MAX)*20;
        P[1] += (double(rand())/RAND_MAX)*20;
        P[2] += (double(rand())/RAND_MAX)*20;

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

        int onedimold = onedimcalls;
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
                algV = NewtonMethod(model,P, start[0], start[1], eps,0.1,100,0);
                if(plot)
                    cout << "Newton's Method ";
                break;
            case 3:
                t_start = chrono::high_resolution_clock::now();
                algV = NewtonMethodProjected(model,P,start[0],start[1],eps,0.01,0);
                t_stop = chrono::high_resolution_clock::now();
                t_duration2 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
                if(preprocess) {
                    avgTime += t_duration1.count() + t_duration2.count();
//                    cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
//                    cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
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
                    //cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
                    //cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
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
            if(onedimold != onedimcalls){
                count1D++;
            }
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
   // cout << "% onedimcalls: " << double(onedimcalls)/iterations << " error by 1D search %: " << double(count1D)/onedimcalls << endl;
}
void TestAlgorithmOnRandomCloseToSurface(Model model, int algorithm, int preprocess, int seed, int plot = 0, int iterations = 1000, double GSAcc = 0.05, double eps = 0.0000001){
    cout << "Close to surface test: Preprocessing: " << preprocess << " ";
    int reliability = 0;
    double avgError = 0;
    int sumIters = 0;
    srand(seed);
    int count1D = 0;
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

        int onedimold = onedimcalls;
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
                algV = NewtonMethod(model,P, start[0], start[1], eps,0.1,100,0);
                if(plot)
                    cout << "Newton's Method ";
                break;
            case 3:
                t_start = chrono::high_resolution_clock::now();
                algV = NewtonMethodProjected(model,P,start[0],start[1],eps,0.01,0);
                t_stop = chrono::high_resolution_clock::now();
                t_duration2 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
                if(preprocess) {
                    avgTime += t_duration1.count() + t_duration2.count();
//                    cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
//                    cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
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
                    //cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
                    //cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
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
            if(onedimold != onedimcalls){
                count1D++;
            }
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
   // cout << "% onedimcalls: " << double(onedimcalls)/iterations << " error by 1D search %: " << double(count1D)/onedimcalls << endl;
}

void TestAlgorithmPreciseCloseToSurface(Model model, int algorithm, int preprocess, int seed, int plot = 0, int iterations = 1000, double GSAcc = 0.05, double eps = 0.0000001){
    cout << "Close to surface test: Preprocessing: " << preprocess << " ";
    int reliability = 0;
    double avgError = 0;
    int sumIters = 0;
    srand(seed);
    int count1D = 0;
    time_t avgTime = 0;

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration1;
    chrono::microseconds t_duration2;


    for (int i = 0; i < iterations; i++) {
        double u = double(rand())/RAND_MAX, v = double(rand())/RAND_MAX;
        vector<double> P = model.fillTop(u,v);
        vector<double> Pu = model.fillTopDerU(u,v);
        vector<double> Pv = model.fillTopDerV(u,v);
        vector<double> N = cross(Pu,Pv);
        N = div(N,magnitude(N));
        double Pdist = (double(rand())/RAND_MAX)*0.5;

        P = add(P,mul(N,Pdist));

        vector<double> start;
        vector<double> algV;
        if(algorithm == 5)
            preprocess = 0;

        if(preprocess){
            t_start = chrono::high_resolution_clock::now();
            start = findStartingPoint(model, P, preprocess);
            t_stop = chrono::high_resolution_clock::now();
            t_duration1 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
        }else
            start = {0.5,0.5};

        int onedimold = onedimcalls;
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
                algV = NewtonMethod(model,P, start[0], start[1], eps,0.1,100,0);
                if(plot)
                    cout << "Newton's Method ";
                break;
            case 3:
                t_start = chrono::high_resolution_clock::now();
                algV = NewtonMethodProjected(model,P,start[0],start[1],eps,0.01,0);
                t_stop = chrono::high_resolution_clock::now();
                t_duration2 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
                if(preprocess) {
                    avgTime += t_duration1.count() + t_duration2.count();
//                    cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
//                    cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
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
                    //cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
                    //cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
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

        if(algV[2] <= Pdist + eps) {
            reliability++;
        }else{
            avgError += (algV[2] - Pdist)/Pdist;
            if(onedimold != onedimcalls){
                count1D++;
            }
            if(plot) {
                cout << "Grid Search: " << Pdist << endl;
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
    // cout << "% onedimcalls: " << double(onedimcalls)/iterations << " error by 1D search %: " << double(count1D)/onedimcalls << endl;
}

void TestAlgorithmPreciseFarFromSurface(Model model, int algorithm, int preprocess, int seed, int plot = 0, int iterations = 1000, double GSAcc = 0.05, double eps = 0.0000001){
    cout << "Close to surface test: Preprocessing: " << preprocess << " ";
    int reliability = 0;
    double avgError = 0;
    int sumIters = 0;
    srand(seed);
    int count1D = 0;
    time_t avgTime = 0;

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration1;
    chrono::microseconds t_duration2;


    for (int i = 0; i < iterations; i++) {

        double Pdist = (double(rand())/RAND_MAX)*20.0 + 3.0;

        CBezier c;
        vector<double> dir;
        double t, u, v;
        vector<double> P, Pu, Pv, N;
        switch (i/(iterations/5)){
            case 0:
                c = model.leafCurve;
                dir = {0,0,1};
                t = double(rand())/RAND_MAX;
                P = add(cubicBezier(t,c), mul(dir,Pdist));
                break;
            case 1:
                c = model.symCurveTop;
                dir = {0,-1,0};
                t = double(rand())/RAND_MAX;
                P = add(cubicBezier(t,c), mul(dir,Pdist));
                break;
            case 2:
                c = model.bendCurve;
                dir = {0,0,-1};
                t = double(rand())/RAND_MAX;
                P = add(cubicBezier(t,c), mul(dir,Pdist));
                break;
            case 3:
                c = raiseDeg(model.sinCurveTop);
                dir = {0,1,0};
                t = double(rand())/RAND_MAX;
                P = add(cubicBezier(t,c), mul(dir,Pdist));
                break;
            case 4:
                u = double(rand())/RAND_MAX, v = double(rand())/RAND_MAX;
                P = model.fillTop(u,v);
                Pu = model.fillTopDerU(u,v);
                Pv = model.fillTopDerV(u,v);
                N = cross(Pu,Pv);
                N = div(N,magnitude(N));
                P = add(P,mul(N,Pdist));
                break;
            default:
                break;

        }

        vector<double> start;
        vector<double> algV;
        if(algorithm == 5)
            preprocess = 0;

        if(preprocess){
            t_start = chrono::high_resolution_clock::now();
            start = findStartingPoint(model, P, preprocess);
            t_stop = chrono::high_resolution_clock::now();
            t_duration1 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
        }else
            start = {0.5,0.5};

        int onedimold = onedimcalls;
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
                algV = NewtonMethod(model,P, start[0], start[1], eps,0.1,100,0);
                if(plot)
                    cout << "Newton's Method ";
                break;
            case 3:
                t_start = chrono::high_resolution_clock::now();
                algV = NewtonMethodProjected(model,P,start[0],start[1],eps,0.01,0);
                t_stop = chrono::high_resolution_clock::now();
                t_duration2 = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);
                if(preprocess) {
                    avgTime += t_duration1.count() + t_duration2.count();
//                    cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
//                    cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
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
                    //cout << "Duration: " << t_duration1.count() + t_duration2.count() << " microseconds " << endl;
                }else {
                    avgTime += t_duration2.count();
                    //cout << "Duration: " << t_duration2.count() << " microseconds " << endl;
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

        if(algV[2] <= Pdist + eps) {
            reliability++;
        }else{
            avgError += (algV[2] - Pdist)/Pdist;
            if(onedimold != onedimcalls){
                count1D++;
            }
            if(plot) {
                cout << "Grid Search: " << Pdist << endl;
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
    // cout << "% onedimcalls: " << double(onedimcalls)/iterations << " error by 1D search %: " << double(count1D)/onedimcalls << endl;
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

    vector<vector<double>> xs,ys,zs;
    drawPart(p1,xs,ys,zs,1,0);
    plt::plot_surface(xs,ys,zs);
    plt::show();
//    TestAlgorithmOnRandomFarFromSurface(p1,6,12,3,0);
//    TestAlgorithmOnRandomCloseToSurface(p1,6,12,3,0);
    TestAlgorithmPreciseCloseToSurface(p1,6,8,3,0);
    TestAlgorithmPreciseFarFromSurface(p1,6,8,3,0);


//    vector<double> A = {4,4,4};
//    vector<vector<double>> u,v,dist;
//    double min_dist = 100000;
//    for (double i = -0; i <= 1; i += 0.05) {
//        vector<double> uRow, vRow, distRow;
//        for (double j = -0; j <= 1; j += 0.05) {
//            uRow.push_back(j);
//            vRow.push_back(i);
//            distRow.push_back(p1.distanceToTopPoint(j,i,A) );
//            if(min_dist > p1.distanceToTopPoint(j,i,A))
//                min_dist = p1.distanceToTopPoint(j,i,A);
//        }
//        u.push_back(uRow);
//        v.push_back(vRow);
//        dist.push_back(distRow);
//    }
//    double scale = 0.10;
//    vector<double> x2,y2,z2,u2,v2,w2;
//    for (double i = -0; i <= 1; i += 0.05) {
//        for (double j = -0; j <= 1; j += 0.05) {
//            x2.push_back(j);
//            y2.push_back(i);
//            z2.push_back(min_dist);
//            double tU, tV;
//            tU = p1.completeDistanceTopDerU(j,i,A);
//            tV = p1.completeDistanceTopDerV(j,i,A);
//            u2.push_back(tU);
//            v2.push_back(tV);
//            w2.push_back(0);
//        }
//    }
//    cout << "Rough expected distance: " << min_dist << endl;
//    vector<vector<double>> CGus, CGvs, CGVs;
//    vector<vector<double>> Nus, Nvs, NVs;
//
//
//    plt::plot_surface_with_vector_field(u,v,dist,x2,y2,z2,u2,v2,w2,0.03);
//    plt::show();

//    for (int i = 8; i <= 12; i+= 4) {
////        TestAlgorithmOnRandomCloseToSurface(p1,1,i,0,0);
////        TestAlgorithmOnRandomFarFromSurface(p1,1,i,0,0);
////        TestAlgorithmOnRandomCloseToSurface(p1,2,i,0,0);
////        TestAlgorithmOnRandomFarFromSurface(p1,2,i,0,0);
//        TestAlgorithmOnRandomCloseToSurface(p1,3,i,0,0);
//        TestAlgorithmOnRandomFarFromSurface(p1,3,i,0,0);
//
//        TestAlgorithmOnRandomCloseToSurface(p1,6,i,0,0);
//        onedimcalls = 0;
//        TestAlgorithmOnRandomFarFromSurface(p1,6,i,0,0);
//        onedimcalls = 0;
//    }
}



