//
// Created by s3179222 on 12/3/19.
//
#include "GeometricProjection.h"

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

vector<double> geometricProjectionTest(){

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

        double k = curvature2D(F1, F2);
        vector<double> c = circleCenterDir2D(F1, F2);
        vector<double> s = add(F,div(c,k));
        vector<double> u = sub(P,s);
        Q = add(s, div(u,k*magnitude(u)));


        drawC(Q0,Q1,Q2,Q3);
        drawLine(F,P);
        drawLine(F,s);
        drawLine(F,add(F,F2));
        plt::show();

        dt = sign2D(F1, sub(Q, F)) * sqrt(abs(2 * det(F1, sub(Q, F)) / det(F1, F2)));

        t += dt;
        i++;
        cout << "dt: " << dt << " iter: " << i << endl;
    }

    drawC(Q0,Q1,Q2,Q3);
    drawLine(F,P);
    plt::show();
    return {B(t,Q0,Q1,Q2,0), B(t,Q0,Q1,Q2,1)};
}