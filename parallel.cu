//
// Created by s3179222 on 11/22/19.
//

#include <iostream>
#include <chrono>
#include "model.h"
#include <cmath>

struct QBezierGPU{
    double C0[3], C1[3], C2[3];
};

struct CBezierGPU{
    double C0[3], C1[3], C2[3], C3[3];
};

struct gpumodel{
    CBezierGPU leafCurve, bendCurve, symCurveTop;
    QBezierGPU sinCurveTop;
    double Q_b_sym[3], Q_l_sym[3], Q_b_sin[3], Q_l_sin[3];
};

__device__
double cubicBezier(double t, int i, CBezierGPU *b){

    return pow((1 - t), 3.0) * b->C0[i] + 3 * pow((1 - t), 2.0) * t * b->C1[i] +
           3 * (1 - t) * pow(t, 2.0) * b->C2[i] + pow(t, 3.0) * b->C3[i];

}

__device__
double cubicBezierDer1(double t, int i, CBezierGPU *b) {

    return 3 * pow((1 - t), 2.0) * (b->C1[i] - b->C0[i]) + 6 * t * (1 - t) * (b->C2[i] - b->C1[i]) +
           3 * pow(t, 2.0) * (b->C3[i] - b->C2[i]);

}

__device__
double cubicBezierDer2(double t, int i, CBezierGPU *b) {

    return 6 * (1 - t) * (b->C2[i] - 2 * b->C1[i] + b->C0[i]) + 6 * t * (b->C3[i] - 2 * b->C2[i] + b->C1[i]);

}

__device__
double quadBezier(double t, int i, QBezierGPU *b){

    return pow((1 - t), 2.0) * b->C0[i] + 2 * (1 - t) * t * b->C1[i] +
           pow(t, 2.0) * b->C2[i];

}

__device__
double quadBezierDer1(double t, int i, QBezierGPU *b) {

    return 2 * (1 - t) * (b->C1[i] - b->C0[i]) + 2 * t * (b->C2[i] - b->C1[i]);

}

__device__
double quadBezierDer2(double t, int i, QBezierGPU *b) {

    return 2*b->C0[i] - 4*b->C1[i] + 2*b->C2[i];

}

__device__
void fillTop(double *point, double u, double v, gpumodel *m){
    for (int i = 0; i < 3; ++i) {
        point[i] = (1 - v) * cubicBezier(u, i, &m->bendCurve) + v * cubicBezier(u, i, &m->leafCurve) +
                   (1 - u) * cubicBezier(1 - v, i, &m->symCurveTop) + u * quadBezier(1 - v, i, &m->sinCurveTop) -
                   ((1 - u) * (1 - v) * (m->Q_b_sym)[i] + u * v * (m->Q_l_sin)[i] +
                    u * (1 - v) * (m->Q_b_sin)[i] + (1 - u) * v * (m->Q_l_sym)[i]);
    }
}

__device__
void fillTopDerU(double *point, double u, double v, gpumodel *m) {
    for (int i = 0; i < 3; ++i) {
        point[i] = (1 - v) * cubicBezierDer1(u, i, &m->bendCurve) + v * cubicBezierDer1(u, i, &m->leafCurve)
                   - cubicBezier(1 - v, i, &m->symCurveTop) + quadBezier(1 - v, i, &m->sinCurveTop) -
                   (- (1 - v) * (m->Q_b_sym)[i] + v * (m->Q_l_sin)[i] +
                    (1 - v) * (m->Q_b_sin)[i] - v * (m->Q_l_sym)[i]);
    }
}

__device__
void fillTopDerUU(double *point, double u, double v, gpumodel *m){
    for (int i = 0; i < 3; ++i) {
        point[i] = (1 - v) * cubicBezierDer2(u, i, &m->bendCurve) + v * cubicBezierDer2(u, i, &m->leafCurve);
    }
}

__device__
void fillTopDerUV(double *point, double u, double v, gpumodel *m) {
    for (int i = 0; i < 3; ++i) {
        point[i] = - cubicBezierDer1(u, i, &m->bendCurve) + cubicBezierDer1(u, i, &m->leafCurve)
                   + cubicBezierDer1(1 - v, i, &m->symCurveTop) - quadBezierDer1(1 - v, i, &m->sinCurveTop) -
                   ( + m->Q_b_sym[i] + m->Q_l_sin[i] - m->Q_b_sin[i] - m->Q_l_sym[i]);
    }
}

__device__
void fillTopDerV(double *point, double u, double v, gpumodel *m){
    for (int i = 0; i < 3; ++i) {
        point[i] = - cubicBezier(u, i, &m->bendCurve) + cubicBezier(u, i, &m->leafCurve) -
                   (1 - u) * cubicBezierDer1(1 - v, i, &m->symCurveTop) - u * quadBezierDer1(1 - v, i, &m->sinCurveTop) -
                   (- (1 - u) * m->Q_b_sym[i] + u * m->Q_l_sin[i]
                    - u * m->Q_b_sin[i] + (1 - u) * m->Q_l_sym[i]);
    }
}

__device__
void fillTopDerVV(double *point, double u, double v, gpumodel *m){
    for (int i = 0; i < 3; ++i) {
        point[i] = (1 - u) * cubicBezierDer2(1 - v, i, &m->symCurveTop) + u * quadBezierDer2(1 - v, i, &m->sinCurveTop);
    }
}


void translateBezier(CBezier b, CBezierGPU *bgpu){
    for (int i = 0; i < 3; ++i) {
        bgpu->C0[i] = b.C0[i];
        bgpu->C1[i] = b.C1[i];
        bgpu->C2[i] = b.C2[i];
        bgpu->C3[i] = b.C3[i];
    }
}

void translateBezier(QBezier b, QBezierGPU *bgpu){
    for (int i = 0; i < 3; ++i) {
        bgpu->C0[i] = b.C0[i];
        bgpu->C1[i] = b.C1[i];
        bgpu->C2[i] = b.C2[i];
    }
}

void translatePoint(vector<double> q, double *qgpu){
    for (int i = 0; i < 3; ++i) {
        qgpu[i] = q[i];
    }
}

void translateModel(Model m, gpumodel *mgpu){

    translateBezier(m.symCurveTop,&mgpu->symCurveTop);
    translateBezier(m.leafCurve,&mgpu->leafCurve);
    translateBezier(m.bendCurve,&mgpu->bendCurve);
    translateBezier(m.sinCurveTop,&mgpu->sinCurveTop);
    translatePoint(m.Q_l_sin,mgpu->Q_l_sin);
    translatePoint(m.Q_b_sin,mgpu->Q_b_sin);
    translatePoint(m.Q_l_sym,mgpu->Q_l_sym);
    translatePoint(m.Q_b_sym,mgpu->Q_b_sym);

}


// -------------------------------------------- Math functions ------------------------------------------------------ //

__device__
void add(double *P, double *Q, double *R){
    R[0] = P[0] + Q[0];
    R[1] = P[1] + Q[1];
    R[2] = P[2] + Q[2];
}

__device__
void sub(double *P, double *Q, double *R){
    R[0] = P[0] - Q[0];
    R[1] = P[1] - Q[1];
    R[2] = P[2] - Q[2];
}

__device__
void mul(double *P, double *Q, double *R){
    R[0] = P[0] * Q[0];
    R[1] = P[1] * Q[1];
    R[2] = P[2] * Q[2];
}

__device__
void mul(double *P, double Q, double *R){
    R[0] = P[0] * Q;
    R[1] = P[1] * Q;
    R[2] = P[2] * Q;
}

__device__
void div(double *P, double Q, double *R){
    R[0] = P[0] / Q;
    R[1] = P[1] / Q;
    R[2] = P[2] / Q;
}

__device__
double magnitude(double *v){
    double sum = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return sqrt(sum);
}

__device__
int factorial(int n){
    int f = 1;
    for (int i = 1; i <= n; ++i) {
        f *= i;
    }
    return f;
}

__device__
int binom(int n, int i){
    if(n < 2){
        return 1;
    }
    return factorial(n)/(factorial(i)*factorial(n-i));
}

__device__
double dot(double *a, double *b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

__device__
void cross(double *a, double *b, double *c){
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
}


__device__
void solve3Dlinear(double *x, double *y, double *z){

//    cout << " Initial linear system: " << endl;
//    cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << endl;
//    cout << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << endl;
//    cout << z[0] << " " << z[1] << " " << z[2] << " " << z[3] << " " << endl;
//    cout << " ---------------------- " << endl;


    if(x[2] == 0){
        if(z[2] != 0){
            double* temp = x;
            x = z;
            z = temp;
        }else{
            double* temp = x;
            x = y;
            y = temp;
        }
    }
    if(y[1] == 0){
        if(z[1] != 0){
            double* temp = y;
            y = z;
            z = temp;
        }else{
            double* temp = y;
            y = x;
            x = temp;
        }
    }
    if(z[0] == 0){
        if(y[0] != 0){
            double* temp = y;
            y = z;
            z = temp;
        }else{
            double* temp = z;
            z = x;
            x = temp;
        }
    }
    double temp = z[2]/x[2];
    for (int i = 0; i < 4; ++i) {
        z[i] -= x[i]*temp;
    }
    temp = y[2]/x[2];
    for (int i = 0; i < 4; ++i) {
        y[i] -= x[i]*temp;
    }
    temp = z[1]/y[1];
    for (int i = 0; i < 4; ++i) {
        z[i] -= y[i]*temp;
    }
    temp = y[0]/z[0];
    for (int i = 0; i < 4; ++i) {
        y[i] -= z[i]*temp;
    }
    temp = x[0]/z[0];
    for (int i = 0; i < 4; ++i) {
        x[i] -= z[i]*temp;
    }
    temp = x[1]/y[1];
    for (int i = 0; i < 4; ++i) {
        x[i] -= y[i]*temp;
    }

//    cout << " Solved linear system: " << endl;
//    cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << endl;
//    cout << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << endl;
//    cout << z[0] << " " << z[1] << " " << z[2] << " " << z[3] << " " << endl;
//    cout << " ---------------------- " << endl;

//    return {z[3]/z[0],y[3]/y[1],x[3]/x[2]};

}

__device__
double sign3D(double *P, double *Q) {
    if(dot(P,Q) > 0){
        return 1.0;
    }else{
        return -1.0;
    }
}

// ------------------------------------------------------------------------------------------------------------------ //

__global__
void orthogonalProjection(gpumodel *mGlobal, double *Parr, double *mindist, double u, double v, double eps = 0.0000001, int iter = 100){
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    printf("at thread %d; threadIdx.x: %d, blockIdx.x %d, blockDim.x %d \n", &index, threadIdx.x, blockIdx.x, blockDim.x);


    gpumodel m = *mGlobal;
    double P[3] = {Parr[3*index + 0], Parr[3*index + 1], Parr[3*index + 2]};

    //printf("%f %f %f\n", &P[0], &P[1], &P[2]);
    //printf("%f %f %f\n", &((m->bendCurve).C0[0]), &((m->bendCurve).C0[1]), &((m->bendCurve).C0[2]));

    double dist = 100000, distOld;
    double lam1 = 1, lam2 = 2;
    double dt = 1;
    int i = 0;

    double p0[3],p0P[3],du[3],duu[3],duv[3],dvv[3],dv[3],n[3],Ps[3],q[3],c1[3],c2[3],cq[3],s[3];
    double g11, g12, g22, h11, h12, h22, g21, h21, k;

    fillTop(p0, u, v, &m);

    sub(P, p0, p0P);

    dist = magnitude(p0P);
    do {
        fillTopDerU(du, u, v, &m);
        fillTopDerV(dv, u, v, &m);
        fillTopDerUU(duu, u, v, &m);
        fillTopDerUV(duv, u, v, &m);
        fillTopDerVV(dvv, u, v, &m);

        cross(du, dv, n);
        div(n, magnitude(n), n);

        double eq1[4] = {du[0], dv[0], n[0], p0P[0]};
        double eq2[4] = {du[1], dv[1], n[1], p0P[1]};
        double eq3[4] = {du[2], dv[2], n[2], p0P[2]};

        solve3Dlinear(eq1, eq2, eq3);

        lam1 = eq3[3]/eq3[0];
        lam2 = eq2[3]/eq2[1];

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


        mul(n, 1.0 / k, n);
        add(p0, n, s);


        sub(P, s, Ps);
        div(Ps, magnitude(Ps) * abs(k), Ps);
        add(s, Ps, q);

        mul(du, lam1, du);
        mul(dv, lam2, dv);
        mul(duu, lam1 * lam1, duu);
        mul(duv, 2 * lam1 * lam2, duv);
        mul(dvv, lam2 * lam2, dvv);
        add(du, dv, c1);
        add(duu, duv, c2);
        add(c2, dvv, c2);


        sub(q, p0, cq);
        double temp[3];
        cross(c1, cq, temp);
        double mag1 = magnitude(temp);
        cross(c1, c2, temp);
        double mag2 = magnitude(temp);
        dt = sign3D(c1, cq) * sqrt(abs( mag1 / mag2 ));

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

        fillTop(p0, u, v, &m);

        sub(P, p0, p0P);
        distOld = dist;
        dist = magnitude(p0P);
    }while(abs(dist - distOld) > eps && i < iter);
    fillTop(p0,u,v,&m);
    sub(P,p0,p0P);
    mindist[index] = magnitude(p0P);
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

    int N = 32;

    gpumodel *m = new gpumodel;

    translateModel(model,m);

    double *Ps, *dists, dists2[N];

    cudaMallocManaged(&Ps, 3*N*sizeof(double));
    cudaMallocManaged(&dists, N*sizeof(double));

    for (int i = 0; i < N; i++) {
        vector<double> P = model.fillTop(double(rand())/RAND_MAX,double(rand())/RAND_MAX);

        P[0] += (double(rand())/RAND_MAX)*0.5;
        P[1] += (double(rand())/RAND_MAX)*0.5;
        P[2] += (double(rand())/RAND_MAX)*0.5;

        Ps[i*3 + 0] = P[0];
        Ps[i*3 + 1] = P[1];
        Ps[i*3 + 2] = P[2];

        dists2[i] = gridSearch(model,P);
    }

    for (int k = 0; k < 1000; k+=100) {
        cout << Ps[k*3 + 0] << " " << Ps[k*3 + 1] << " " << Ps[k*3 + 2] << endl;
    }

    orthogonalProjection<<<1, 1>>>(m,Ps,dists,0.5,0.5);

    cudaDeviceSynchronize();

    for (int k = 0; k < 1000; k+=100) {
        if(dists[k] < dists2[k])
            cout << "True ";
        else
            cout << "False ";
        cout << "Dist " << k << ": " << dists[k] << endl;
    }
    cudaFree(Ps);
    cudaFree(dists);
    delete m;
}