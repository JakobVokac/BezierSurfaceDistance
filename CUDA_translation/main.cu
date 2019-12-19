#include <iostream>
#include <chrono>
#include "../model.h"
#include "full_optimizer.cu"
#include "top_parametric_functions.cu"
#include <cmath>





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