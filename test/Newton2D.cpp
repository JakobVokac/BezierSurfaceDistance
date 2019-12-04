//
// Created by s3179222 on 12/3/19.
//
#include "Newton2D.h"

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
double NewtonMethodTest2D(double x, double y, double eps=0.01, double sigma = 0.1, int plot = 0){

    double V = f(x,y);
    double V_old;
    int iter = 0;

    vector<double> xs, ys;
    if(plot) {
        xs.push_back(x);
        ys.push_back(y);
    }
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
        if(plot) {
            xs.push_back(x);
            ys.push_back(y);
        }
    }while(abs(V_old - V) > eps);
    if(plot)
        plt::plot(xs,ys);
    return V;
}