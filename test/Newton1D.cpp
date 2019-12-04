//
// Created by s3179222 on 12/3/19.
//
#include "Newton1D.h"

double func(double x){
    return pow((x-2.5),2.0) + 2*x + pow(x, 3.0);
}
double funcdx(double x) {
    return 2 * (x-2.5) + 2 + 3*pow(x, 2.0);
}
double funcdxdx(double x){
    return 2 + 6*x;
};
double NewtonMethodTest1D(double x, double eps, int plot){
    double V = func(x);
    double V_old;
    double t = 0;

    vector<double> xs,ys, ys2, ys3;
    if(plot) {
        cout << "Fletcher Reeves Newton Test" << endl;
        for (double i = 0; i < 1; i += 0.01) {
            xs.push_back(i);
            ys.push_back(func(i));
            ys2.push_back(funcdx(i));
            ys3.push_back(funcdxdx(i));
        }
        plt::plot(xs, ys, "g");
        plt::plot(xs, ys2, "b");
        plt::plot(xs, ys3, "r");

        plt::show();
    }
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