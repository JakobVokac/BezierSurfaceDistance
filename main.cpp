#include <iostream>
#include <vector>
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


double gradientDescent(Model model, vector<double> A, double u, double v, double eps=0.01, double gamma = 0.1){
    double dist = model.distanceToTopPoint(u,v,A);
    double dist_old;

    int iter = 0;

    vector<double> us, vs;
    us.push_back(u);
    vs.push_back(v);

    do{
//        cout << "u: " << u << " v: " << v << " A: x:" << A[0] << " y:" << A[1] << " z:" << A[2] << endl;
//        cout << "Iter: " << iter << " Delta u: " << model.distanceToTopPointDerU(u,v,A) << " Delta v: " << model.distanceToTopPointDerV(u,v,A) << endl;
        u -= gamma * model.completeDistanceTopDerU(u,v,A);
        v -= gamma * model.completeDistanceTopDerV(u,v,A);

        dist_old = dist;
        dist = model.distanceToTopPoint(u,v,A);
        us.push_back(u);
        vs.push_back(v);
        iter++;
    }while(abs(dist_old - dist) > eps);
    plt::plot(us,vs);
//    cout << "Gradient descent iterations: " << iter << " u: " << u << " v: " << v << endl;
    return dist;
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

double NewtonForLinearCG(Model model, vector<double> A, double u, double v, double dU, double dV, double eps, double HModifier = 0){
    double V = model.distanceToTopPoint(u,v,A);
    double V_old;
    double t = 0;

//    vector<double> xs,ys, ys2, ys3;
/*
    cout << "Fletcher Reeves Newton" << endl;
    cout << u << " " << dU << " " << -u/dU << " " << (1-u)/dU << (1-2*u)/dU * 0.01 << endl;
    for(double i = 0; i < 1; i+= 0.01 ){
        xs.push_back(i);
        ys.push_back( model.distanceToTopPoint(u + dU*((1-2*u)/dU * i), v + dV*((1-2*u)/dU * i), A));
        ys2.push_back((dU * model.completeDistanceTopDerU(u + dU*((1-2*u)/dU * i),v + dV*((1-2*u)/dU * i),A) + dV * model.completeDistanceTopDerV(u + dU*((1-2*u)/dU * i),v + dV*((1-2*u)/dU * i),A))/100);
        ys3.push_back((dU*dU*model.completeDistanceTopDerUU2(u + dU*((1-2*u)/dU * i),v + dV*((1-2*u)/dU * i),A)
                       + 2*dU*dV*model.completeDistanceTopDerUV2(u + dU*((1-2*u)/dU * i),v + dV*((1-2*u)/dU * i),A)
                       + dV*dV*model.completeDistanceTopDerVV2(u + dU*((1-2*u)/dU * i),v + dV*((1-2*u)/dU * i),A))/10000);
        cout << model.distanceToTopPoint(u + dU*((1-2*u)/dU * i), v + dV*((1-2*u)/dU * i), A) << endl;
    }
    plt::plot(xs,ys,"g");
    plt::plot(xs,ys2,"b");
    plt::plot(xs,ys3,"r");

    plt::show();
*/
    do{
        t -= (dU * model.completeDistanceTopDerU(u + dU*t,v + dV*t,A) + dV * model.completeDistanceTopDerV(u + dU*t,v + dV*t,A))
                /(dU*dU*model.completeDistanceTopDerUU2(u + dU*t,v + dV*t,A)
                + 2*dU*dV*model.completeDistanceTopDerUV2(u + dU*t,v + dV*t,A)
                + dV*dV*model.completeDistanceTopDerVV2(u + dU*t,v + dV*t,A));
//        cout << "t: " << t << endl;
        V_old = V;
        V = model.distanceToTopPoint(u + dU*t, v + dV*t, A);
    }while(abs(V-V_old) > eps && 0 <= u + dU*t && u + dU*t <= 1 && v + dV*t <= 1 && 0 <= v + dV*t);

    return t;
}

vector<double> FletcherReevesCG(Model model, vector<double> A, double u, double v, double eps=0.0001, int plot = 0, double HModifier = 0){
    double dU, dV, alpha, beta, dU2, dV2, u0 = u, v0 = v;
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
        alpha = NewtonForLinearCG(model,A,u,v,dU,dV,0.0000001);
//        cout << "alpha: " << alpha << endl;
//        cout << "dU: " << dU << " dV: " << dV << endl;
        if(plot){
            cout << "u: " << u << " v: " << v << " dU: " << dU << " dV: " << dV << " alpha: " << alpha << endl;
        }
        u += alpha * dU;
        v += alpha * dV;

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
    return {u,v,model.distanceToTopPoint(u,v,A), u0, v0};
}

vector<double> findStartingPoint(Model model, vector<double> A){
    double ul = 0, ur = 1, vl = 0, vr = 1;
    double dist = model.distanceToTopPoint(0.5,0.5,A);
    double um, vm;
    for (int i = 0; i < 12; ++i) {
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
    return {um,vm};
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
//    drawPart(p2, x, y, z);
//
//    drawPart(p3, x, y, z);
//    drawPart(p4, x, y, z);
//
//    drawPart(p5, x, y, z);
//    drawPart(p6, x, y, z);

    plt::plot_surface(x,y,z);
    plt::xlim(-20,20);
    plt::ylim(-20,20);

    plt::show();


    vector<double> A = {8,6.3,10};
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

           /* res = FletcherReevesCG(p1,A,j,i,0.00001);
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
            }*/
        }
        /*CGus.push_back(CGusR);
        CGvs.push_back(CGvsR);
        CGVs.push_back(CGVsR);*/
        Nus.push_back(NusR);
        Nvs.push_back(NvsR);
        NVs.push_back(NVsR);
    }
    plt::plot_surface(Nus,Nvs,NVs);
    plt::show();
   /* plt::plot_surface(CGus,CGvs,CGVs);
    plt::show();*/
    cout << "Distance via Newton's Method for 2D: " << NewtonMethod(p1,A,0.8,0.2,0.00001,.5)[2] << endl;
    plt::plot_surface_with_vector_field(u,v,dist,x2,y2,z2,u2,v2,w2,0.03);
    plt::show();
    cout << "Distance via grid search: " << min_dist << endl;
    cout << "Distance via Fletcher-Reeves CG: " << FletcherReevesCG(p1,A,0.8,0.2,0.00001)[2] << endl;


    plt::xlim(0,1);
    plt::ylim(0,1);
    plt::show();

    vector<double> P = {0,0,0};

    for(double i = 4; i < 10; i += 0.2){
        for(double j = 4; j < 10; j += 0.2){
            for(double k = 4; k < 10; k += 0.2){
                P = {i,j,k};
                cout << "Point: " << i << ", " << j << ", " << k << endl;
                vector<double> start = findStartingPoint(p1,P);
                vector<double> res;
                res = NewtonMethod(p1,P,start[0],start[1],0.00001,.5);
                cout << "Newton u: " << res[0] << " v: " << res[1] << " V: " << res[2] << " u0: " << res[3] << " v0: " << res[4] << " iterations: " << res[5] << endl;

            }
        }
    }
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



