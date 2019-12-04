//
// Created by s3179222 on 12/3/19.
//
#include "GeometricProjectionSur.h"

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

//        cout << "curvature2D: " << k << " radius: " << 1.0/k << endl;

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