//
// Created by s3179222 on 12/10/19.
//
/*
 * This class is the implementation of the Geometric method for the numerical optimizer.
 *
 * Source: Shi-min Hu and Johannes Wallner,
 *        "A second order algorithm for orthogonal projection onto curves and surfaces",
 *        http://www.geometrie.tugraz.at/wallner/sproj.pdf
 */
#include "Geometric.h"
using namespace std;

OptState2D Geometric::doStep(surface &sur, vec3d &P, OptState2D &loc) {

    double u = loc.u, v = loc.v;

    vec3d p0{},p0P{},du{},duu{},duv{},dvv{},dv{},n{},Ps{},q{},c1{},c2{},cq{},s{},coefs{};
    double g11, g12, g22, h11, h12, h22, g21, h21, k, lam1, lam2, dt;

    p0 = sur.at(u, v);
    p0P = P - p0;

    du = sur.atDerU(u, v);
    dv = sur.atDerV(u, v);
    duu = sur.atDerUU(u, v);
    dvv = sur.atDerVV(u, v);
    duv = sur.atDerUV(u, v);

    n = du.cross(dv);
    n /= n.mag();

    coefs = solve3Dlinear(du,dv,n,p0P);
    lam1 = coefs.getx();
    lam2 = coefs.gety();

    g11 = du.dot(du);
    g12 = du.dot(dv);
    g21 = dv.dot(du);
    g22 = dv.dot(dv);
    h11 = duu.dot(n);
    h12 = duv.dot(n);
    h21 = h12;
    h22 = dvv.dot(n);

    k = (h11 * lam1 * lam1 + h12 * lam1 * lam2 + h21 * lam2 * lam1 + h22 * lam2 * lam2) /
        (g11 * lam1 * lam1 + g12 * lam1 * lam2 + g21 * lam2 * lam1 + g22 * lam2 * lam2);

    s = p0 + n / k;

    Ps = P - s;
    q = s + Ps / (Ps.mag() * abs(k));

    c1 = du * lam1 + dv * lam2;
    c2 = duu * lam1 * lam1 + duv * 2 * lam1 * lam2 + dvv * lam2 * lam2;

    cq = q - p0;
    dt = sign(c1, cq) * sqrt(abs(2 * (c1.cross(cq)).mag() / (c1.cross(c2)).mag()));

    loc.u += dt * lam1;
    loc.v += dt * lam2;

    p0 = sur.at(loc.u, loc.v);
    p0P = P - p0;
    loc.dist = p0P.mag();

    return loc;
}

OptState1D Geometric::doStep(curve &crv, vec3d &P, OptState1D &loc) {
    double dt = 0;
    double t = loc.t;

    vec3d p0{}, p0P{}, q{}, c1{}, c2{}, cq{}, c{}, s{}, u{};
    double k;

    p0 = crv.f(t);

    c1 = crv.df(t);
    c2 = crv.ddf(t);

    k = curvature(c1, c2);
    c = circleCenterDir(c1, c2);
    s = p0 + c / k;
    u = P - s;
    q = s + u / (k * u.mag());
    cq = q - p0;
    dt = sign(c1, cq) * sqrt(abs(2 * (c1.cross(cq)).mag() / (c1.cross(c2)).mag()));

    t += dt;

    p0 = crv.f(t);

    p0P = P - p0;

    loc.t = t;
    loc.dist = p0P.mag();

    return loc;

}
