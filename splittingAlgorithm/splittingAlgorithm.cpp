//
// Created by s3179222 on 12/16/19.
//

#include "splittingAlgorithm.h"

OptState2D splittingAlgorithm::optimize(bicubicsrf s, vec3d P, double eps, int recursiveLimit) {

    bicubicsrf sf;
    double uf, vf;
    double dist, distMin;
    if(s.closestPointInPatch(P)){
        surfaceSplitting(1, 0, 1, 0, 1, uf, vf, s, sf, P, eps, 0, recursiveLimit);
        distMin = s.distTo(uf,vf,P);
    }else{
        cubiccrv *u0,*v0,*u1,*v1, cf;
        u0 = dynamic_cast<cubiccrv*>(&s.edgeU0());
        v0 = dynamic_cast<cubiccrv*>(&s.edgeV0());
        u1 = dynamic_cast<cubiccrv*>(&s.edgeU1());
        v1 = dynamic_cast<cubiccrv*>(&s.edgeV1());
        if(u0->closestPointInCurve(P)){
            curveSplitting(0,1,vf,*u0,cf,P,eps,0,recursiveLimit);
            uf = 0;
            distMin = u0->distTo(vf,P);
        }else if(v0->closestPointInCurve(P)){
            curveSplitting(0,1,uf,*v0,cf,P,eps,0,recursiveLimit);
            vf = 0;
            distMin = v0->distTo(uf,P);
        }else if(u1->closestPointInCurve(P)){
            curveSplitting(0,1,vf,*u1,cf,P,eps,0,recursiveLimit);
            uf = 1;
            distMin = u1->distTo(vf,P);
        }else if(v1->closestPointInCurve(P)){
            curveSplitting(0,1,uf,*v1,cf,P,eps,0,recursiveLimit);
            vf = 1;
            distMin = v1->distTo(uf,P);
        }else{
            distMin = s.distTo(0,0,P);
            uf = 0, vf = 0;
            dist = s.distTo(0,1,P);
            if(dist < distMin) {
                distMin = dist;
                uf = 0, vf = 1;
            }
            dist = s.distTo(1,0,P);
            if(dist < distMin) {
                distMin = dist;
                uf = 1, vf = 0;
            }
            dist = s.distTo(1,1,P);
            if(dist < distMin) {
                distMin = dist;
                uf = 1, vf = 1;
            }
        }
    }
    return {uf,vf,distMin};
}

void splittingAlgorithm::surfaceSplitting(int dir, double ul, double ur, double vl, double vr, double &uf, double &vf,
                                          bicubicsrf s, bicubicsrf &sf, vec3d P, double eps, int i, int limit) {
    bool valid = s.hasValidControlNet();
    if(valid){
        if(!s.closestPointInPatch(P)){
            return;
        }else{
            vec3d diag = s.ctrlP(0) - s.ctrlP(15);
            if(diag.mag() < eps || i > limit){
                sf = s;
                uf = (ul+ur)/2;
                vf = (vl+vr)/2;
                return;
            }
        }
    }
    bicubicsrf s1, s2;

    s.subdivideInDir(dir, 0.5, s1, s2);

//    cout << "uvlr: " << ul << " " << ur << " " << vl << " " << vr << " i: " << i << endl;
//    plotSurface(s1,0.1);
//    plotSurface(s2,0.1);
    if(!dir){
        surfaceSplitting((dir + 1)%2,ul,(ul+ur)/2,vl,vr,uf,vf,s1,sf,P,eps,i+1,limit);
        surfaceSplitting((dir + 1)%2,(ul+ur)/2,ur,vl,vr,uf,vf,s2,sf,P,eps,i+1,limit);
    }else{
        surfaceSplitting((dir + 1)%2,ul,ur,vl,(vl+vr)/2,uf,vf,s1,sf,P,eps,i+1,limit);
        surfaceSplitting((dir + 1)%2,ul,ur,(vl+vr)/2,vr,uf,vf,s2,sf,P,eps,i+1,limit);
    }
}

void splittingAlgorithm::curveSplitting(double tl, double tr, double &tf, cubiccrv c, cubiccrv &cf, vec3d P, double eps,
                                        int i, int limit) {
    bool valid = c.hasConvexPolygon();
    if(valid){
        if(!c.closestPointInCurve(P)){
            return;
        }else{
            vec3d c0c3 = c.getCtrlP(3) - c.getCtrlP(0);
            if(c0c3.mag() < eps || i > limit){
                cf = c;
                tf = (tl+tr)/2;
                return;
            }
        }
    }
    cubiccrv c1,c2;
    c.subdivide(0.5, c1, c2);
    curveSplitting(tl, (tl+tr)/2, tf, c1, cf, P, eps, i + 1, limit);
    curveSplitting((tl+tr)/2, tr, tf, c2, cf, P, eps, i + 1, limit);
}


