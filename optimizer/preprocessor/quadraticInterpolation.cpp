//
// Created by s3179222 on 12/10/19.
//

#include "quadraticInterpolation.h"

OptState2D quadraticInterpolation::preprocess(surface &sur, const vec3d &P) {
    double u1 = 0.3, u2 = 0.6, u3 = 0.9, v1 = 0, v2 = 0.5, v3 = 1;
    double h1,h2,h3,h4;
    double unew = 0.5, vnew = 0.5;
    for (int i = 0; i < iterations; ++i) {

        if(unew != 0 && unew != 1) {

            h1 = sur.sqDistTo(u1, vnew, P);
            h2 = sur.sqDistTo(u2, vnew, P);
            h3 = sur.sqDistTo(u3, vnew, P);

            unew = (u1*u1*(h2-h3) + u2*u2*(h3-h1) + u3*u3*(h1-h2))/(2.0*(u1*(h2-h3) + u2*(h3-h1) + u3*(h1-h2)));

            if (unew < 0) {
                unew = 0;
            } else if (unew > 1) {
                unew = 1;
            }

            h4 = sur.sqDistTo(unew, vnew, P);

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

        if(vnew != 0 && vnew != 1) {

            h1 = sur.sqDistTo(unew, v1, P);
            h2 = sur.sqDistTo(unew, v2, P);
            h3 = sur.sqDistTo(unew, v3, P);

            vnew = 0.5 * ((v2 * v2 - v3 * v3) * h1 + (v3 * v3 - v1 * v1) * h2 + (v1 * v1 - v2 * v2) * h3) /
                   ((v2 - v3) * h1 + (v3 - v1) * h2 + (v1 - v2) * h3);

            if (vnew < 0) {
                vnew = 0;
            } else if (vnew > 1) {
                vnew = 1;
            }

            h4 = sur.sqDistTo(unew, vnew, P);

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
    }
    return {unew,unew,sur.distTo(unew,vnew,P)};
}

OptState1D quadraticInterpolation::preprocess(curve &crv, const vec3d &P) {
    double h1, h2, h3, hn;
    double t1 = 0, t2 = 0.5, t3 = 1, tn = 0;

    for (int i = 0; i < iterations; i++) {

        h1 = crv.sqDistTo(t1,P);
        h2 = crv.sqDistTo(t2,P);
        h3 = crv.sqDistTo(t3,P);

        if(t1 == t2 || t2 == t3 || t3 == t1){
            break;
        }

        tn = (t1*t1*(h2-h3) + t2*t2*(h3-h1) + t3*t3*(h1-h2))/(2.0*(t1*(h2-h3) + t2*(h3-h1) + t3*(h1-h2)));

        if (tn < 0) {
            tn = 0;
        } else if (tn > 1) {
            tn = 1;
        }

        hn = crv.sqDistTo(tn,P);

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

    if(tn > 1)
        tn = 1;
    if(tn < 0)
        tn = 0;

    return {tn,crv.distTo(tn,P)};
}
