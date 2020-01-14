#include "vec3d.h"

struct QBezierGPU{
    vec3d C0, C1, C2;
};

struct CBezierGPU{
    vec3d C0, C1, C2, C3;
};

struct gpumodeltop{
    CBezierGPU leafCurve, bendCurve, symCurveTop;
    QBezierGPU sinCurveTop;
    vec3d Q_b_sym, Q_l_sym, Q_b_sin, Q_l_sin;
};