//
// Created by s3179222 on 12/4/19.
//
/*
 * This is a simple implementation of a vector with doubles for higher precision.
 * It overrides all of the basic operators and adds the dot product, cross product, an output function for cout.
 *
 * There is also an implementation for solving a system of 3 linear equations, which is useful for finding the
 * parameters in a linear combination of vectors for another vector.
 */

#ifndef HEARTVALVEMODEL_VEC3D_H
#define HEARTVALVEMODEL_VEC3D_H

#include <cassert>
#include <cmath>
#include <ostream>

class vec3d {
public:
    double x,y,z;

    vec3d () = default;

    vec3d(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    vec3d& operator += (const vec3d& b);

    vec3d& operator -= (const vec3d& b);

    vec3d& operator *= (const vec3d& b);

    vec3d& operator /= (vec3d b);

    vec3d& operator += (const double& b);

    vec3d& operator -= (const double& b);

    vec3d& operator *= (const double& b);

    vec3d& operator /= (const double& b);
    vec3d operator + ();
    vec3d operator - ();
    double dot (vec3d b);
    vec3d cross (vec3d b);
    double mag ();
    double sum ();
    double sqdist (vec3d b);
    double dist (vec3d b);
    void non_zero();

    vec3d RotAxisZ(double rad);
    vec3d MirrorY();

    double getx();
    double gety();
    double getz();

    bool operator == (const vec3d& b) const;
    bool operator != (const vec3d& b) const;

    friend std::ostream & operator << (std::ostream &out, vec3d &v);
};

vec3d operator + (vec3d a, const vec3d& b);
vec3d operator - (vec3d a, const vec3d& b);
vec3d operator * (vec3d a, const vec3d& b);
vec3d operator / (vec3d a, const vec3d& b);
vec3d operator + (vec3d a, const double& b);
vec3d operator - (vec3d a, const double& b);
vec3d operator * (vec3d a, const double& b);
vec3d operator / (vec3d a, const double& b);
vec3d operator + (const double& b, vec3d a);
vec3d operator - (const double& b, vec3d a);
vec3d operator * (const double& b, vec3d a);
vec3d operator / (const double& b, vec3d a);

double curvature(vec3d c1, vec3d c2);
vec3d circleCenterDir(vec3d c1,vec3d c2);
double sign(vec3d P, vec3d Q);
vec3d solve3Dlinear(vec3d v1, vec3d v2, vec3d v3, vec3d v4);


#endif //HEARTVALVEMODEL_VEC3D_H
