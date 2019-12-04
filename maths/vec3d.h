//
// Created by s3179222 on 12/4/19.
//

#ifndef HEARTVALVEMODEL_VEC3D_H
#define HEARTVALVEMODEL_VEC3D_H

#include <cassert>
#include <cmath>

class vec3d {
private:
    double x,y,z;
public:
    vec3d(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    vec3d operator + (vec3d b);
    vec3d operator - (vec3d b);
    vec3d operator * (vec3d b);
    vec3d operator / (vec3d b);
    vec3d operator + (double b);
    vec3d operator - (double b);
    vec3d operator * (double b);
    vec3d operator / (double b);
    int operator == (vec3d b);
    vec3d operator - ();
    double dot (vec3d b);
    vec3d cross (vec3d b);
    double mag ();
    double sqdist (vec3d b);
    double dist (vec3d b);
    void non_zero();
};


#endif //HEARTVALVEMODEL_VEC3D_H
