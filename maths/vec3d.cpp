//
// Created by s3179222 on 12/4/19.
//

#include "vec3d.h"

vec3d operator + (vec3d a, const vec3d& b){
    a += b;
    return a;
};

vec3d& vec3d::operator += (const vec3d& b)
{
    this->x += b.x;
    this->y += b.y;
    this->z += b.z;

    return *this;
}
vec3d operator - (vec3d a, const vec3d& b){
    a -= b;
    return a;
};

vec3d& vec3d::operator -= (const vec3d& b)
{
    this->x -= b.x;
    this->y -= b.y;
    this->z -= b.z;

    return *this;
}

vec3d operator * (vec3d a, const vec3d& b){
    a *= b;
    return a;
};

vec3d& vec3d::operator *= (const vec3d& b)
{
    this->x *= b.x;
    this->y *= b.y;
    this->z *= b.z;

    return *this;
}

vec3d operator / (vec3d a, const vec3d& b){
    a /= b;
    return a;
};

vec3d& vec3d::operator /= (const vec3d& b)
{
    this->x /= b.x;
    this->y /= b.y;
    this->z /= b.z;

    return *this;
}

vec3d operator + (vec3d a, const double& b){
    a += b;
    return a;
};
vec3d operator + (const double& b, vec3d a){
    a += b;
    return a;
};
vec3d& vec3d::operator += (const double& b)
{
    this->x += b;
    this->y += b;
    this->z += b;

    return *this;
}
vec3d operator - (vec3d a, const double& b){
    a -= b;
    return a;
};
vec3d operator - (const double& b, vec3d a){
    a -= b;
    return a;
};
vec3d& vec3d::operator -= (const double& b)
{
    this->x -= b;
    this->y -= b;
    this->z -= b;

    return *this;
}
vec3d operator * (vec3d a, const double& b){
    a *= b;
    return a;
};
vec3d operator * (const double& b, vec3d a){
    a *= b;
    return a;
};
vec3d& vec3d::operator *= (const double& b)
{
    this->x *= b;
    this->y *= b;
    this->z *= b;

    return *this;
}
vec3d operator / (vec3d a, const double& b){
    a /= b;
    return a;
};
vec3d operator / (const double& b, vec3d a){
    a /= b;
    return a;
};
vec3d& vec3d::operator /= (const double& b)
{
    this->x /= b;
    this->y /= b;
    this->z /= b;

    return *this;
}
vec3d vec3d::operator - (){
    return {-x, -y, -z};

};
double vec3d::dot (vec3d b){
    return x * b.x + y * b.y + z * b.y;

};
vec3d vec3d::cross (vec3d b){
    return {y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.z};

};
double vec3d::mag (){
    return sqrt(x*x+y*y+z*z);
};
double vec3d::sqdist (vec3d b){
    return x*b.x+y*b.y+z*b.z;
};

double vec3d::dist (vec3d b){
    return sqrt(sqdist(b));
};

void vec3d::non_zero(){
    assert(x != 0);
    assert(y != 0);
    assert(z != 0);
};