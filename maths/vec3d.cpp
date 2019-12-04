//
// Created by s3179222 on 12/4/19.
//

#include "vec3d.h"

vec3d vec3d::operator + (vec3d b){
    return {x + b.x, y + b.y, z + b.z};
};
vec3d vec3d::operator - (vec3d b){
    return {x - b.x, y - b.y, z - b.z};

};
vec3d vec3d::operator * (vec3d b){
    return {x * b.x, y * b.y, z * b.z};

};
vec3d vec3d::operator / (vec3d b){
    b.non_zero();
    return {x / b.x, y / b.y, z / b.z};

};
vec3d vec3d::operator + (double b){
    return {x + b, y + b, z + b};
};
vec3d vec3d::operator - (double b){
    return {x - b, y - b, z - b};

};
vec3d vec3d::operator * (double b){
    return {x * b, y * b, z * b};

};
vec3d vec3d::operator / (double b){
    return {x / b, y / b, z / b};
};
int vec3d::operator == (vec3d b){
    return (x == b.x && y == b.y && z == b.z);
};
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