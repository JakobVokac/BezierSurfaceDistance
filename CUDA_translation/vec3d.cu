//
// Created by s3179222 on 12/4/19.
//

#include <cstdlib>
#include "vec3d.h"

__host__ __device__ vec3d operator + (vec3d a, const vec3d& b){
    a += b;
    return a;
};

__host__ __device__ vec3d& vec3d::operator += (const vec3d& b)
{
    this->x += b.x;
    this->y += b.y;
    this->z += b.z;

    return *this;
}
__host__ __device__ vec3d operator - (vec3d a, const vec3d& b){
    a -= b;
    return a;
};

__host__ __device__ vec3d& vec3d::operator -= (const vec3d& b)
{
    this->x -= b.x;
    this->y -= b.y;
    this->z -= b.z;

    return *this;
}

__host__ __device__ vec3d operator * (vec3d a, const vec3d& b){
    a *= b;
    return a;
};

__host__ __device__ vec3d& vec3d::operator *= (const vec3d& b)
{
    this->x *= b.x;
    this->y *= b.y;
    this->z *= b.z;

    return *this;
}

__host__ __device__ vec3d operator / (vec3d a, const vec3d& b){
    a /= b;
    return a;
};

__host__ __device__ vec3d& vec3d::operator /= (vec3d b)
{
    b.non_zero();
    this->x /= b.x;
    this->y /= b.y;
    this->z /= b.z;

    return *this;
}

__host__ __device__ vec3d operator + (vec3d a, const double& b){
    a += b;
    return a;
};
__host__ __device__ vec3d operator + (const double& b, vec3d a){
    a += b;
    return a;
};
__host__ __device__ vec3d& vec3d::operator += (const double& b)
{
    this->x += b;
    this->y += b;
    this->z += b;

    return *this;
}
__host__ __device__ vec3d operator - (vec3d a, const double& b){
    a -= b;
    return a;
};
__host__ __device__ vec3d operator - (const double& b, vec3d a){
    a -= b;
    return a;
};
__host__ __device__ vec3d& vec3d::operator -= (const double& b)
{
    this->x -= b;
    this->y -= b;
    this->z -= b;

    return *this;
}
__host__ __device__ vec3d operator * (vec3d a, const double& b){
    a *= b;
    return a;
};
__host__ __device__ vec3d operator * (const double& b, vec3d a){
    a *= b;
    return a;
};
__host__ __device__ vec3d& vec3d::operator *= (const double& b)
{
    this->x *= b;
    this->y *= b;
    this->z *= b;

    return *this;
}
__host__ __device__ vec3d operator / (vec3d a, const double& b){
    a /= b;
    return a;
};
__host__ __device__ vec3d operator / (const double& b, vec3d a){
    a /= b;
    return a;
};
__host__ __device__ vec3d& vec3d::operator /= (const double& b){
    assert(b != 0);
    this->x /= b;
    this->y /= b;
    this->z /= b;

    return *this;
}
__host__ __device__ vec3d vec3d::operator - (){
    return {-x, -y, -z};

};
__host__ __device__ vec3d vec3d::operator + (){
    return *this;

};
__host__ __device__ double vec3d::dot (vec3d b){
    return x * b.x + y * b.y + z * b.z;

};
__host__ __device__ vec3d vec3d::cross (vec3d b){
    return {y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x};
};

__host__ __device__ double vec3d::mag (){
    return sqrt(x*x+y*y+z*z);
};
__host__ __device__ double vec3d::sqdist (vec3d b){
    return (x-b.x)*(x-b.x)+(y-b.y)*(y-b.y)+(z-b.z)*(z-b.z);
};

__host__ __device__ double vec3d::dist (vec3d b){
    return sqrt(sqdist(b));
};

__host__ __device__ void vec3d::non_zero(){
    assert(x != 0);
    assert(y != 0);
    assert(z != 0);
};

__host__ __device__ double vec3d::sum(){
    return x+y+z;
}

__host__ __device__ vec3d vec3d::RotAxisZ(double rad){
    return {
        x * cos(rad) - y * sin(rad) + 0,
        x * sin(rad) + y * cos(rad) + 0,
        0            + 0            + z
    };
}

__host__ __device__ vec3d vec3d::MirrorY() {
    return {x,-y,z};
}

__host__ __device__ double vec3d::getx() {
    return x;
}
__host__ __device__ double vec3d::gety() {
    return y;
}
__host__ __device__ double vec3d::getz() {
    return z;
}

__host__ __device__ bool vec3d::operator == (const vec3d& b) const{
    return x == b.x && y == b.y && z == b.z;
}

__host__ __device__ bool vec3d::operator != (const vec3d& b) const{
    return !(*this == b);
}

__host__ __device__ double curvature(vec3d c1, vec3d c2) {
    double mag = c1.mag();
    return (c1.cross(c2)).mag()/(mag*mag*mag);
}

__host__ __device__ vec3d circleCenterDir(vec3d c1, vec3d c2) {
    vec3d c = c1.cross(c2);
    c = c.cross(c1);
    return c/c.mag();
}

__host__ __device__ double sign(vec3d P, vec3d Q) {
    if(P.dot(Q) > 0){
        return 1.0;
    }else{
        return -1.0;
    }
};

__host__ __device__ vec3d solve3Dlinear(vec3d v1, vec3d v2, vec3d v3, vec3d v4){

    double *x, *y, *z;
    x = new double[4]{v1.getx(),v2.getx(),v3.getx(),v4.getx()};
    y = new double[4]{v1.gety(),v2.gety(),v3.gety(),v4.gety()};
    z = new double[4]{v1.getz(),v2.getz(),v3.getz(),v4.getz()};

    if(x[2] == 0){
        if(z[2] != 0){
            double *temp = x;
            x = z;
            z = temp;
        }else{
            double *temp = x;
            x = y;
            y = temp;
        }
    }
    if(y[1] == 0){
        if(z[1] != 0){
            double *temp = y;
            y = z;
            z = temp;
        }else{
            double *temp = y;
            y = x;
            x = temp;
        }
    }
    if(z[0] == 0){
        if(y[0] != 0){
            double *temp = y;
            y = z;
            z = temp;
        }else{
            double *temp = z;
            z = x;
            x = temp;
        }
    }
    double temp;

    temp = z[2]/x[2];
    for(int i = 0; i < 4; i++){
        z[i] -= x[i]*temp;
    }
    temp = y[2]/x[2];
    for(int i = 0; i < 4; i++){
        y[i] -= x[i]*temp;
    }
    temp = z[1]/y[1];
    for(int i = 0; i < 4; i++){
        z[i] -= y[i]*temp;
    }
    temp = y[0]/z[0];
    for(int i = 0; i < 4; i++){
        y[i] -= z[i]*temp;
    }
    temp = x[0]/z[0];
    for(int i = 0; i < 4; i++){
        x[i] -= z[i]*temp;
    }
    temp = x[1]/y[1];
    for(int i = 0; i < 4; i++){
        x[i] -= y[i]*temp;
    }

    vec3d res = {z[3]/z[0],y[3]/y[1],x[3]/x[2]};

    delete x;
    delete y;
    delete z;

    return res;

}

__host__ std::ostream &operator<<(std::ostream &out, vec3d &v) {
    out << "x: " << v.x << " y: " << v.y << " z: " << v.z;
    return out;
}
