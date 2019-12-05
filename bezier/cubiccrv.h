//
// Created by s3179222 on 12/4/19.
//

#ifndef HEARTVALVEMODEL_CUBICCRV_H
#define HEARTVALVEMODEL_CUBICCRV_H

#include "../maths/vec3d.h"

class cubiccrv {
private:
    vec3d c0{},c1{},c2{},c3{};

public:

    cubiccrv(vec3d c0, vec3d c1, vec3d c2, vec3d c3){
        this->c0 = c0;
        this->c1 = c1;
        this->c2 = c2;
        this->c3 = c3;
    }

    cubiccrv(vec3d c0, vec3d c1, vec3d c2){
        this->c0 = c0;
        this->c1 = c0 * 1.0/3.0 + c1 * 2.0/3.0;
        this->c2 = c1 * 2.0/3.0 + c2 * 1.0/3.0;
        this->c3 = c2;
    }

    vec3d f(double t);
    vec3d df(double t);
    vec3d ddf(double t);
};


#endif //HEARTVALVEMODEL_CUBICCRV_H
