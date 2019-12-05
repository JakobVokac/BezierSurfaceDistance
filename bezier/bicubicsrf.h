//
// Created by s3179222 on 12/4/19.
//

#ifndef HEARTVALVEMODEL_BICUBICSRF_H
#define HEARTVALVEMODEL_BICUBICSRF_H

#include "../maths/vec3d.h"

class bicubicsrf {
private:
    vec3d ctrl[16]{};
public:

    explicit bicubicsrf(vec3d ctrl[]){
        for (int i = 0; i < 16; ++i) {
            this->ctrl[i] = ctrl[i];
        }
    }

    vec3d f(double u, double v);
    vec3d fdu(double u, double v);
    vec3d fdv(double u, double v);
    vec3d fdudu(double u, double v);
    vec3d fdvdv(double u, double v);
    vec3d fdudv(double u, double v);
};


#endif //HEARTVALVEMODEL_BICUBICSRF_H
