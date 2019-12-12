//
// Created by s3179222 on 12/4/19.
//

#ifndef HEARTVALVEMODEL_BICUBICSRF_H
#define HEARTVALVEMODEL_BICUBICSRF_H

#include "surface.h"

class bicubicsrf : public surface{
private:
    vec3d ctrl[16]{};
    cubiccrv U0, V0, U1, V1;
public:
    bicubicsrf() = default;
    explicit bicubicsrf(vec3d ctrl[16]){
        for (int i = 0; i < 16; ++i) {
            this->ctrl[i] = ctrl[i];
        }
        this->U0 = cubiccrv(ctrl[0],ctrl[1],ctrl[2],ctrl[3]);
        this->V0 = cubiccrv(ctrl[0],ctrl[4],ctrl[8],ctrl[12]);
        this->V1 = cubiccrv(ctrl[3],ctrl[7],ctrl[11],ctrl[15]);
        this->U1 = cubiccrv(ctrl[12],ctrl[13],ctrl[14],ctrl[15]);

    }
    ~bicubicsrf() = default;

    vec3d at(double u, double v) override;
    vec3d atDerU(double u, double v) override;
    vec3d atDerV(double u, double v) override;
    vec3d atDerUU(double u, double v) override;
    vec3d atDerVV(double u, double v) override;
    vec3d atDerUV(double u, double v) override;

    curve & edgeU0() override;
    curve & edgeU1() override;
    curve & edgeV0() override;
    curve & edgeV1() override;
};


#endif //HEARTVALVEMODEL_BICUBICSRF_H
