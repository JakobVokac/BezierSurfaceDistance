//
// Created by s3179222 on 12/9/19.
//
/*
 * This class is the implementation of the Newton-Raphson method for the numerical optimizer.
 */
#ifndef HEARTVALVEMODEL_NEWTON_H
#define HEARTVALVEMODEL_NEWTON_H


#include "step.h"

class Newton : public step{
private:
    double sigma;
    double Gu, Gv;
    double Huu, Huv, Hvv;
public:
    explicit Newton(double sigma){
        this->sigma = sigma;
        Gu = 0;
        Gv = 0;
        Huu = 0;
        Huv = 0;
        Hvv = 0;
    }
    ~Newton() = default;
    OptState2D doStep(surface &sur, vec3d &P, OptState2D &loc) override;
    OptState1D doStep(curve &crv, vec3d &P, OptState1D &loc) override;
};


#endif //HEARTVALVEMODEL_NEWTON_H
