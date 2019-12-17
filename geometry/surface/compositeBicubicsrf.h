//
// Created by s3179222 on 12/16/19.
//

#ifndef HEARTVALVEMODEL_COMPOSITEBICUBICSRF_H
#define HEARTVALVEMODEL_COMPOSITEBICUBICSRF_H

#include "bicubicsrf.h"
#include "../../measurements/measurements.h"


class compositeBicubicsrf {
private:
    bicubicsrf srf;
    compositeBicubicsrf *topLeft = nullptr, *topRight = nullptr, *bottomLeft = nullptr, *bottomRight = nullptr;
public:
    compositeBicubicsrf() = default;
    explicit compositeBicubicsrf(bicubicsrf &sur){
        srf = sur;
        std::cout << srf.hasValidControlNet() << std::endl;
        //plotSurface(srf);
        if(!srf.hasValidControlNet()){
            bicubicsrf tl, tr, bl, br;
            srf.subdivide(tl,tr,bl,br);
            topLeft = new compositeBicubicsrf(tl);
            topRight = new compositeBicubicsrf(tr);
            bottomLeft = new compositeBicubicsrf(bl);
            bottomRight = new compositeBicubicsrf(br);
        }
    }
    ~compositeBicubicsrf(){
        if(topLeft != nullptr){
            delete topLeft;
            delete topRight;
            delete bottomLeft;
            delete bottomRight;
        }
    };

    vec3d at(double u, double v);
    vec3d atDerU(double u, double v);
    vec3d atDerV(double u, double v);
    vec3d atDerUU(double u, double v);
    vec3d atDerVV(double u, double v);
    vec3d atDerUV(double u, double v);

    curve & edgeU0();
    curve & edgeU1();
    curve & edgeV0();
    curve & edgeV1();


    bool hasValidControlNet();
    bool closestPointInPatch(vec3d P);
};

#endif //HEARTVALVEMODEL_COMPOSITEBICUBICSRF_H
