//
// Created by s3179222 on 12/16/19.
//

#ifndef HEARTVALVEMODEL_BEZIERSURFACE_H
#define HEARTVALVEMODEL_BEZIERSURFACE_H


#include "surface.h"

class BezierSurface : public surface {
public:
    virtual ~BezierSurface() = 0;
    virtual bool hasValidControlNet() = 0;
    virtual void controlNetSubdivide() = 0;
    virtual bool closestPointInPatch(vec3d P) = 0;
    virtual void subdivide(BezierSurface *tl,  BezierSurface *tr, BezierSurface *bl, BezierSurface*br) = 0;
};


#endif //HEARTVALVEMODEL_BEZIERSURFACE_H
