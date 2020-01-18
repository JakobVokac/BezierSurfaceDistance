//
// Created by s3179222 on 12/9/19.
//
/*
 * The structs OptState1D and OptState2D represent states in the optimization, which include the current location on
 * the surface (u,v or t) and the current distance to the point were optimizing for (dist). OptState2D is used for the
 * inside of a surface. Optstate1D is used for edges.
 */
#ifndef HEARTVALVEMODEL_STATES_H
#define HEARTVALVEMODEL_STATES_H


struct OptState2D{
    double u,v;
    double dist;
};

struct OptState1D{
    double t;
    double dist;
};

#endif //HEARTVALVEMODEL_STATES_H
