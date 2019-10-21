//
// Created by s3179222 on 10/12/19.
//

#ifndef HEARTVALVEMODEL_PARAMETERS_H
#define HEARTVALVEMODEL_PARAMETERS_H

#include <cmath>

double root_radius = 12.0;
double leaflet_angle = 59.5/180 * M_PI;
double leaflet_height = 11.4;
double leaflet_radius_out = 14.4;
double leaflet_radius_in = 3.5;
double leaflet_power = 1.35;

double bending_angle = 50.0/180 * M_PI;
double bending_height = 0.5 * leaflet_height;
double bending_radius_out = 16.8;
double bending_radius_in = 15.5;
double bending_power = 1.35;

double s_sin = bending_height / leaflet_height;
double s_sym = bending_height / leaflet_height;
double xQ1_sym = -0.2;


#endif //HEARTVALVEMODEL_PARAMETERS_H
