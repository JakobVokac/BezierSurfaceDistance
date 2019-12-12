//
// Created by s3179222 on 12/3/19.
//
#include "vectorMath.h"

vector<double> add(vector<double> P, vector<double> Q){
    assert(P.size() == Q.size());
    for (int i = 0; i < P.size(); ++i) {
        P[i] += Q[i];
    }
    return P;
}
vector<double> sub(vector<double> P, vector<double> Q){
    assert(P.size() == Q.size());
    for (int i = 0; i < P.size(); ++i) {
        P[i] -= Q[i];
    }
    return P;
}
vector<double> mul(vector<double> P, double q){
    for (int i = 0; i < P.size(); ++i) {
        P[i] *= q;
    }
    return P;
}

vector<double> div(vector<double> P, double q){
    for (int i = 0; i < P.size(); ++i) {
        P[i] /= q;
    }
    return P;
}

double det(vector<double> P, vector<double> Q){
    return P[0]*Q[1] - P[1]*Q[0];
}

double sign2D(vector<double> P, vector<double> Q) {
    if(P[0]*Q[0] + P[1]*Q[1] > 0){
        return 1.0;
    }else{
        return -1.0;
    }
}

double magnitude(vector<double> v){
    double sum = 0;
    for (int i = 0; i < v.size(); ++i) {
        sum += v[i]*v[i];
    }
    return sqrt(sum);
}


double dot(vector<double> a, vector<double> b){
    double sum = 0;
    assert(a.size() == b.size());
    for (int i = 0; i < a.size(); ++i) {
        sum += a[i]*b[i];
    }
    return sum;
}

vector<double> cross(vector<double> a, vector<double> b){
    assert(a.size() == 3);
    assert(b.size() == 3);
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}

double curvature3D(vector<double> c1, vector<double> c2){
    double mag = magnitude(c1);
    return magnitude(cross(c1,c2))/(mag*mag*mag);
}
vector<double> circleCenterDir3D(vector<double> c1, vector<double> c2) {
    vector<double> c = cross(c1,c2);
    c = cross(c,c1);
    return div(c, magnitude(c));
}

vector<double> solve3Dlinear(vector<double> x, vector<double> y, vector<double> z){
    assert(x.size() == 4);
    assert(y.size() == 4);
    assert(z.size() == 4);

    if(x[2] == 0){
        if(z[2] != 0){
            vector<double> temp = x;
            x = z;
            z = temp;
        }else{
            vector<double> temp = x;
            x = y;
            y = temp;
        }
    }
    if(y[1] == 0){
        if(z[1] != 0){
            vector<double> temp = y;
            y = z;
            z = temp;
        }else{
            vector<double> temp = y;
            y = x;
            x = temp;
        }
    }
    if(z[0] == 0){
        if(y[0] != 0){
            vector<double> temp = y;
            y = z;
            z = temp;
        }else{
            vector<double> temp = z;
            z = x;
            x = temp;
        }
    }
    z = sub(z,mul(x,z[2]/x[2]));
    y = sub(y,mul(x,y[2]/x[2]));
    z = sub(z,mul(y,z[1]/y[1]));
    y = sub(y,mul(z,y[0]/z[0]));
    x = sub(x,mul(z,x[0]/z[0]));
    x = sub(x,mul(y,x[1]/y[1]));

    return {z[3]/z[0],y[3]/y[1],x[3]/x[2]};

}

double sign3D(vector<double> P, vector<double> Q) {
    if(dot(P,Q) > 0){
        return 1.0;
    }else{
        return -1.0;
    }
}


vector<double> circleCenterDir2D(vector<double> c1, vector<double> c2){
    double sign = 0;
    if(det(c1,c2) > 0){
        sign = 1;
    }else{
        sign = -1;
    }
    return div({c1[1] * (-sign),c1[0] * sign},magnitude(c1));
}

double curvature2D(vector<double> c1, vector<double> c2){
    double mag = magnitude(c2);
    return det(c1,c2)/(mag*mag*mag);
}