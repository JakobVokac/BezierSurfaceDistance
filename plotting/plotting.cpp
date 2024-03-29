//
// Created by s3179222 on 12/3/19.
//
/*
 * This class contains functions that use the matplotlib extension to draw the model's surfaces and distance functions.
 */
#include "plotting.h"

void DrawDist(Model &model, const vec3d &P) {
    vector<vector<double>> u, v, dist;
    double min_dist = 100000;
    for (double i = -0.5; i < 1.5; i += 0.15) {
        vector<double> uRow, vRow, distRow;
        for (double j = -0.5; j < 1.5; j += 0.15) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(model.distanceToTopPoint(j, i, P));
            if (min_dist > model.distanceToTopPoint(j, i, P)) {
                min_dist = model.distanceToTopPoint(j, i, P);
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    vector<double> x2, y2, z2, u2, v2, w2;
    for (double i = -0.5; i < 1.5; i += 0.15) {
        for (double j = -0.5; j < 1.5; j += 0.15) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.completeDistanceTopDerU(j, i, P);
            tV = model.completeDistanceTopDerV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);
    plt::show();
}

double constraintFunc(double u, double v){
    return pow((u < 0 ? -u : 0),3.0)
           + pow((u > 1 ? u-1 : 0),3.0)
           + pow((v < 0 ? -v : 0),3.0)
           + pow((v > 1 ? v-1 : 0),3.0);
}

//void plotSurface(vec3d P, Model m, double ul = -.2, double vl = -.2, double ur = 1.2, double vr = 1.2, double dU = 0.1, double dV = 0.1, int withConst = 0){
//    vector<vector<double>> u,v,dist;
//    double min_dist = 100000;
//    for (double i = ul-dU; i < ur+dU; i += (ur+dU - (ul-dU))/10) {
//        vector<double> uRow, vRow, distRow;
//        for (double j = vl-dV; j < vr+dV; j += (vr+dV - (vl-dV))/10) {
//            uRow.push_back(j);
//            vRow.push_back(i);
//            double c = 0;
//            if(withConst)
//                c = 100*constraintFunc(j,i);
//            distRow.push_back(m.distanceToTopPoint(j,i,P) + c);
//            if(min_dist > m.distanceToTopPoint(j,i,P) + c)
//                min_dist = m.distanceToTopPoint(j,i,P) + c;
//        }
//        u.push_back(uRow);
//        v.push_back(vRow);
//        dist.push_back(distRow);
//    }
//    plt::plot_surface(u,v,dist);
//}

void plotSurface(surface &s, double eps){
    if(eps >= 1){
        cout << "surface plotting: eps needs to be a fraction of 1" << endl;
        return;
    }
    vector<vector<double>> x,y,z;
    for (double i = 0; i <= 1; i+=eps) {
        vector<double> xr,yr,zr;
        for (double j = 0; j <= 1; j+=eps) {
            vec3d p = s.at(j,i);
            xr.push_back(p.getx());
            yr.push_back(p.gety());
            zr.push_back(p.getz());
        }
        x.push_back(xr);
        y.push_back(yr);
        z.push_back(zr);
    }
//    plt::plot_surface_with_line_and_axes_set(x,y,z,{-1,5,-1,5,-1,5},{},{},{},{},{},{},0);
    plt::plot_surface(x,y,z);
    plt::show();
}

void drawFillTopWithLine(Model m, vector<double> line, vector<double> axes){
    assert(line.size() == 6);
    vector<vector<double>> xs, ys, zs;
    for (double u = 0; u <= 1.0; u+=0.1) {
        vector<double> xr, yr, zr;
        for (double v = 0; v <= 1.0; v+=0.1) {
            vec3d p = m.fillTop(u,v);
            xr.push_back(p.getx());
            yr.push_back(p.gety());
            zr.push_back(p.getz());
        }
        xs.push_back(xr);
        ys.push_back(yr);
        zs.push_back(zr);
    }
    double len = sqrt((line[3])*(line[3]) + (line[4])*(line[4]) + (line[5])*(line[5]));
    plt::plot_surface_with_line_and_axes_set(xs,ys,zs,axes,{line[0]},{line[1]},{line[2]},{line[3]},{line[4]},{line[5]},len);
}

void plotEdgeDistSq(curve &c, vec3d P){
    vector<double> xs, ys, yds, ydds;
    for (double i = 0; i <= 1; i += 0.01) {
        xs.push_back(i);
        ys.push_back(c.sqDistTo(i,P));
        yds.push_back(c.sqDistToDer1(i,P));
        ydds.push_back(c.sqDistToDer2(i,P));

    }

    plt::plot(xs,ys);
    plt::show();

    plt::plot(xs,yds);
    plt::show();

    plt::plot(xs,ydds);

    plt::show();
}

void drawPart(Model &model, vector<vector<double>> &x, vector<vector<double>> &y, vector<vector<double>> &z, int top, int bottom) {

    if(top) {
        for (double i = 1; i >= 0; i -= 0.1) {
            vector<double> x_row, y_row, z_row;
            for (double j = 0; j <= 1; j += 0.1) {
                vec3d p = model.fillTop(j, i);
                x_row.push_back(p.getx());
                y_row.push_back(p.gety());
                z_row.push_back(p.getz());
            }
            x.push_back(x_row);
            y.push_back(y_row);
            z.push_back(z_row);
        }
    }

    if(bottom) {
        for (double i = 0.1; i <= 1; i += 0.1) {
            vector<double> x_row, y_row, z_row;
            for (double j = 0; j <= 1; j += 0.1) {
                vec3d p = model.fillBottom(j, i);
                x_row.push_back(p.getx());
                y_row.push_back(p.gety());
                z_row.push_back(p.getz());
            }
            x.push_back(x_row);
            y.push_back(y_row);
            z.push_back(z_row);
        }
    }
}

void DrawDistWithDirection(Model &model, const vec3d &P) {
    vector<vector<double>> u, v, dist;

    double min_dist = 100000;
    double max_dist = -100000;
    for (double i = 0; i <= 1.05; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j <= 1.05; j += 0.05) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(model.squaredTopDist(j, i, P));
            if (min_dist > model.squaredTopDist(j, i, P)) {
                min_dist = model.squaredTopDist(j, i, P);
            }
            if (max_dist < model.squaredTopDist(j, i, P)) {
                max_dist = model.squaredTopDist(j, i, P);
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    vector<double> x2, y2, z2, u2, v2, w2;
    for (double i = 0; i <= 1.05; i += 0.05) {
        for (double j = 0; j <= 1.05; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.squaredTopDistDerU(j, i, P);
            tV = model.squaredTopDistDerV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);

    plt::show();

    u.clear(), v.clear(), dist.clear();

    min_dist = 100000;
    max_dist = -100000;
    for (double i = 0; i <= 1.05; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j <= 1.05; j += 0.05) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(model.squaredTopDistDerU(j, i, P));
            if (min_dist > model.squaredTopDistDerU(j, i, P)) {
                min_dist = model.squaredTopDistDerU(j, i, P);
            }
            if (max_dist < model.squaredTopDistDerU(j, i, P)) {
                max_dist = model.squaredTopDistDerU(j, i, P);
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    x2.clear(), y2.clear(), z2.clear(), u2.clear(), v2.clear(), w2.clear();
    for (double i = 0; i <= 1.05; i += 0.05) {
        for (double j = 0; j <= 1.05; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.squaredTopDistDerUU(j, i, P);
            tV = model.squaredTopDistDerUV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);

    plt::show();

    u.clear(), v.clear(), dist.clear();

    min_dist = 100000;
    max_dist = -100000;
    for (double i = 0; i <= 1.05; i += 0.05) {
        vector<double> uRow, vRow, distRow;
        for (double j = 0; j <= 1.05; j += 0.05) {
            uRow.push_back(j);
            vRow.push_back(i);
            distRow.push_back(model.squaredTopDistDerV(j, i, P));
            if (min_dist > model.squaredTopDistDerV(j, i, P)) {
                min_dist = model.squaredTopDistDerV(j, i, P);
            }
            if (max_dist < model.squaredTopDistDerV(j, i, P)) {
                max_dist = model.squaredTopDistDerV(j, i, P);
            }
        }
        u.push_back(uRow);
        v.push_back(vRow);
        dist.push_back(distRow);
    }
    x2.clear(), y2.clear(), z2.clear(), u2.clear(), v2.clear(), w2.clear();
    for (double i = 0; i <= 1.05; i += 0.05) {
        for (double j = 0; j <= 1.05; j += 0.05) {
            x2.push_back(j);
            y2.push_back(i);
            z2.push_back(min_dist);
            double tU, tV;
            tU = model.squaredTopDistDerUV(j, i, P);
            tV = model.squaredTopDistDerVV(j, i, P);
            u2.push_back(tU);
            v2.push_back(tV);
            w2.push_back(0);
        }
    }
    plt::plot_surface_with_vector_field(u, v, dist, x2, y2, z2, u2, v2, w2, 0.1);

    plt::show();

}