#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "model.h"
#include "optimizer/optimizer.h"
#include "optimizer/preprocessor/bisection.h"
#include "optimizer/step/Newton.h"
#include "optimizer/step/Geometric.h"
#include "optimizer/preprocessor/quadraticInterpolation.h"
#include "measurements/measurements.h"


double compareSurfaces(TopParametric top, bicubicsrf bez){
    double totalDist = 0;
    int iterations = 0;
    for (double i = 0; i <= 1; i += 0.01) {
        for (double j = 0; j <= 1; j += 0.01) {
            vec3d p1, p2;
            p1 = top.at(j,i);
            p2 = bez.at(j,i);
            totalDist += p1.dist(p2);
            iterations++;
        }
    }
    return totalDist/iterations;
}

double compareSurfaceEdges(TopParametric top, bicubicsrf bez){
    double totalDist = 0;
    int iterations = 0;
    for (double i = 0; i <= 1; i += 0.01) {
        vec3d p1, p2;
        p1 = top.at(0,i);
        p2 = bez.at(0,i);
        totalDist += p1.dist(p2);
        p1 = top.at(1,i);
        p2 = bez.at(1,i);
        totalDist += p1.dist(p2);
        p1 = top.at(i,0);
        p2 = bez.at(i,0);
        totalDist += p1.dist(p2);
        p1 = top.at(i,1);
        p2 = bez.at(i,1);
        totalDist += p1.dist(p2);
        iterations += 4;
    }
    return totalDist/iterations;
}
int main() {

    std::vector<double> inputPoints;
    ifstream inputFile("input.txt");        // Input file stream object

    // Check if exists and then open the file.
    if (inputFile.good()) {
        // Push items into a vector
        double current_number = 0;
        while (inputFile >> current_number){
            inputPoints.push_back(current_number);
        }

        // Close the file.
        inputFile.close();

        cout << endl;
    }else {
        cout << "Error reading input file!" << endl;
        cout << "Input file must have name \"input.txt\" and must consist of only numbers, no text!" << endl;
        cout << "The program will read the numbers in triples and each 3 numbers will be interpreted as one point." << endl;

        _exit(0);
    }

    if(inputPoints.size() % 3 != 0){
        cout << "Number of input numbers should be divisible by 3 (3 dimensions per point)!" << endl;
        cout << "Ignoring last few numbers." << endl;

        int size = inputPoints.size();
        int truncate = size % 3;
        for (int i = 0; i < truncate; ++i) {
            inputPoints.pop_back();
        }
    }

    Model model = Model(
            12.0,
            59.5/180 * M_PI,
            11.4,
            14.4,
            10,
            1.2,
            50.0/180*M_PI,
            7.2,
            16.8,
            3.5,
            1.35,
            -0.2,
            -0.2,
            0.01,
            1.0,
            6.5);

    Model p0 = Model::getPart(model,0),
          p1 = Model::getPart(model,1),
          p2 = Model::getPart(model,2),
          p3 = Model::getPart(model,3),
          p4 = Model::getPart(model,4),
          p5 = Model::getPart(model,5);

    preprocessor *bisect;
    preprocessor *quad;
    step *newton;
    step *geometric;
    optimizer *opt;
    quad = dynamic_cast<preprocessor*>(new quadraticInterpolation(8));
    bisect = dynamic_cast<preprocessor*>(new bisection(6));
    geometric = dynamic_cast<step*>(new Geometric());
    newton = dynamic_cast<step*>(new Newton(1));

    TopParametric top1 = p0.getTopParametric();
    TopParametric top2 = p3.getTopParametric();
    TopParametric top3 = p2.getTopParametric();
    TopParametric top4 = p5.getTopParametric();
    TopParametric top5 = p4.getTopParametric();
    TopParametric top6 = p1.getTopParametric();

    BottomParametric bot1 = p0.getBottomParametric();
    BottomParametric bot2 = p3.getBottomParametric();
    BottomParametric bot3 = p2.getBottomParametric();
    BottomParametric bot4 = p5.getBottomParametric();
    BottomParametric bot5 = p4.getBottomParametric();
    BottomParametric bot6 = p1.getBottomParametric();


    bicubicsrf bez = model.getTopBezier();
    surface *sur = dynamic_cast<surface*>(&top1);
    opt = new optimizer(
            *bisect,
            *quad,
            *newton,
            *newton,
            sur,
            {0,0,0},
            0.00000001,
            20,
            2
    );

    chrono::time_point<chrono::high_resolution_clock> t_start;
    chrono::time_point<chrono::high_resolution_clock> t_stop;
    chrono::microseconds t_duration;

    double sum = 0;

    t_start = chrono::high_resolution_clock::now();
    for (int i = 0; i < inputPoints.size()/3; i++) {

        vec3d P = {inputPoints[i*3 + 0],inputPoints[i*3 + 1],inputPoints[i*3 + 2]};
        if(P.x == 0){
            if(P.y == 0){
                opt->setSur(&top1);
                opt->setSur(&bot1);
            }
            else if(P.y > 0){
                opt->setSur(&top2);
                opt->setSur(&bot2);
            }
            else if(P.y < 0){
                opt->setSur(&top5);
                opt->setSur(&bot5);
            }
        }

        if(P.x > 0 && P.y > 0){
            if(P.y/P.x < sin(M_PI/6)/cos(M_PI/6)){
                opt->setSur(&top1);
                opt->setSur(&bot1);
            }else if(P.y/P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
                opt->setSur(&top1);
                opt->setSur(&bot1);
            }else{
                opt->setSur(&top2);
                opt->setSur(&bot2);
            }
        }else if(P.x < 0 && P.y < 0){
            if(P.y/-P.x < sin(M_PI/6)/cos(M_PI/6)){
                opt->setSur(&top3);
                opt->setSur(&bot3);
            }else if(P.y/-P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
                opt->setSur(&top3);
                opt->setSur(&bot3);
            }else{
                opt->setSur(&top2);
                opt->setSur(&bot2);
            }
        }else if(P.x < 0 && P.y < 0){
            if(-P.y/-P.x < sin(M_PI/6)/cos(M_PI/6)){
                opt->setSur(&top4);
                opt->setSur(&bot4);
            }else if(-P.y/-P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
                opt->setSur(&top4);
                opt->setSur(&bot4);
            }else{
                opt->setSur(&top5);
                opt->setSur(&bot5);
            }
        }else if(P.x > 0 && P.y < 0){
            if(-P.y/P.x < sin(M_PI/6)/cos(M_PI/6)){
                opt->setSur(&top6);
                opt->setSur(&bot6);
            }else if(-P.y/P.x < sin(2*M_PI/6)/cos(2*M_PI/6)){
                opt->setSur(&top6);
                opt->setSur(&bot6);
            }else{
                opt->setSur(&top5);
                opt->setSur(&bot5);
            }
        }
        sum += opt->optimizeForPoint(P).dist;
    }
    t_stop = chrono::high_resolution_clock::now();
    t_duration = chrono::duration_cast<chrono::microseconds>(t_stop - t_start);

    cout << "Computation time: " << t_duration.count() << " microseconds" << endl;

    cout << endl << "Total distance sum for N = " << inputPoints.size()/3 << " points: " << sum << ", average distance: " << sum/(inputPoints.size()/3) << endl;


    delete quad;
    delete geometric;
    delete opt;
    delete bisect;
    delete newton;
}



