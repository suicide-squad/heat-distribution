//
// Created by lenferd on 28.11.16.
//

#include "SparseMatrix.h"

#include <iostream>
#include <cmath>
#include <omp.h>

using std::string;

double norm2Vect(double *&vect1, double *&vect2, int size) {
    double norm = 0;
    norm = fabs(vect1[0] - vect2[0]);
    for (int h = 0; h < size; h++) {
        if (fabs(vect1[h] - vect2[h]) > norm)
            norm = fabs(vect1[h] - vect2[h]);
    }
    return norm;
}

int main() {


    // Timing variables
    double time_S, time_E;  // Time for allocate memory

    // File open

    string filename = "../../../../initial/INPUT.txt";
    FILE *infile = fopen(filename.c_str(), "r");

    if (infile == NULL) {
        printf("File reading error. Try to relocate input file\n");
        exit(0);
    }

    // Init variables
    double xStart = 0.0, xEnd = 0.0;
    double sigma = 0.0;

    int bc = 0; // Not use

    int nX = 1;
    double tStart = 0.0, tFinal = 0.0;
    double dt = 0.0;

    int i, j;
    int sizeTime = 0;
    int currTime, prevTime;

    double step = 0.0;

    double eps = 1e-7;
    //  File reading

    fscanf(infile, "XSTART=%lf\n", &xStart);    // start coordinate
    fscanf(infile, "XEND=%lf\n", &xEnd);        // end coordinate
    fscanf(infile, "SIGMA=%lf\n", &sigma);      // coef of heat conduction
    fscanf(infile, "NX=%d\n", &nX);             // count of initial elements?
    fscanf(infile, "TSTART=%lf\n", &tStart);    // start time
    fscanf(infile, "TFINISH=%lf\n", &tFinal);   // finish time
    fscanf(infile, "dt=%lf\n", &dt);            // delta of time difference
    fscanf(infile, "BC=%d\n", &bc);         // Not using right now

    printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinal %lf; dt %lf;\n",
           xStart, xEnd, sigma, nX, tStart, tFinal, dt);

    //  Memory allocation

    double** vect = new double*[2];
    vect[0] = new double[nX + 2];
    vect[1] = new double[nX + 2];

    // Read file
    for (int i = 1; i <= nX; i++) {
        fscanf(infile, "%lf\n", &vect[0][i]);
    }
    fclose(infile);

    //  Prev val calculating
    step = (fabs(xStart) + fabs(xEnd)) / nX;      // calculate step

    prevTime = 0;
    currTime = 1;

    vect[0][0] = vect[0][1];
    vect[0][nX+1] = vect[0][nX];
    // Sparse Matrix fill

    SparseMatrix spMatrix;
    spMatrixInit(spMatrix, nX*3 + 2, nX + 2);
//    double expr1 = -sigma * dt / (step * step);
//    double expr2 = 1 + 2 * (sigma * dt / (step * step));
    double dopexprUp = (sigma * dt) / (step * step);
    double dopexprDw = 1 + 2 * dopexprUp;
    double expr1 = dopexprUp / dopexprDw;

    double expr3 = 1 / dopexprDw;
//  double expr2 = 1 / (1 + 2 * (sigma * dt / (step * step)));
    double expr2 = 0;
    fillMatrix2Expr(spMatrix, nX + 2, expr1, expr2);

    std::cout << expr1 << std::endl;
    std::cout << expr2 << std::endl;
    std::cout << expr3 << std::endl;

    // x1 = ax0 + b
    prevTime = 0;
    currTime = 1;

    time_S = omp_get_wtime();

    int iCOUNT = 0;

    // prevTime is b vect now
    double* const_vect = new double[nX + 2];

    double* vis_vect = new double[nX + 2];


    for (double j = 0; j < tFinal; j += dt) {
        for (int i = 0; i < nX + 2; i++) {
//            vis_vect[i] = vect[prevTime][i];
            const_vect[i] = vect[prevTime][i] * expr3 ;
//            vect[prevTime][i] = 1;
//            std::cout << const_vect[i] << std::endl;
        }

        int iterationCounter = 0;
        do {
            multiplicateVector(spMatrix, vect[prevTime], vect[currTime], nX + 2);
            for (int i = 0; i < nX + 2; i++) {
                vect[currTime][i] += const_vect[i];
            }

            prevTime = (prevTime + 1) % 2;
            currTime = (currTime + 1) % 2;

//            if (iCOUNT % 1 == 0)
//               std::cout << norm2Vect(vis_vect, vect[currTime], nX + 2) << "======"  << j << std::endl;
//            iCOUNT++;
//            std::cout << "====" ><< j << std::endl;
            iterationCounter++;
        } while (norm2Vect(vect[prevTime], vect[currTime], nX + 2) > eps);  // vect[prevTime] = k + 1 now
        if (iCOUNT % 100 == 0)
            printf("time %lf, iterations %d\n", j, iterationCounter);
        iCOUNT++;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);


    // Output
    FILE *outfile = fopen("OUTPUT.txt", "w");

    for (int i = 1; i <= nX; i++) {
        fprintf(outfile, "%2.15le\n", vect[prevTime][i]); }

}

