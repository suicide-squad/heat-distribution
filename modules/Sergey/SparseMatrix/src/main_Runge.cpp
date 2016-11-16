//
// Created by lenferd on 10.11.16.
//

#include <iostream>
#include <omp.h>
#include <cmath>
#include "SparseMatrix.h"
using std::string;


int main(int argc, char** argv) {


    // Timing variables
    double time_S, time_E;  // Time for allocate memory

    // File open

    string filename = "../../../../../initial/INPUT.txt";
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

    double* v_k1 = new double[nX + 2];
    double* v_k2 = new double[nX + 2];
    double* v_k3 = new double[nX + 2];
    double* v_k4 = new double[nX + 2];


    SparseMatrix sm_k1;
    spMatrixInit(sm_k1, nX * 3 + 2, nX + 2);
    double expression1 = sigma / (step * step);
    double expression2 = -2.0 * expression1;
    fillMatrix2Expr(sm_k1, nX + 2, expression1, expression2);

    // KURWA
    SparseMatrix sm_k2;
    spMatrixInit(sm_k2, nX * 3 + 2, nX + 2);
    double k2expr1 = dt * expression1 * 0.5;
    double k2expr2 = 1 - 2.0 * k2expr1;
    fillMatrix2Expr(sm_k2, nX+2, k2expr1, k2expr2);

    SparseMatrix sm_k3;
    spMatrixInit(sm_k3, nX * 3 + 2, nX + 2);
    double k3expr1 = k2expr1;
    double k3expr2 = k2expr2;
    fillMatrix2Expr(sm_k3, nX+2, k3expr1, k3expr2);

    SparseMatrix sm_k4;
    spMatrixInit(sm_k4, nX * 3 + 2, nX + 2);
    double k4expr1 = dt * expression1;
    double k4expr2 = 1 - 2.0 * k4expr1;
    fillMatrix2Expr(sm_k4, nX+2, k4expr1, k4expr2);


    // Calculating

    time_S = omp_get_wtime();
    double expressionResult = dt / 6;
    for (double j = 0; j < tFinal; j += dt) {
        multiplicateVector(sm_k1, vect[prevTime], v_k1, nX+2);
        multiplicateVector(sm_k2, v_k1, v_k2, nX+2);
        multiplicateVector(sm_k3, v_k2, v_k3, nX+2);
        multiplicateVector(sm_k4, v_k3, v_k4, nX+2);

        // Fill result vector
        for (int i = 1; i <= nX; i++) {
            vect[currTime][i] = vect[prevTime][i] +
                    expressionResult * (v_k1[i] + 2.0 * v_k2[i] + 2.0 * v_k3[i] + v_k4[i]);
        }

        vect[currTime][0] = vect[currTime][1];
        vect[currTime][nX+1] = vect[currTime][nX];

        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);


    // Output
    FILE *outfile = fopen("OUTPUT.txt", "w");

    for (int i = 1; i <= nX; i++) {
        fprintf(outfile, "%2.15le\n", vect[prevTime][i]); }
}

