//
// Created by lenferd on 27.10.16.
//

#include <iostream>
#include <omp.h>
#include "SparseMatrix.h"
using std::string;

vector<double> fillVect();
void fillMatrix(double** &matrix, string filename, int &size) {
    FILE *infile = fopen(filename.c_str(), "r");

    if (infile == NULL) {
        printf("File reading error. Try to relocate input file\n");
        exit(0);
    }

    // Scan size of matrix.
    fscanf(infile, "size=%d", &size);
    if (size == 0 || size < 0) {
        printf("Error, wrong size");
        exit(0);
    }

    matrix = new double*[size];
    for (int i = 0; i < size; ++i) {
        matrix[i] = new double[size];
        for (int j = 0; j < size; ++j) {
            fscanf(infile, "%lf", &matrix[i][j]);
        }
    }
    fclose(infile);
}



int main(int argc, char** argv) {

    // Timing variables
    double time_memoryS, time_memoryE;  // Time for allocate memory
    // double t0 = omp_get_wtime(); // Thanks Petrov

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
    time_memoryS = omp_get_wtime();

    double** vect = new double*[2];
    /*vect[0] = new double[1000000000];
    vect[1] = new double[1000000000];*/
    vect[0] = (double *) malloc((1000000000) * sizeof(double));
    vect[1] = (double *) malloc((1000000000) * sizeof(double));

    time_memoryE = omp_get_wtime();
    printf("Run time %.15lf\n", time_memoryE - time_memoryS);
    /*// Read file
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

    double expression = (sigma * dt) / (step * step);

    // Sparse Matrix fill

    SparseMatrix matrix;
    matrix.testEuler(nX+2, expression);

    // Calculating

    for (double j = 0; j < tFinal; j += dt) {
        matrix.multiplicateVector(vect[prevTime], vect[currTime], nX+2);
        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }



    // Output
    FILE *outfile = fopen("OUTPUT.txt", "w");

    for (int i = 1; i <= nX; i++) {
        fprintf(outfile, "%2.15le\n", vect[prevTime][i]);
    }*/
}
