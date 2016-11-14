#include <iostream>
#include <cmath>
#include <omp.h>

using std::string;

const int ENABLE_PARALLEL = 1;

int main(int argc, char** argv) {

    // Timing variables
    double time_S, time_E;  // Time for allocate memory

    string filename = "../../../../../initial/INPUT.txt";
    FILE *infile = fopen(filename.c_str(), "r");

    if (infile == NULL) {
        printf("File reading error. Try to relocate input file\n");
        exit(0);
    }

    double xStart = 0.0, xEnd = 0.0;
    double sigma = 0.0;

    int bc = 0; // Not use

    int nX = 1;
    double tStart = 0.0, tFinal = 0.0;
    double dt = 0.0;

    int i, j;
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

    double** vect = new double*[2];
    vect[0] = new double[nX + 2];
    vect[1] = new double[nX + 2];

    // Read file
    for (int i = 1; i <= nX; i++) {
        fscanf(infile, "%lf\n", &vect[0][i]);
    }
    fclose(infile);

    time_S = omp_get_wtime();

    step = (fabs(xStart) + fabs(xEnd)) / nX;      // calculate step
    prevTime = 0;
    currTime = 1;

    vect[0][0] = vect[0][1];
    vect[0][nX+1] = vect[0][nX];

    double* k1Vect = new double[nX + 2];
    double* k2Vect = new double[nX + 2];
    double* k3Vect = new double[nX + 2];
    double* k4Vect = new double[nX + 2];

    double expression = sigma / (step * step);

    printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinal %lf; dt %lf;\n",
           xStart, xEnd, sigma, nX, tStart, tFinal, dt);

    int timeStep = (int) ((tFinal - tStart) / dt);

    for (double j = 0; j < timeStep; j += 1) {
        omp_set_num_threads(3);
        // Fill k1 vect
        #pragma omp parallel for if (ENABLE_PARALLEL)
        for (int i = 1; i <= nX; i++) {
            k1Vect[i] = (vect[prevTime][i + 1] - 2.0 * vect[prevTime][i] + vect[prevTime][i - 1]) * expression;
        }
        //bounders k1
        k1Vect[0] = k1Vect[1];
        k1Vect[nX + 1] = k1Vect[nX];

        // Fill k2.0 vect
        #pragma omp parallel for if (ENABLE_PARALLEL)
        for (int i = 1; i <= nX; i++) {
            k2Vect[i] = expression * ((vect[prevTime][i + 1]) + (k1Vect[i + 1] * dt / 2.0)
                        - (2.0 * vect[prevTime][i] + k1Vect[i] * dt)
                        + vect[prevTime][i - 1] + k1Vect[i - 1] * dt / 2.0);
        }
        // bounders
        k2Vect[0] = k2Vect[1];
        k2Vect[nX + 1] = k2Vect[nX];

        // Fill k3 vect
        #pragma omp parallel for if (ENABLE_PARALLEL)
        for (int i = 1; i <= nX; i++) {
            k3Vect[i] = expression * ((vect[prevTime][i + 1]) + (k2Vect[i + 1] * dt / 2.0)
                        - (2.0 * vect[prevTime][i] + k2Vect[i] * dt)
                        + vect[prevTime][i - 1] + k2Vect[i - 1] * dt / 2.0);
        }
        // bounders
        k3Vect[0] = k3Vect[1];
        k3Vect[nX + 1] = k3Vect[nX];

        // Fill k4 vect
        #pragma omp parallel for if (ENABLE_PARALLEL)
        for (int i = 1; i <= nX; i++) {
            k4Vect[i] = expression * ((vect[prevTime][i + 1]) + (k3Vect[i + 1] * dt)
                        - (2.0 * vect[prevTime][i] + k3Vect[i] * 2.0 * dt)
                        + vect[prevTime][i - 1] + k3Vect[i - 1] * dt);
        }
        // bounders
        k4Vect[0] = k4Vect[1];
        k4Vect[nX + 1] = k4Vect[nX];

        // FIll result vector
        #pragma omp parallel for if (ENABLE_PARALLEL)
        for (int i = 1; i <= nX; i++) {
            vect[currTime][i] = vect[prevTime][i] + ((dt / 6) * (k1Vect[i] + 2.0*k2Vect[i] + 2.0*k3Vect[i] + k4Vect[i]));

        }

        // boundary conditions
        vect[currTime][0] = vect[currTime][1];
        vect[currTime][nX+1] = vect[currTime][nX];


        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }

    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);

    FILE *outfile = fopen("OUTPUT_Runge.txt", "w");

    for (int i = 1; i <= nX; i++) {
        fprintf(outfile, "%2.15le\n", vect[prevTime][i]);
    }

}
