//
// Created by lenferd on 04.04.17.
//

#include <iostream>
#include <omp.h>
#include <cmath>
#include <mpi.h>

const int ENABLE_PARALLEL = 0;

using std::string;

int main(int argc, char** argv) {
    /***
     * Initial
     */

    const int ROOT = 0;
    // Timing variables
    double time_S, time_E;  // Time for allocate memory

    string filename = "../initial/INPUT.txt";
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

//    printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinal %lf; dt %.10lf;\n",
//           xStart, xEnd, sigma, nX, tStart, tFinal, dt);

    double** vect = new double*[2];
    vect[0] = new double[nX + 2];
    vect[1] = new double[nX + 2];

    // Read file
    for (int i = 1; i <= nX; i++) {
        fscanf(infile, "%lf\n", &vect[0][i]);
    }
    fclose(infile);


    /***
     * Calculating
     */

    time_S = omp_get_wtime();

    step = (fabs(xStart) + fabs(xEnd)) / nX;      // calculate step

    prevTime = 0;
    currTime = 1;

    vect[0][0] = vect[0][1];
    vect[0][nX+1] = vect[0][nX];

    double timesize = (tFinal - tStart) / dt;

    printf("1/timesize:\t %.10lf\n", 1/timesize);
    printf("timesize:\t %.0f\n", timesize);
    printf("dt:\t\t %.2e\n", 1/timesize * (tFinal - tStart));

    double expr = (sigma * 1/timesize * (tFinal - tStart)) / (step * step);


    /***
     * MPI Stuff
     */
    int sizeP, rankP;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    // First scatter
    if (rankP == ROOT) {
        if (sizeP % 2 != 0 || sizeP == 1) {
            perror("Process counter must be divided by 2 or equal 1");
            exit(0);
        }
        // Calculate size of block for each process
        int sizeblock = nX / sizeP;

        // But we should send an two cell more each side
        

    }


    for (double j = 0; j < timesize; j += 1) {

        omp_set_num_threads(4);
        {
            #pragma omp parallel for if (ENABLE_PARALLEL)
            for (int i = 1; i <= nX; i++) {
                vect[currTime][i] = expr * (vect[prevTime][i + 1] - 2 * vect[prevTime][i] + vect[prevTime][i - 1])
                                    + vect[prevTime][i];
            }
        }
        // boundary conditions
        vect[currTime][0] = vect[currTime][1];
        vect[currTime][nX+1] = vect[currTime][nX];

        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;

    }

    time_E = omp_get_wtime();
    printf("Run time:\t %.15lf\n", time_E-time_S);

    FILE *outfile;

    string outfilename = "../result/Sergey/Sergey_MPI_SparseEuler1D.txt";
    outfile = fopen(outfilename.c_str(), "w");

    for (int i = 1; i <= nX; i++) {
        fprintf(outfile, "%2.15le\n", vect[prevTime][i]);
    }

    MPI_Finalize();
}
