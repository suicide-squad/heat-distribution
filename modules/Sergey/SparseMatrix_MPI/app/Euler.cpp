//
// Created by lenferd on 04.04.17.
//

#include <iostream>
#include <cmath>
#include <mpi.h>
#include "SparseMatrix.h"

using std::string;

int main(int argc, char** argv) {
    /***
     * Initial
     */

    double startTime, endTime;
    const int ROOT = 0;
    const int ADD_CELL = 2;  // Additional cell (in each size)

    /***
    * MPI Stuff
    */
    int sizeP = 0, rankP = 0;
    MPI_Status *status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    printf("%d %d\n", rankP, sizeP);


    double xStart = 0.0, xEnd = 0.0;
    double sigma = 0.0;

    int bc = 0; // Not use

    int nX;
    double tStart = 0.0, tFinal = 0.0;
    double dt = 0.0;

    int i, j;
    int sizeTime = 0;

    double timeStep = 0.0;

    double **vect;

    if (rankP == ROOT) {
        string filename = "../initial/INPUT.txt";
        FILE *infile = fopen(filename.c_str(), "r");

        if (infile == NULL) {
            printf("File reading error. Try to relocate input file\n");
            exit(0);
        }

        //  File reading
        fscanf(infile, "XSTART=%lf\n", &xStart);    // start coordinate
        fscanf(infile, "XEND=%lf\n", &xEnd);        // end coordinate
        fscanf(infile, "SIGMA=%lf\n", &sigma);      // coef of heat conduction
        fscanf(infile, "NX=%d\n", &nX);             // count of initial elements?
        fscanf(infile, "TSTART=%lf\n", &tStart);    // start time
        fscanf(infile, "TFINISH=%lf\n", &tFinal);   // finish time
        fscanf(infile, "dt=%lf\n", &dt);            // delta of time difference
        fscanf(infile, "BC=%d\n", &bc);         // Not using right now

        printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinal %lf; dt %.10lf;\n",
           xStart, xEnd, sigma, nX, tStart, tFinal, dt);

        vect = new double *[2];
        // if we use add cell, we need add more that for our boundaries conditionals
        int new_size = nX + ADD_CELL * 2;
        vect[0] = new double[new_size];
        vect[1] = new double[new_size];

        // Fill zero
        for (int k = 0; k < new_size; ++k) {
            vect[0][k] = 0;
            vect[1][k] = 0;
        }

        // ..[....]..
        // if we add some additional cell, we should change start reading position
        for (int i = ADD_CELL; i <= nX + ADD_CELL; i++) {
            fscanf(infile, "%lf\n", &vect[0][i]);
        }
        fclose(infile);
    }

    // now we should scatter data for all proc
    MPI_Bcast(&sizeTime, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&nX, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    // CAREFUL. Now all process know only about:
    //  sizeTime
    //  nX          size of original vector
    //  dt

    int block_size = nX / sizeP;    // block size wo boundaries addition
    int block_size_add = block_size + ADD_CELL * 2;

    // Generate data of postions
    int *displs = new int[sizeP];
    int *sendcounts = new int[sizeP];

    for (int l = 0; l < sizeP; ++l) {
        displs[l] = block_size * l;
        sendcounts[l] = block_size_add;
    }

    // Init recv
    double **procVect = new double*[2];
    procVect[0] = new double[block_size_add];
    procVect[1] = new double[block_size_add];


    // Scatter vector
    MPI_Scatterv(vect, sendcounts, displs, MPI_DOUBLE, procVect[0], block_size_add, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    for (int i = 0; i < block_size_add; ++i) {
//        printf("%d %lf\n", rankP, procVect[0][i]);
    }

//
//
////    MPI_Scatterv(vect, block_size_add, MPI_DOUBLE, Ur + 1, block, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
//
//
//    /***
//     * Calculating
//     */
//
//    timeStep = (fabs(xStart) + fabs(xEnd)) / nX;      // calculate step
//
//    int prevTime = 0;
//    int currTime = 1;
//
//    vect[0][0] = vect[0][1];
//    vect[0][nX+1] = vect[0][nX];
//
//    double timesize = (tFinal - tStart) / dt;
//
////    printf("1/timesize:\t %.10lf\n", 1/timesize);
////    printf("timesize:\t %.0f\n", timesize);
////    printf("dt:\t\t %.2e\n", 1/timesize * (tFinal - tStart));
////    printf("%s", argv[1]);
////    double expr = (sigma * 1/timesize * (tFinal - tStart)) / (step * step);
//
//
////    printf("%d", rankP);
//
//    // First scatter
//    // Calculate size of block for each process
//    int sizeblock = nX / sizeP;
//
//    // But we should send an two cell more each side
//    int sizeVect = sizeblock + 4;
//
//    double **procVect = new double*[2];
//    procVect[0] = new double[sizeVect];
//    procVect[1] = new double[sizeVect];
//
//    // Some pointers
//    double* sendData;
//
////    printf("СУКА 1 %d\n", sizeVect);
////    if (rankP == ROOT) {
////        if (sizeP % 2 != 0 && sizeP != 1) {
////            perror("Process counter must be divided by 2 or equal 1");
////            exit(0);
////        } else if (sizeP == 1){
////            exit(0);
////        }
////        // Preparation vectors
////        sendData = new double[sizeP * sizeVect];
////        for (int i = 0; i < sizeP; ++i) {
////            // Another fill for boundaries
////            if (i == 0) {
////                sendData[i * sizeVect + 0] = 0;
////                sendData[i * sizeVect + 1] = 0;
////                for (int k = 2; k < sizeVect; ++k) {
////                    sendData[i * sizeVect + k] = vect[prevTime][k - 2];   // Don't have first
////                }
////            } else if ( i + 1 == sizeP) {
////                for (int k = 0; k < sizeVect - 2; ++k) {
////                    sendData[i * sizeVect + k] = vect[prevTime][i * sizeblock + k];
////                }
////                sendData[i * sizeVect + sizeVect - 1] = 0;
////                sendData[i * sizeVect + sizeVect - 2] = 0;
////            } else {
////                for (int k = 0; k < nX; ++k) {
////                    sendData[i * sizeVect + k] = vect[prevTime][i * sizeblock + k];
////                }
////            }
////        }
////    }
////
////    int *displs = new int[sizeP];
////    int *scounts = new int[sizeP];
////    for (i=0; i< sizeP; ++i) {
////        displs[i] = i * sizeVect;
////        scounts[i] = sizeVect;
////    }
////    // After preparation, send vectors
////    MPI_Scatterv(sendData, scounts, displs, MPI_DOUBLE, procVect[0], sizeVect, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
////
////
////    /**
////     * Sparse matrix
////     */
////    printf("СУКА 2 %d\n", rankP);
////    double expression = (sigma * dt) / (timeStep * timeStep);
////    double expression2 = 1.0 - 2.0 * expression;
////
////    // Sparse Matrix fill
////
////    printf("СУКА 11 %d\n", rankP);
////    SparseMatrix matrix;
////    printf("СУКА 11 %d\n", rankP);
////
////    spMatrixInit(matrix, sizeVect * 3, sizeVect);
////    printf("СУКА 11 %d\n", sizeVect);
////    fillMatrix2ExprWithoutBoundaries(matrix, sizeVect, expression, expression2);
////    printf("СУКА 14 %d\n", rankP);
////
////    startTime = MPI_Wtime();
////    int fuckCounter = 0;
////    printf("СУКА 12 %d\n", rankP);
////
////    double *vectRightSend = new double[2];
////    double *vectRightRecv = new double[2];
////    double *vectLeftSend = new double[2];
////    double *vectLeftRecv = new double[2];
////
////    printf("СУКА 3 %d\n", rankP);
////    for (double l = 0; l < tFinal; l += dt) {
////        fuckCounter++;
////        multiplicateVector(matrix, procVect[prevTime], procVect[currTime], sizeVect);
////
////        printf("СУКА 4 %d %d\n", rankP, l);
////        if (fuckCounter == 3) {
////            if (rankP == ROOT) {
////                vectRightSend[0] = procVect[currTime][sizeVect - 4];
////                vectRightSend[1] = procVect[currTime][sizeVect - 3];
////                MPI_Sendrecv(vectRightSend, 2, MPI_DOUBLE, rankP + 1, 0,
////                             vectRightRecv, 2, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, status);
////                procVect[currTime][sizeVect - 4] = vectRightRecv[0];
////                procVect[currTime][sizeVect - 3] = vectRightRecv[1];
////            } else if (rankP == sizeP - 1) {
////                vectLeftSend[0] = procVect[currTime][2];
////                vectLeftSend[1] = procVect[currTime][3];
////                MPI_Sendrecv(vectLeftSend, 2, MPI_DOUBLE, rankP - 1, 0,
////                             vectLeftRecv, 2, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, status);
////                procVect[currTime][2] = vectLeftRecv[0];
////                procVect[currTime][3] = vectLeftRecv[1];
////            } else {
////                vectRightSend[0] = procVect[currTime][sizeVect - 4];
////                vectRightSend[1] = procVect[currTime][sizeVect - 3];
////                vectLeftSend[0] = procVect[currTime][2];
////                vectLeftSend[1] = procVect[currTime][3];
////                MPI_Sendrecv(vectRightSend, 2, MPI_DOUBLE, rankP + 1, 0,
////                             vectRightRecv, 2, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, status);
////                MPI_Sendrecv(vectLeftSend, 2, MPI_DOUBLE, rankP - 1, 0,
////                             vectLeftRecv, 2, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, status);
////                procVect[currTime][sizeVect - 4] = vectRightRecv[0];
////                procVect[currTime][sizeVect - 3] = vectRightRecv[1];
////                procVect[currTime][2] = vectLeftRecv[0];
////                procVect[currTime][3] = vectLeftRecv[1];
////            }
////
////        }
////
////        // boundaries
////        if (rankP == ROOT) {
////            procVect[currTime][1] = procVect[currTime][2];
////        } else if (rankP == sizeP - 1) {
////            procVect[currTime][sizeVect - 3] = procVect[currTime][sizeVect - 2];
////        }
////        prevTime = (prevTime + 1) % 2;
////        currTime = (currTime + 1) % 2;
////    }
//
//
//
////        step = (fabs(xStart) + fabs(xEnd)) / nX;      // calculate step
////    // Calculating
////    time_S = omp_get_wtime();
////
////    for (double j = 0; j < tFinal; j += dt) {
////        multiplicateVector(matrix, vect[prevTime], vect[currTime], nX+2);
////        prevTime = (prevTime + 1) % 2;
////        currTime = (currTime + 1) % 2;
////    }
////    time_E = omp_get_wtime();
////    printf("Run time %.15lf\n", time_E-time_S);
//
////
////    for (double j = 0; j < timesize; j += 1) {
////
////        omp_set_num_threads(4);
////        {
////            #pragma omp parallel for if (ENABLE_PARALLEL)
////            for (int i = 1; i <= nX; i++) {
////                vect[currTime][i] = expr * (vect[prevTime][i + 1] - 2 * vect[prevTime][i] + vect[prevTime][i - 1])
////                                    + vect[prevTime][i];
////            }
////        }
////        // boundary conditions
////        vect[currTime][0] = vect[currTime][1];
////        vect[currTime][nX+1] = vect[currTime][nX];
////
////        prevTime = (prevTime + 1) % 2;
////        currTime = (currTime + 1) % 2;
////
////    }
////
////    time_E = omp_get_wtime();
////    printf("Run time:\t %.15lf\n", time_E-time_S);
//
//    if (rankP == ROOT) {
////        FILE *outfile;
////
////        string outfilename = "../result/Sergey/Sergey_MPI_SparseEuler1D.txt";
////        outfile = fopen(outfilename.c_str(), "w");
////
////        for (int i = 1; i <= nX; i++) {
////            fprintf(outfile, "%2.15le\n", vect[prevTime][i]);
////        }
//    }
    MPI_Finalize();
}
