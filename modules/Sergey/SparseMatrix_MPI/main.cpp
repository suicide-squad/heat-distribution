//
// Created by lenferd on 04.04.17.
//

#include <mpi.h>
#include <iostream>
#include <cmath>
#include "SparseMatrix.h"

using std::string;

int printlf(double value){
    printf("%lf\n", value);
}

int printint(int value){
    printf("%d\n", value);
}

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
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    const int LAST_PROC = sizeP - 1;

    printf("Rank %d Size %d\n", rankP, sizeP);


    double xStart = 0.0, xEnd = 0.0;
    double sigma = 0.0;

    int bc = 0; // Not use

    int nX;
    double tStart = 0.0, tFinal = 0.0;
    double dt = 0.0;

    int i, j;
    int sizeTime = 0;

//    double timeStep = 0.0;
    double step = 0.0;

    double *vect;

    if (rankP == ROOT) {
        string filename = "../../initial/INPUT.txt";
//        string filename = "../../initial_test/INPUT.txt";
//        string filename = "../../initial_test/INPUT2.txt";
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

        // if we use add cell, we need add more that for our boundaries conditionals
        int new_size = nX + ADD_CELL * 2;
        vect = new double[new_size];
        // Fill zero
        for (int k = 0; k < new_size; ++k) {
            vect[k] = 0;
        }

        // ..[....]..
        // if we add some additional cell, we should change start reading position
        for (int i = ADD_CELL; i <= nX + ADD_CELL; i++) {
            fscanf(infile, "%lf\n", &vect[i]);
        }
        fclose(infile);

        step = (fabs(xStart) + fabs(xEnd)) / nX;      // calculate step
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // now we should scatter data for all proc
//    MPI_Bcast(&sizeTime, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&nX, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&step, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&tFinal, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&tStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&sigma, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

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

//    if (rankP == LAST_PROC) {
//        for (int i = 0; i < block_size_add; ++i) {
//            printf("%d %lf\n", rankP, procVect[0][i]);
//        }
//    }

    /***
     * Calculating
     */

    int prevTime = 0;
    int currTime = 1;


    double expression = (sigma * dt) / (step * step);
//    double expression = 0;
    double expression2 = 1.0 - 2.0 * expression;
//    double expression2 = 1;

    SparseMatrix sparseMatrix;
    spMatrixInit(sparseMatrix, block_size_add * 3, block_size_add);
    fillMatrix2ExprWithoutBoundaries(sparseMatrix, block_size_add, expression, expression2);


    int fuck_count = 0;

    // Prepare vect for sharing
    double *vectRightSend = new double[ADD_CELL];
//    double *vectRight = new double[ADD_CELL];
    double *vectRightRecv = new double[ADD_CELL];
//    double *vectLeft = new double[ADD_CELL];
    double *vectLeftSend = new double[ADD_CELL];
    double *vectLeftRecv = new double[ADD_CELL];
//
//

    if (rankP == ROOT) {
        startTime = MPI_Wtime();
    }
    for (double t = tStart; t < tFinal; t += dt) {
//         if info on the add cell get old
        if (fuck_count == ADD_CELL) {
            if (sizeP == 1) {
                ;
            } else if (rankP == ROOT) {
                for (int k = 0; k < ADD_CELL; ++k) {
                    vectRightSend[k] = procVect[currTime][block_size + k];
                }
//                MPI_Send(vectRightSend, ADD_CELL, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD);
                MPI_Sendrecv(vectRightSend, ADD_CELL, MPI_DOUBLE, rankP + 1, 0,
                             vectRightRecv, ADD_CELL, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, &status);
                for (int k = 0; k < ADD_CELL; ++k) {
                    procVect[currTime][block_size+ADD_CELL+k] = vectRightRecv[k];
                }
            } else if (rankP == LAST_PROC) {
                for (int k = 0; k < ADD_CELL; ++k) {
                    vectLeftSend[k] = procVect[currTime][ADD_CELL + k];
                }
//                printint(rankP);
//                MPI_Recv(vectRight, ADD_CELL, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, status);
                MPI_Sendrecv(vectLeftSend, ADD_CELL, MPI_DOUBLE, rankP - 1, 0,
                             vectLeftRecv, ADD_CELL, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, &status);
                for (int k = 0; k < ADD_CELL; ++k) {
                    procVect[currTime][k] = vectLeftRecv[k];
                }
            } else {
            for (int k = 0; k < ADD_CELL; ++k) {
                vectRightSend[k] = procVect[currTime][block_size + k];
                vectLeftSend[k] = procVect[currTime][ADD_CELL + k];
            }
            MPI_Sendrecv(vectRightSend, ADD_CELL, MPI_DOUBLE, rankP + 1, 0,
                         vectRightRecv, ADD_CELL, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, &status);
            MPI_Sendrecv(vectLeftSend, ADD_CELL, MPI_DOUBLE, rankP - 1, 0,
                         vectLeftRecv, ADD_CELL, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, &status);
            for (int k = 0; k < ADD_CELL; ++k) {
                procVect[currTime][block_size+ADD_CELL+k] = vectRightRecv[k];
                procVect[currTime][k] = vectLeftRecv[k];
            }
        }
        fuck_count = 0;
        }

        fuck_count++;
        multiplicateVector(sparseMatrix, procVect[prevTime], procVect[currTime], block_size_add);

//        if (rankP == LAST_PROC) {
//            for (int k = 0; k < block_size_add; ++k) {
//                printf("%2.15le\n", procVect[prevTime][k]);
//            }
//        }
//        printf("===========");

        if (rankP == ROOT) {
            procVect[currTime][ADD_CELL-1] = procVect[currTime][ADD_CELL];
        } else if (rankP == LAST_PROC) {
            procVect[currTime][block_size + ADD_CELL] =  procVect[currTime][block_size + ADD_CELL - 1];
        }

        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;

    }

    if (rankP == LAST_PROC) {
        for (int k = ADD_CELL; k < block_size_add - ADD_CELL; ++k) {
//            printf("%2.15le\n", procVect[prevTime][k]);
        }
    }


    int result_vector_size = block_size * sizeP;
    double* result_vector;

    if (rankP == ROOT){
        result_vector = new double[result_vector_size];
    }

//    std::cout << result_vector_size << std::endl;
    MPI_Gather(procVect[0] + ADD_CELL, block_size, MPI_DOUBLE, result_vector, block_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
//
//    if (rankP == ROOT) {
//        for (int m = 0; m < result_vector_size; ++m) {
//            printf("%d\n", m);
//            printf("%lf\n", result_vector[m]);
//        }
//    }
    if (rankP == ROOT & sizeP != 1){
//        printf("%d", result_vector_size);

        endTime = MPI_Wtime();
        printf("Run time:\t %.15lf\n", endTime-startTime);

        FILE *outfile;

        string outfilename = "../../result/Sergey_MPI_SparseEuler1D.txt";
        outfile = fopen(outfilename.c_str(), "w");

        for (int i = 0; i < result_vector_size; i++) {
            fprintf(outfile, "%2.15le\n", result_vector[i]);
        }
    }
    MPI_Finalize();
}
