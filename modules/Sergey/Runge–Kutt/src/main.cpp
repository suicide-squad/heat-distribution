#include <iostream>
#include <cmath>

using std::string;

int main(int argc, char** argv) {
    string filename = "INPUT.txt";
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

    printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinal %lf; dt %lf;\n",
           xStart, xEnd, sigma, nX, tStart, tFinal, dt);

    double** vect = new double*[2];
    vect[0] = new double[nX + 2];
    vect[1] = new double[nX + 2];

    // Read file
    for (int i = 1; i <= nX; i++) {
        fscanf(infile, "%lf\n", &vect[0][i]);
    }
    fclose(infile);

    step = (fabs(xStart) + fabs(xEnd)) / nX;      // calculate step

    prevTime = 0;
    currTime = 1;

    vect[0][0] = vect[0][1];
    vect[0][nX+1] = vect[0][nX];

    double* k1Vect = new double[nX + 2];
    double* k2Vect = new double[nX + 2];
    double* k3Vect = new double[nX + 2];
    double* k4Vect = new double[nX + 2];

    double expression = sigma / step * step;
    double expression2 = sigma / (step + step / 2) * (step + step / 2); // h + h/2
    double expression3 = sigma / (step + step) * (step + step);

    for (double j = 0; j < tFinal; j += dt) {
        // Fill k1 vect
        for (int k = 1; k <= nX; ++k) {
            k1Vect[i] = expression * (vect[prevTime][i + 1] - 2 * vect[prevTime][i] + vect[prevTime][i - 1]);
        }
        //bounders k1
        k1Vect[0] = k1Vect[1];
        k1Vect[nX + 1] = k1Vect[nX];

        // Fill k2 vect
        for (int k = 1; k <= nX; ++k) {
            k2Vect[i] = expression2 * (vect[prevTime][i + 1]) + (k1Vect[i + 1] * step / 2)
                        - 2 * vect[prevTime][i] + k1Vect[i] * step / 2
                        + vect[prevTime][i - 1] + k1Vect[i - 1] * step / 2;
        }
        // bounders
        k2Vect[0] = k2Vect[1];
        k2Vect[nX + 1] = k2Vect[nX];

        // Fill k3 vect
        for (int k = 1; k <= nX; ++k) {
            k3Vect[i] = expression2 * (vect[prevTime][i + 1]) + (k2Vect[i + 1] * step / 2)
                        - 2 * vect[prevTime][i] + k2Vect[i] * step / 2
                        + vect[prevTime][i - 1] + k2Vect[i - 1] * step / 2;
        }
        // bounders
        k3Vect[0] = k3Vect[1];
        k3Vect[nX + 1] = k3Vect[nX];

        // Fill k4 vect
        for (int k = 1; k <= nX; ++k) {
            k4Vect[i] = expression3 * (vect[prevTime][i + 1]) + (k3Vect[i + 1] * step)
                        - 2 * vect[prevTime][i] + k3Vect[i] * step
                        + vect[prevTime][i - 1] + k3Vect[i - 1] * step;
        }
        // bounders
        k3Vect[0] = k3Vect[1];
        k3Vect[nX + 1] = k3Vect[nX];

        // FIll result vector
        for (int k = 1; k <= nX; ++k) {
            vect[currTime][i] = vect[prevTime][i] + step / 6 * (k1Vect[i] + 2*k2Vect[i] + 2*k3Vect[i] + k4Vect[i]);
        }

        // boundary conditions
        vect[currTime][0] = vect[currTime][1];
        vect[currTime][nX+1] = vect[currTime][nX];

        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }


    /*for (double j = 0; j < tFinal; j += dt) {
        for (int i = 1; i <= nX; i++) {
            vect[currTime][i] = expression * (vect[prevTime][i + 1] - 2 * vect[prevTime][i] + vect[prevTime][i - 1])
                                + vect[prevTime][i];
        }

        // boundary conditions
        vect[currTime][0] = vect[currTime][1];
        vect[currTime][nX+1] = vect[currTime][nX];

        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;

    }*/

    FILE *outfile = fopen("OUTPUT_Runge.txt", "w");

    for (int i = 1; i <= nX; i++) {
        fprintf(outfile, "%2.15le\n", vect[prevTime][i]);
    }
}