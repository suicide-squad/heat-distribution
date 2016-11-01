//
// Created by lenferd on 27.10.16.
//

#include <iostream>
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
    double** original_matrix = nullptr;
    int matrixSize = 0;

    string filename = "InputMatrix.txt";
    fillMatrix(original_matrix, filename, matrixSize);

    SparseMatrix matrix;
    matrix.fillMatrix(original_matrix, matrixSize, matrixSize);
    matrix.printVectors();
/*    // Vector 2
    vector<double> vect = fillVect();

    matrix.multiplicateVector(vect);
    matrix.printVectors();*/



}

vector<double> fillVect() {
    vector<double> vect;
    vect.push_back(2);
    vect.push_back(3);
}