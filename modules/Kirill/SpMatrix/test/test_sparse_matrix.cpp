#include <gtest/gtest.h>

#include "sp_mat.h"

TEST(can_mult_matrix_and_vector, spMatrix) {
  // Arrange
  const int N = 6;
  const int NZ = 9;
  spMatrix A;
  initSpMat(&A, NZ, N);
  TYPE value[] = {1, 2, 3, 4, 8, 5, 7, 1, 6};
  int col[] = {0, 4, 2, 3, 3, 5, 1, 2, 5};
  int rowIndex[] = {0, 2, 4, 4, 6, 6, 9};

  for (int i = 0; i < N + 1; i++)
    A.rowIndex[i] = rowIndex[i];

  for (int i = 0; i < NZ; i++) {
    A.col[i] = col[i];
    A.value[i] = value[i];
  }

  TYPE v[] = {1, 2, 4, 1, 0, 3};

  // Act
  TYPE* result = (TYPE *)malloc(sizeof(TYPE) * N);
  multMV(&result, A, v);

  // Assert
  TYPE expected_v[] = {1, 16, 0, 23, 0, 36};
  for (int i = 0; i < N; i++)
    EXPECT_EQ(expected_v[i], result[i]);
  freeSpMat(&A);
  free(result);
}

TEST(can_sum, vectors) {
  // Arrange
  const int N = 4;
  const double h = 0.5;

  TYPE U[] = {1, 2, 4, 1};
  TYPE k1[] = {1, 1, 5, 0};
  TYPE k2[] = {2, 1, 3, 1};
  TYPE k3[] = {0, 1, 1, 1};
  TYPE k4[] = {3, 3, 1, 0};


  // Act
  TYPE* result = (TYPE *)malloc(sizeof(TYPE) * N);
  sumV(N, h, &result, U, k1, k2, k3, k4);

  // Assert
  TYPE expected_v[] = {5, 6, 11, 3};
  for (int i = 0; i < N; i++)
    EXPECT_EQ(expected_v[i], result[i]);
  free(result);
}