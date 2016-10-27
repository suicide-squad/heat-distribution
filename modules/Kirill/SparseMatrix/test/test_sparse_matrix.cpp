#include <gtest/gtest.h>
#include <sparse_matrix.h>

#include "sparse_matrix.h"

TEST(can_init_matrix, SpareMatrix) {
  // Arrange
  TYPE value[] = {1, 2, 3, 4, 8, 5, 7, 1, 6};
  int col[] = {0, 2, 4, 3, 3, 5, 1, 2, 5};
  int rowIndex[] = {0, 2, 4, 4, 6, 6, 9};

  // Act
  SpareMatrix s(value, col, rowIndex, 9, 6);

  // Assert
  //EXPECT_EQ(*value, *s.getValue());
}


TEST(can_mult_matrix_and_vector, SpareMatrix) {
  // Arrange
  const int N = 6;
  TYPE value[] = {1, 2, 3, 4, 8, 5, 7, 1, 6};
  int col[] = {0, 4, 2, 3, 3, 5, 1, 2, 5};
  int rowIndex[] = {0, 2, 4, 4, 6, 6, 9};
  SpareMatrix s(value, col, rowIndex, 9, 6);

  vector v = {1, 2, 4, 1, 0, 3};

  // Act
  vector res = s*v;

  // Assert
  vector expected_v = {1, 16, 0, 23, 0, 36};
  EXPECT_EQ(expected_v, res);
}

TEST(can_mult_matrix_m_n_and_vector, SpareMatrix) {
  // Arrange
  const int m = 5;
  const int n = 4;
  const int nz = 6;
  TYPE value[nz] = {3, 3, 1, 1, 1, 1};
  int col[nz] = {0, 1, 2, 3, 0, 3};
  int rowIndex[m + 1] = {0, 1, 3, 3, 4, 6};
  SpareMatrix s(value, col, rowIndex, nz, m);

  vector v = {1, 2, 0, 1};

  // Act
  vector res = s*v;

  // Assert
  vector expected_v = {3, 6, 0, 1, 2};
  EXPECT_EQ(expected_v, res);
}

TEST(can_add_vectors, vector) {
  // Arrange
  vector v1 = {1, 2, 0, 1};
  vector v2 = {4, 5, 7, 1};

  // Act
  vector res = v1 + v2;

  // Assert
  vector expected_v = {5, 7, 7, 2};
  EXPECT_EQ(expected_v, res);
}


TEST(can_mult_vectors_and_number, vector) {
  // Arrange
  vector v = {1, 2, 0, 1};
  TYPE d = 3;

  // Act
  vector res = v * d;

  // Assert
  vector expected_v = {3, 6, 0, 3};
  EXPECT_EQ(expected_v, res);
}
