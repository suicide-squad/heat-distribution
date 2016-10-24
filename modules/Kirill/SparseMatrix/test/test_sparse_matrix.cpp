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

  TYPE vect[N] = {1, 2, 4, 1, 0, 3};

  // Act
  TYPE* res = s*vect;

  // Assert
  TYPE expected_vect[N] = {1, 16, 0, 23, 0, 36};
  EXPECT_TRUE(std::equal(res, res + N, expected_vect, expected_vect + N));
}