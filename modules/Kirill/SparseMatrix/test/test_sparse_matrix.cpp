#include <gtest/gtest.h>
#include <sparse_matrix.h>

#include "sparse_matrix.h"

TEST(first_test, SpareMatrix) {
  SpareMatrix* s = nullptr;
  s = new SpareMatrix();
  EXPECT_NE(s, nullptr);
}