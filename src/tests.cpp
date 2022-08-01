#include <gtest/gtest.h>
#include "s21_matrix_oop.h"

TEST(Constructor, DefaultConstructor) {
  S21Matrix firstMatrix;
  EXPECT_EQ(firstMatrix.GetRows(), 1);
  EXPECT_EQ(firstMatrix.GetCols(), 1);
}

TEST(Constructor, ConstructorWithTwoArguments) {
  S21Matrix firstMatrix(1, 1);
  EXPECT_EQ(firstMatrix.GetRows(), 1);
  EXPECT_EQ(firstMatrix.GetCols(), 1);

  S21Matrix secondMatrix(22, 22);
  EXPECT_EQ(secondMatrix.GetRows(), 22);
  EXPECT_EQ(secondMatrix.GetCols(), 22);

  S21Matrix thirdMatrix(47, 41);
  EXPECT_EQ(thirdMatrix.GetRows(), 47);
  EXPECT_EQ(thirdMatrix.GetCols(), 41);
}

TEST(Constructor, ConstructorCopy) {
  S21Matrix firstMatrix(4, 4);
  S21Matrix secondMatrix(firstMatrix);

  EXPECT_EQ(firstMatrix.GetRows(), 4);
  EXPECT_EQ(secondMatrix.GetRows(), 4);

  EXPECT_EQ(firstMatrix.GetCols(), 4);
  EXPECT_EQ(secondMatrix.GetCols(), 4);

  for (int i = 0; i < firstMatrix.GetRows(); i++) {
      for (int j = 0; j < firstMatrix.GetCols(); j++) {
          EXPECT_EQ(firstMatrix(i, j), secondMatrix(i, j));
      }
  }
}

TEST(Constructor, ConstructorMove) {
  S21Matrix firstMatrix(4, 4);
  S21Matrix secondMatrix(firstMatrix);

  EXPECT_EQ(firstMatrix.GetRows(), 4);
  EXPECT_EQ(secondMatrix.GetRows(), 4);

  EXPECT_EQ(firstMatrix.GetCols(), 4);
  EXPECT_EQ(secondMatrix.GetCols(), 4);

  for (int i = 0; i < firstMatrix.GetRows(); i++) {
      for (int j = 0; j < firstMatrix.GetCols(); j++) {
          EXPECT_EQ(firstMatrix(i, j), secondMatrix(i, j));
      }
  }
}

TEST(Accessor, GetRows) {
  S21Matrix firstMatrix(5, 2);
  EXPECT_EQ(firstMatrix.GetRows(), 5);
}

TEST(Accessor, GetColumns) {
  S21Matrix firstMatrix(5, 2);
  EXPECT_EQ(firstMatrix.GetCols(), 2);
}

TEST(Accessor, GetMatrix) {
  S21Matrix firstMatrix(3, 3);
  EXPECT_EQ(firstMatrix(0, 1), 0);
}

TEST(Mutator, SetRows) {
  S21Matrix firstMatrix(3, 3);
  firstMatrix.SetRows(2);
  EXPECT_EQ(firstMatrix.GetRows(), 2);
  EXPECT_EQ(firstMatrix.GetCols(), 3);
}

TEST(Mutator, SetColumns) {
  S21Matrix firstMatrix(3, 3);
  firstMatrix.SetColumns(2);
  EXPECT_EQ(firstMatrix.GetCols(), 2);
  EXPECT_EQ(firstMatrix.GetRows(), 3);
}

TEST(EqMatrix, EqTest1) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(3, 3);
  EXPECT_TRUE(firstMatrix.eq_matrix(secondMatrix));
}

TEST(EqMatrix, EqTest2) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(3, 2);
  EXPECT_FALSE(firstMatrix.eq_matrix(secondMatrix));
}

TEST(SumMatrix, SumTest1) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(3, 3);
  firstMatrix.sum_matrix(secondMatrix);
  EXPECT_EQ(firstMatrix(0, 0), 0.0);
  EXPECT_EQ(firstMatrix(0, 1), 0.0);
}

TEST(SumMatrix, SumTest2) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(3, 2);
  EXPECT_THROW(firstMatrix.sum_matrix(secondMatrix), std::invalid_argument);
}

TEST(SubMatrix, SubTest1) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(3, 3);
  firstMatrix.sum_matrix(secondMatrix);
  EXPECT_EQ(firstMatrix(0, 0), 0.0);
  EXPECT_EQ(firstMatrix(0, 1), 0.0);
}

TEST(SubMatrix, SubTest2) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(3, 2);
  EXPECT_THROW(firstMatrix.sum_matrix(secondMatrix), std::logic_error);
}

TEST(MulMatrix, MulTest1) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(3, 3);
  firstMatrix.mul_matrix(secondMatrix);
  EXPECT_EQ(firstMatrix(0, 1), 0.0);
  EXPECT_EQ(firstMatrix.GetRows(), 3);
  EXPECT_EQ(firstMatrix.GetCols(), 3);
}

TEST(MulMatrix, MulTest2) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(2, 3);
  EXPECT_THROW(firstMatrix.mul_matrix(secondMatrix), std::logic_error);
}

TEST(MulMatrix, MulTest3) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(3, 2);
  firstMatrix.mul_matrix(secondMatrix);
  EXPECT_EQ(firstMatrix(0, 1), 0.0);
  EXPECT_EQ(firstMatrix.GetRows(), 3);
  EXPECT_EQ(firstMatrix.GetCols(), 2);
}

TEST(MulNumber, MulNumTest) {
  S21Matrix firstMatrix(3, 3);
  firstMatrix.mul_number(2);
  EXPECT_EQ(firstMatrix(0, 1), 0.0);
}

TEST(Transpose, TransposeTest) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix(firstMatrix.transpose());
  EXPECT_EQ(secondMatrix(0, 1), 0.0);
}

TEST(Determinant, DeterminantTest) {
  S21Matrix firstMatrix(3, 3);
  double determinant = firstMatrix.determinant();
  ASSERT_EQ(determinant, 0);
}

TEST(CalcComplements, ComplTest1) {
  S21Matrix firstMatrix(3, 3);
  S21Matrix secondMatrix(firstMatrix.calc_complements());
  EXPECT_EQ(secondMatrix(2, 2), 0.0);
}

TEST(CalcComplements, ComplTest2) {
  S21Matrix firstMatrix(3, 2);
  EXPECT_THROW(firstMatrix.calc_complements(), std::logic_error);
}

TEST(InverseMatrix, InverseTest) {
  S21Matrix firstMatrix(3, 3);
  EXPECT_THROW(S21Matrix secondMatrix(firstMatrix.inverse_matrix()), std::logic_error);
}

TEST(Operator, Plus) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix(2, 2);
  S21Matrix thirdMatrix = firstMatrix + secondMatrix;
  EXPECT_EQ(thirdMatrix(0, 0), 0.0);
}

TEST(Operator, Minus) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix(2, 2);
  S21Matrix thirdMatrix = firstMatrix - secondMatrix;
  EXPECT_EQ(thirdMatrix(0, 0), 0.0);
}

TEST(Operator, MulMatrix) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix(2, 2);
  S21Matrix thirdMatrix = firstMatrix * secondMatrix;
  EXPECT_EQ(thirdMatrix(0, 0), 0.0);
  EXPECT_EQ(thirdMatrix.GetRows(), 2);
  EXPECT_EQ(thirdMatrix.GetCols(), 2);
}

TEST(Operator, MulNum) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix = firstMatrix * 2;
  EXPECT_EQ(secondMatrix(0, 0), 0.0);
}

TEST(Operator, Equal) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix(2, 2);
  EXPECT_TRUE(firstMatrix == secondMatrix);
}

TEST(Operator, PlusEqual) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix(2, 2);
  firstMatrix += secondMatrix;
  EXPECT_EQ(firstMatrix(0, 0), 0.0);
}

TEST(Operator, MinusEqual) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix(2, 2);
  firstMatrix -= secondMatrix;
  EXPECT_EQ(firstMatrix(0, 0), 0.0);
}

TEST(Operator, MulNumEqual) {
  S21Matrix firstMatrix(2, 2);
  firstMatrix *= 2;
  EXPECT_EQ(firstMatrix(0, 0), 0.0);
}

TEST(Operator, MulMatEqual) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix(2, 2);
  firstMatrix *= secondMatrix;
  EXPECT_EQ(firstMatrix(0, 0), 0.0);
  EXPECT_EQ(firstMatrix.GetRows(), 2);
  EXPECT_EQ(firstMatrix.GetCols(), 2);
}

TEST(Operator, Ravno) {
  S21Matrix firstMatrix(2, 2);
  S21Matrix secondMatrix = firstMatrix;
  EXPECT_EQ(secondMatrix.GetRows(), 2);
  EXPECT_EQ(secondMatrix.GetCols(), 2);
  EXPECT_EQ(secondMatrix(0, 0), 0.0);
}
