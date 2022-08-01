#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <iostream>
#include <cmath>

class S21Matrix {
 private:
    int _rows, _cols;
    double** _matrix;

    S21Matrix get_matrix_minor(const S21Matrix& OtherMatrix, int rows, int cols);
    void init_matrix(int rows, int cols);
    void free_matrix();
    bool valid_matrix(const S21Matrix& other_matrix);
    bool compare_two_matrix(const S21Matrix& other_matrix);
    bool is_matrix_square(const S21Matrix& other_matrix);
    void null_object_field();

 public:
    S21Matrix();
    S21Matrix(int rows, int cols);
    S21Matrix(const S21Matrix& other_matrix);
    S21Matrix(S21Matrix&& other_matrix);
    ~S21Matrix();

    bool eq_matrix(const S21Matrix& other_matrix);
    void sum_matrix(const S21Matrix& other_matrix);
    void sub_matrix(const S21Matrix& other_matrix);
    void mul_number(const double num);
    void mul_matrix(const S21Matrix& other_matrix);

    double determinant();
    S21Matrix calc_complements();
    S21Matrix inverse_matrix();
    S21Matrix transpose();

    int GetRows();
    int GetCols();
    void SetRows(int rows);
    void SetColumns(int cols);

    void operator+=(const S21Matrix& other_matrix);
    void operator-=(const S21Matrix& other_matrix);
    void operator*=(const S21Matrix& other_matrix);
    void operator*=(double num);
    void operator=(S21Matrix&& other_matrix);
    bool operator==(const S21Matrix& other_matrix);
    double& operator()(int rows, int cols);

    S21Matrix operator+(const S21Matrix& other_matrix);
    S21Matrix operator-(const S21Matrix& other_matrix);
    S21Matrix operator*(const S21Matrix& other_matrix);
    S21Matrix operator*(double num);
};
#endif  // SRC_S21_MATRIX_OOP_H_
