#include "s21_matrix_oop.h"

/**
 * @brief Construct a new S21Matrix::S21Matrix object
 * 
 */
S21Matrix::S21Matrix() :
    _rows(1),
    _cols(1) {
    init_matrix(_rows, _cols);
}

/**
 * @brief Destroy the S21Matrix::S21Matrix object
 * 
 */
S21Matrix::~S21Matrix() {
    free_matrix();
}

/**
 * @brief Construct a new S21Matrix::S21Matrix object
 * 
 * @param rows Count of rows
 * @param cols Count of columns
 */
S21Matrix::S21Matrix(int rows, int cols) :
    _rows(rows),
    _cols(cols) {
    if (rows < 1 || cols < 1) {
        throw std::logic_error("\nWrong value of rows or columns\n");
    }
    init_matrix(_rows, _cols);
}

/**
 * @brief Construct a new S21Matrix::S21Matrix object
 * 
 * @param other_matrix Other Matrix object for copy
 */
S21Matrix::S21Matrix(const S21Matrix& other_matrix) :
    _rows(other_matrix._rows), 
    _cols(other_matrix._cols) {
    init_matrix(_rows, _cols);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _cols; j++) {
            _matrix[i][j] = other_matrix._matrix[i][j];
        }
    }
}

/**
 * @brief Construct a new S21Matrix::S21Matrix object
 * 
 * @param other_matrix Other Matrix object for move
 */
S21Matrix::S21Matrix(S21Matrix&& other_matrix) {
    _rows = other_matrix._rows;
    _cols = other_matrix._cols;
    _matrix = other_matrix._matrix;
    other_matrix._cols = other_matrix._rows = 0;
    other_matrix._matrix = nullptr;
}

/**
 * @brief Initialize new matrix
 * 
 * @param rows Count of rows
 * @param cols Count of columns
 */
void S21Matrix::init_matrix(int rows, int cols) {
    _matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        _matrix[i] = new double[cols]();
    }
}

/**
 * @brief Clear matrix
 * 
 */
void S21Matrix::free_matrix() {
    if (_matrix) {
        for (int i = 0; i < _rows; i++) {
            delete[] _matrix[i];
        }
        delete[] _matrix;
        _matrix = nullptr;
    }
}

/**
 * @brief Compare two matrices for identity
 * 
 * @param other_matrix  Other matrix for compare
 * @return True if the matrices are identical
 * @return False if the matrices are different
 */
bool S21Matrix::eq_matrix(const S21Matrix& other_matrix) {
    static const double EPS = 0.0000001;
    if (valid_matrix(other_matrix) && valid_matrix(*this)
    && other_matrix._rows == _rows && other_matrix._cols == _cols) {
        for (int i = 0; i < _rows; i++) {
            for (int j = 0; j < _cols; j++) {
                if (fabs(_matrix[i][j] - other_matrix._matrix[i][j]) > EPS) {
                    return false;
                }
            }
        }
    } else {
        return false;
    }
    return true;
}

/**
 * @brief Sum two matrices
 * 
 * @param other_matrix Other matrix for sum
 */
void S21Matrix::sum_matrix(const S21Matrix& other_matrix) {
    if (_rows != other_matrix._rows || _cols != other_matrix._cols) {
        throw std::invalid_argument("\nThe number of rows and columns must match\n");
    } else {
        for (int i = 0; i < _rows; i++) {
            for (int j = 0; j < _cols; j++) {
                _matrix[i][j] += other_matrix._matrix[i][j];
            }
        }
    }
}

/**
 * @brief Sub two matrices
 * 
 * @param other_matrix Other matrix for sub
 */
void S21Matrix::sub_matrix(const S21Matrix& other_matrix) {
    if (valid_matrix(other_matrix) && valid_matrix(*this) && compare_two_matrix(other_matrix)) {
        for (int i = 0; i < _rows; i++) {
            for (int j = 0; j < _cols; j++) {
                _matrix[i][j] -= other_matrix._matrix[i][j];
            }
        }
    }
}

/**
 * @brief Multiplies a matrix by a specific number
 * 
 * @param Num Number value
 */
void S21Matrix::mul_number(const double num) {
    if (valid_matrix(*this)) {
        for (int i = 0; i < _rows; i++) {
            for (int j = 0; j < _cols; j++) {
                _matrix[i][j] *= num;
            }
        }
    }
}

/**
 * @brief Multiplies a matrix by another matrix
 * 
 * @param other_matrix Other matrix for multiply
 */
void S21Matrix::mul_matrix(const S21Matrix& other_matrix) {
    if ((_cols != other_matrix._rows) || !valid_matrix(*this) || !valid_matrix(other_matrix)) {
        throw std::logic_error("\nWrong count of rows or columns\n");
    }
    int row = _rows;
    int cols = other_matrix._cols;
    S21Matrix tmpMatrix(_rows, other_matrix._cols);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < other_matrix._cols; j++) {
            for (int k = 0; k < _cols; k++) {
                tmpMatrix._matrix[i][j] += _matrix[i][k] * other_matrix._matrix[k][j];
            }
        }
    }
    free_matrix();
    _matrix = tmpMatrix._matrix;
    _rows = row;
    _cols = cols;
    tmpMatrix.null_object_field();
}

/**
 * @brief Matrix transpose
 * 
 * @return S21Matrix result matrix
 */
S21Matrix S21Matrix::transpose() {
    S21Matrix resultMatrix(_rows, _cols);
    if (valid_matrix(*this)) {
        for (int i = 0; i < resultMatrix._rows; i++) {
            for (int j = 0; j < resultMatrix._cols; j++) {
                resultMatrix._matrix[i][j] = _matrix[j][i];
            }
        }
    }
    return resultMatrix;
}

/**
 * @brief Creates a matrix of algebraic complements
 * 
 * @return S21Matrix returns the finished matrix
 */
S21Matrix S21Matrix::calc_complements() {
    S21Matrix resultMatrix(_rows, _cols);
    S21Matrix tmpMatrix(_rows - 1, _cols - 1);
    if (valid_matrix(*this) && is_matrix_square(*this)) {
        for (int i = 0; i < _rows; i++) {
            for (int j = 0; j < _cols; j++) {
                tmpMatrix.get_matrix_minor(*this, i, j);
                resultMatrix._matrix[i][j] = tmpMatrix.determinant() * pow(-1, i + j);
            }
        }
    }
    return resultMatrix;
}

/**
 * @brief Finds the determinant
 * 
 * @return Double determinant
 */
double S21Matrix::determinant() {
    double result = 0.0;
    if (valid_matrix(*this) && is_matrix_square(*this)) {
        for (int j = 0; j < _cols; j++) {
            if (_rows == 1) {
                result = _matrix[0][0];
            } else if (_rows == 2) {
                result = _matrix[0][0] * _matrix[1][1] -
                    _matrix[1][0] * _matrix[0][1];
            } else {
                S21Matrix tmpMatrix(_rows - 1, _cols - 1);
                tmpMatrix.get_matrix_minor(*this, 0, j);
                result += _matrix[0][j] * pow(-1, j) * tmpMatrix.determinant();
            }
        }
    }
    return result;
}

/**
 * @brief Creates an inverse matrix
 * 
 * @return S21Matrix Returns the finished matrix
 */
S21Matrix S21Matrix::inverse_matrix() {
    if (_rows != _cols) {
    throw std::logic_error("\nRows and columns must match\n");
    } else if (determinant() == 0.0) {
    throw std::logic_error(
        "\ndeterminant value can't be equal to 0\n");
    }
    const double determinant = 1 / this->determinant();
    S21Matrix tmpMatrix(calc_complements());
    S21Matrix resultMatrix(tmpMatrix.transpose());
    resultMatrix.mul_number(determinant);
    return resultMatrix;
}

/**
 * @brief Gets the truncated matrix
 * 
 * @param other_matrix Other matrix
 * @param rows Count of rows
 * @param cols Count of columns
 * @return S21Matrix 
 */
S21Matrix S21Matrix::get_matrix_minor(const S21Matrix& other_matrix, int rows,
                                          int cols) {
  for (int i = 0; i < other_matrix._rows; i++) {
    for (int j = 0; j < other_matrix._cols; j++) {
      if (i != rows && j != cols) {
        if (j > cols && i > rows) {
          _matrix[i - 1][j - 1] = other_matrix._matrix[i][j];
        } else if (j < cols && i < rows) {
          _matrix[i][j] = other_matrix._matrix[i][j];
        } else if (j > cols && i < rows) {
          _matrix[i][j - 1] = other_matrix._matrix[i][j];
        } else if (j < cols && i > rows) {
          _matrix[i - 1][j] = other_matrix._matrix[i][j];
        }
      }
    }
  }
  return *this;
}

/**
 * @brief Is the matrix square
 * 
 * @param other_matrix Other matrix for check 
 * @return true if matrix is square
 * @return false if matrix is not square
 */
bool S21Matrix::is_matrix_square(const S21Matrix& other_matrix) {
    if (other_matrix._rows != other_matrix._cols) {
        throw std::logic_error("\nMatrix is not square\n");
    }
    return true;
}

/**
 * @brief Check matrices for correct values
 * 
 * @param other_matrix Other matrix for check
 * @return True if matrices are correct
 * @return False if matrices are incorrect
 */
bool S21Matrix::valid_matrix(const S21Matrix& other_matrix) {
    if ((other_matrix._matrix == nullptr) || (other_matrix._rows <= 0) || (other_matrix._cols <= 0)
    || (other_matrix._rows == 1 && other_matrix._cols == 1)) {
        throw std::logic_error("\nWrong value of some class field\n");
    }
    return true;
}

/**
 * @brief Null object
 * 
 * @param other_matrix Object
 */
void S21Matrix::null_object_field() {
    _rows = _cols = 0;
    _matrix = nullptr;
}

/**
 * @brief Compare two matrices for identity
 * 
 * @param other_matrix Other matrix for compare
 * @return True if matrices are identical
 * @return False if matrices are different
 */
bool S21Matrix::compare_two_matrix(const S21Matrix& other_matrix) {
    if (other_matrix._rows != _rows && other_matrix._cols != _cols) {
        throw std::logic_error("\nMatrices are non-identical\n");
    }
    return true;
}

/**
 * @brief Get count of matrix rows
 * 
 * @return Int matrix rows value
 */
int S21Matrix::GetRows() {
    return _rows;
}

/**
 * @brief Set rows value
 * 
 * @param rows Rows value
 */
void S21Matrix::SetRows(int rows) {
    if (rows < 1) {
        throw std::logic_error("\nRows value can't be less than 1\n");
    }
    S21Matrix tmpMatrix(rows, _cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < _cols; j++) {
            if (i >= _rows) {
                tmpMatrix._matrix[i][j] = 0.0;
            } else {
                tmpMatrix._matrix[i][j] = _matrix[i][j];
            }
        }
    }
    free_matrix();
    _matrix = tmpMatrix._matrix;
    _rows = rows;
    tmpMatrix.null_object_field();
}

/**
 * @brief Set columns value
 * 
 * @param cols Columns
 */
void S21Matrix::SetColumns(int cols) {
    if (cols < 1) {
        throw std::logic_error("\nCols value can't be less than 1\n");
    }
    S21Matrix tmpMatrix(_rows, cols);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (j >= _cols) {
                tmpMatrix._matrix[i][j] = 0.0;
            } else {
                tmpMatrix._matrix[i][j] = _matrix[i][j];
            }
        }
    }
    free_matrix();
    _matrix = tmpMatrix._matrix;
    _cols = cols;
    tmpMatrix.null_object_field();
}

/**
 * @brief Get count of matrix columns
 * 
 * @return Int matrix columns value
 */
int S21Matrix::GetCols() {
    return _cols;
}

/**
 * @brief Operator sum equals sign overload
 * 
 * @param other_matrix Matrix object
 */
void S21Matrix::operator+=(const S21Matrix& other_matrix) {
    sum_matrix(other_matrix);
}

/**
 * @brief Operator sub equals sign overload
 * 
 * @param other_matrix Matrix object
 */
void S21Matrix::operator-=(const S21Matrix& other_matrix) {
    sub_matrix(other_matrix);
}

/**
 * @brief Operator mul equals sign overload
 * 
 * @param other_matrix Matrix object
 */
void S21Matrix::operator*=(const S21Matrix& other_matrix) {
    mul_matrix(other_matrix);
}

/**
 * @brief Operator mul equals sign overload
 * 
 * @param num Number value
 */
void S21Matrix::operator*=(double num) {
    mul_number(num);
}

/**
 * @brief Operator equals sign oerload
 * 
 * @param other_matrix Matrix object
 */
void S21Matrix::operator=(S21Matrix&& other_matrix) {
    std::swap(*this, other_matrix);
}

/**
 * @brief Operator parentheses overload
 * 
 * @param rows Rows value
 * @param cols Columns value
 * @return Double& matrtx value
 */
double& S21Matrix::operator()(int rows, int cols) {
    if (rows >= _rows || cols >= _cols) {
        throw std::logic_error("\nIndex out of range\n");
    }
    return _matrix[rows][cols];
}

/**
 * @brief Operator double equals sign overload
 * 
 * @param other_matrix Matrix object
 * @return True if matrices are idntity
 * @return False if matrices are different
 */
bool S21Matrix::operator==(const S21Matrix& other_matrix) {
    return eq_matrix(other_matrix);
}

/**
 * @brief Operator sum overload
 * 
 * @param other_matrix Matrix object
 * @return S21Matrix result of matrix
 */
S21Matrix S21Matrix::operator+(const S21Matrix& other_matrix) {
    S21Matrix resultMatrix(*this);
    resultMatrix.sum_matrix(other_matrix);
    return resultMatrix;
}

/**
 * @brief Operator sub overload
 * 
 * @param other_matrix Matrix object
 * @return S21Matrix result of matrix
 */
S21Matrix S21Matrix::operator-(const S21Matrix& other_matrix) {
    S21Matrix resultMatrix(*this);
    resultMatrix.sub_matrix(other_matrix);
    return resultMatrix;
}

/**
 * @brief Operator mul overload
 * 
 * @param other_matrix Matrix object
 * @return S21Matrix result of matrix
 */
S21Matrix S21Matrix::operator*(const S21Matrix& other_matrix) {
    S21Matrix resultMatrix(*this);
    resultMatrix.mul_matrix(other_matrix);
    return resultMatrix;
}

/**
 * @brief Operator mul overload
 * 
 * @param num Number value
 * @return S21Matrix result of Matrix
 */
S21Matrix S21Matrix::operator*(double num) {
    S21Matrix resultMatrix(*this);
    resultMatrix.mul_number(num);
    return resultMatrix;
}
