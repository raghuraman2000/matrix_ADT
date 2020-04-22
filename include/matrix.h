/*
Author: R Raghu Raman
*/
#ifndef _MATRIX_H
#define _MATRIX_H

#include <vector>
typedef unsigned long long ull;

#define TOLERANCE 10e-15 /* Tolerance for numeric data types to be considered 0*/

/* All the error handling classes are defined below*/
class INVALID_MATRIX_SUBSCRIPT{};
class INVALID_MATRIX_ADDITION{};
class INVALID_MATRIX_SUBTRACTION{};
class INVALID_MATRIX_MULTIPLICATION{};
class INVALID_MATRIX_DIMENSION{};
class INVALID_MATRIX_BASE_ASSIGNMENT{};
class SINGULAR_MATRIX_ERROR
{
public:
  ull rank;
  SINGULAR_MATRIX_ERROR(ull n = 0):rank(n){};
};

template <class U>
class matrix
{
protected:
  /* A 2D array to store the matrix entries*/
  std::vector <std::vector <U> > matrixgrid;

  /* Two unsigned variables to store the matrix row and column size*/
  ull matrow,matcol;

public:
  /*Various method to create a matrix object*/
  matrix();
  matrix(const ull&);
  matrix(const ull& ,const ull&);
  matrix(const std::vector <std::vector <U> >& );
  matrix(const matrix& );

  /* The subscript operator to access matrix entries*/
  U& operator ()(const ull&,const ull&);
  const U& operator ()(const ull&,const ull&) const;

  /*The basic arithmetic operations on two matrices*/
  matrix operator +(const matrix& ) const;
  matrix& operator +=(const matrix&);
  matrix operator -(const matrix& ) const;
  matrix& operator -=(const matrix& );
  matrix operator *(const matrix& ) const;
  matrix& operator *=(const matrix& );
  matrix operator -() const;

  /*Scalar Multiplication of matrix*/
  matrix operator *(const U&) const;
  matrix& operator *=(const U&);

  /*Pseudo inverse of a matrix with full column rank*/
  matrix pinv() const;

  /*Least squares polynomial fit of a given degree*/
  static std::vector <U> Poly_fit(std::vector < std::pair<U,U> >&, const ull&);

  /*Decomposition of a matrix into an orthogonal matrix and an upper triangular matrix*/
  void QR(matrix&, matrix&,const U tol = TOLERANCE) const;

  /*To calculate the column rank of a matrix*/
  ull col_rank() const;

  /*The transpose of a matrix*/
  matrix transpose() const;

  /*Resizing matrices*/
  void resize(const ull& ,const ull&);
  void resize(const ull&);

  /*To get the row size and column sizes of the matrix*/
  ull row_size() const;
  ull col_size() const;
};

#endif
