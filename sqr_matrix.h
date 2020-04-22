/*
Author: R Raghu Raman
*/
#ifndef _SQR_MATRIX_H
#define _SQR_MATRIX_H

#include "matrix.h"
#include <vector>
typedef unsigned long long ull;

template <class U>
class sqr_matrix : public matrix<U>
{
public:
  /*Various ways to initialize a square matrix*/
  sqr_matrix();
  sqr_matrix(const ull&);
  sqr_matrix(const std::vector <std::vector <U> >& );
  sqr_matrix(const sqr_matrix& );

  /*Converting a base class matrix to a square matrix if compatible*/
  sqr_matrix& operator =(const matrix<U>&);

  /*To create an instance of the identitiy matrix*/
  static sqr_matrix iden(const ull&);

  /*To get a non-negative integer exponent of a square matrix*/
  sqr_matrix exp(const ull&) const;

  /*Decomposition of a square matrix into a lower triangular matrix,
   an upper triangular matrix and a permutation matrix*/
  void LUP(sqr_matrix&, sqr_matrix&, std::vector <ull>&, const U tol = TOLERANCE) const;

  /*Solves a system of equations given the LUP Decomposition of the coefficient matrix*/
  std::vector<U> LUP_solve(sqr_matrix&, sqr_matrix&, std::vector <ull>&, std::vector<U>&) const;

  /*Solves a system of equations given the coefficient matrix*/
  std::vector<U> eq_solve(std::vector<U>&, const U tol = TOLERANCE) const;

  /*Gets the inverse of a non singular matrix*/
  sqr_matrix inv() const;

  /*Gets the determinant of the matrix*/
  U det() const;

  /*Resizing the square matrix*/
  void resize(const ull&);

  /*To get the size of the square matrix*/
  ull size() const;
};

#endif
