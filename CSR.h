/*
Author: R Raghu Raman
*/
#ifndef _CSR_H
#define _CSR_H

#include "matrix.h"
#include <vector>
typedef unsigned long long ull;


template <class U>
class CSR
{
protected:
  /*Storing the entries of the sparse matrix in 3 arrays as per the CSR format*/
  std::vector <U> _Entry;
  std::vector <ull> _IEntry,_JEntry;
  /*Storing the row and column sizes*/
  ull matrow,matcol;
public:
  /*various ways to initialize a Compressed sparse row Format sparse matrix*/
  CSR();
  CSR(const ull&);
  CSR(const ull& ,const ull&);
  CSR(const CSR&);
  CSR(const matrix<U>& );

  /* The subscript operator to access matrix entries*/
  const U operator ()(const ull& ,const ull&) const;
  void insert(const ull& ,const ull& ,const U&);

  /*The basic arithmetic operations on two matrices*/
  CSR operator +(const CSR& ) const;
  CSR& operator +=(const CSR&);
  CSR operator -(const CSR& ) const;
  CSR& operator -=(const CSR& );
  CSR operator *(const CSR& ) const;
  CSR& operator *=(const CSR&);
  CSR operator -() const;

  /*Scalar Multiplication of matrix*/
  CSR operator *(const U&) const;
  CSR& operator *=(const U&);

  /*To create an instance of the identitiy matrix*/
  static CSR iden(const ull&);

  /*The transpose of a matrix*/
  CSR transpose() const;

  /*Resizing matrices*/
  void resize(const ull&,const ull&);
  void resize(const ull&);

  /*To get the row size and column sizes of the matrix*/
  ull row_size() const;
  ull col_size() const;
};

#endif
