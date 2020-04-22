/*
Author: R Raghu Raman
*/
#ifndef _COO_H
#define _COO_H

#include "matrix.h"
#include <map>
#include <utility>
typedef unsigned long long ull;


/*The comparison functor for sorting the matrix entries - first by row, then by column*/
class LeftBottom
{
public:
  const bool operator()(const std::pair <ull,ull>& a, const std::pair <ull,ull>& b) const
  {
    return ((a.first < b.first) || ( (a.first == b.first) && (a.second < b.second) ) );
  }
};

template <class U>
class COO
{
protected:
  /*Storing the entries of the sparse matrix in a map (Self balancing red black tree)*/
  std::map < std::pair <ull, ull> , U, LeftBottom > Coordlist;

  /*Storing the row and column sizes*/
  ull matrow,matcol;
public:
  /*various ways to initialize a Coordinate Format sparse matrix*/
  COO();
  COO(const ull&);
  COO(const ull& ,const ull&);
  COO(const COO&);
  COO(const matrix<U>&);

  /* The subscript operator to access matrix entries*/
  U& operator ()(const ull&,const ull&);
  const U operator ()(const ull&,const ull&) const;

  /*The basic arithmetic operations on two matrices*/
  COO operator +(const COO& ) const;
  COO& operator +=(const COO&);
  COO operator -(const COO& ) const;
  COO& operator -=(const COO& );
  COO operator *(const COO& ) const;
  COO& operator *=(const COO&);
  COO operator -() const;

  /*Scalar Multiplication of matrix*/
  COO operator *(const U&) const;
  COO& operator *=(const U&);

  /*To create an instance of the identitiy matrix*/
  static COO iden(const ull&);

  /*The transpose of a matrix*/
  COO transpose() const;

  /*Resizing matrices*/
  void resize(const ull& ,const ull&);
  void resize(const ull&);

  /*To get the row size and column sizes of the matrix*/
  ull row_size() const;
  ull col_size() const;
};

#endif
