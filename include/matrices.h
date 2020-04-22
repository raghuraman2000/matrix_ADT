#ifndef _MATRICES_H
#define _MATRICES_H

#include "matrix.h"
#include "sqr_matrix.h"
#include "COO.h"
#include "CSR.h"
#include <cmath>

typedef unsigned long long ull;

/*
*The Base Class Matrix is defined below
*/

template <class U>
matrix<U>::matrix():matrow(0),matcol(0)
{

}

template <class U>
matrix<U>::matrix(const ull& n):matrixgrid(n),matrow(n),matcol(n)
{
  for(ull i = 0;i < n;++i)
    matrixgrid[i].resize(n);
}

template <class U>
matrix<U>::matrix(const ull& m ,const ull& n):matrixgrid(m),matrow(m),matcol(n)
{
  for(ull i = 0;i < m;++i)
    matrixgrid[i].resize(n);
}

template <class U>
matrix<U>::matrix(const std::vector <std::vector <U> >& arg1)
{
  for(ull i = 0; i < arg1.size(); ++i)
    if(arg1[i].size() != arg1[0].size())
      throw INVALID_MATRIX_DIMENSION();
  matrow = arg1.size();
  matcol = 0;
  if( matrow > 0)
    matcol = matrixgrid[0].size();
  matrixgrid = arg1;
}

template <class U>
matrix<U>::matrix(const matrix& arg1)
{
  matrixgrid = arg1.matrixgrid;
  matrow = arg1.matrow;
  matcol = arg1.matcol;
}

template <class U>
inline
U& matrix<U>::operator () (const ull& r,const ull& c)
{
  if(r >= matrow || c >= matcol)
    throw INVALID_MATRIX_SUBSCRIPT();
  return matrixgrid[r][c];
}

template <class U>
inline
const U& matrix<U>::operator () (const ull& r,const ull& c) const
{
  if(r >= matrow || c >= matcol)
    throw INVALID_MATRIX_SUBSCRIPT();
  return matrixgrid[r][c];
}

template <class U>
matrix<U> matrix<U>::operator +(const matrix& arg1) const
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    matrix tempmat(matrow,matcol);
    for(ull i = 0; i < matrow; ++i)
      for(ull j = 0; j < matcol; ++j)
        tempmat(i,j) = matrixgrid[i][j] + arg1.matrixgrid[i][j];
    return tempmat;
  }
  else
    throw INVALID_MATRIX_ADDITION();
}

template <class U>
matrix<U>& matrix<U>::operator +=(const matrix& arg1)
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    for(ull i = 0; i < matrow; ++i)
      for(ull j = 0; j < matcol; ++j)
        matrixgrid[i][j] += arg1.matrixgrid[i][j];
  }
  else
    throw INVALID_MATRIX_ADDITION();
  return *this;
}

template <class U>
matrix<U> matrix<U>::operator -(const matrix& arg1) const
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    matrix tempmat(matrow,matcol);
    for(ull i = 0; i < matrow; ++i)
      for(ull j = 0; j < matcol; ++j)
        tempmat(i,j) = matrixgrid[i][j] - arg1.matrixgrid[i][j];
    return tempmat;
  }
  else
    throw INVALID_MATRIX_SUBTRACTION();
}

template <class U>
matrix<U>& matrix<U>::operator -=(const matrix& arg1)
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    for(ull i = 0; i < matrow; ++i)
      for(ull j = 0; j < matcol; ++j)
        matrixgrid[i][j] -= arg1.matrixgrid[i][j];
  }
  else
    throw INVALID_MATRIX_ADDITION();
  return *this;
}

template <class U>
matrix<U> matrix<U>::operator *(const matrix& arg1) const
{
  if(matcol == arg1.matrow)
  {
    matrix tempmat(matrow,arg1.matcol);
    for(ull i = 0; i < matrow; ++i)
    {
      for(ull j = 0; j < arg1.matcol; ++j)
      {
          tempmat(i,j) = matrixgrid[i][0]*(arg1.matrixgrid[0][j]);
          for(ull k = 1; k < matcol; ++k)
            tempmat(i,j) += matrixgrid[i][k]*(arg1.matrixgrid[k][j]);
      }
    }
    return tempmat;
  }
  else
    throw INVALID_MATRIX_MULTIPLICATION();
}

template <class U>
matrix<U>& matrix<U>::operator *=(const matrix& arg1)
{
  return (*this = (*this) * arg1);
}

template <class U>
matrix<U> matrix<U>::operator -() const
{
  matrix tempmat(matrow,matcol);
  for(ull i = 0; i < matrow; ++i)
    for(ull j = 0; j < matcol; ++j)
      tempmat.matrixgrid[i][j] = -matrixgrid[i][j];
  return tempmat;
}

template <class U>
matrix<U> matrix<U>::pinv() const
{
  matrix <U> tempmat = transpose();
  sqr_matrix <U> temp;
  temp = (tempmat)*(*this);
  return temp.inv()*(tempmat);
}

template <class U>
std::vector<U> matrix<U>::Poly_fit(std::vector < std::pair <U,U> >& points, const ull& deg)
{
  ull m = points.size();
  matrix<U> tempmat(m,deg+1);
  for(ull i = 0; i < m; ++i)
  {
    tempmat(i,0) = (U)1;
    for(ull j = 1; j < deg + 1; ++j)
      tempmat(i,j) = tempmat(i,j-1)*(points[i].first);
  }
  tempmat = tempmat.pinv();
  std::vector <U> temp(deg + 1);
  for(ull i = 0; i < deg + 1; ++i)
  {
    temp[i] = (U)0;
    for(ull j = 0; j < m; ++j)
      temp[i] += tempmat(i,j)*(points[j].second);
  }
  return temp;
}

template <class U>
void matrix<U>::QR(matrix& Orth, matrix& Upper,const U tol) const
{
  std::vector <ull> pos;
  Orth.resize(matrow,matcol);
  ull n = matrow, m = matcol;
  for(ull i = 0; i < m; ++i)
  {
    for(ull j = 0; j < n; ++j)
      Orth(j,i) = matrixgrid[j][i];
    for(const ull& x: pos)
    {
      U num = (U)0, denom = (U)0;
      for(ull j = 0; j < n; ++j)
      {
        num += (Orth(j,x)*Orth(j,i));
        denom += (Orth(j,x)*Orth(j,x));
      }
      for(ull j = 0; j < n; ++j)
        Orth(j,i) -= (Orth(j,x)*num)/denom;
    }
    U p = (U)0;
    for(ull j = 0; j < n; ++j)
      p += (Orth(j,i)*Orth(j,i));
    if(p > tol)
    {
      p = sqrt(p);
      for(ull j = 0; j < n; ++j)
        Orth(j,i) /= p;
      pos.push_back(i);
    }

  }
  if(pos.size() < m)
    throw SINGULAR_MATRIX_ERROR(pos.size());
  Upper.resize(matcol,matcol);
  Upper = (Orth.transpose())*(*this);
}

template <class U>
ull matrix<U>::col_rank() const
{
  matrix Orth, Upper;
  try
  {
    QR(Orth,Upper);
    return matcol;
  }
  catch(SINGULAR_MATRIX_ERROR& excep)
  {
    return excep.rank;
  }
}

template <class U>
matrix<U> matrix<U>::transpose() const
{
  matrix tempmat(matcol,matrow);
  for(ull i = 0; i < matrow; ++i)
    for(ull j = 0; j < matcol; ++j)
        tempmat(j,i) = matrixgrid[i][j];
  return tempmat;
}

template <class U>
matrix<U> matrix<U>::operator *(const U& scalar) const
{
  matrix tempmat(matrow,matcol);
  for(ull i = 0; i < matrow; ++i)
    for(ull j = 0; j < matcol; ++j)
      tempmat(i,j) = matrixgrid[i][j]*scalar;
  return tempmat;
}

template <class U>
matrix<U>& matrix<U>::operator *=(const U& scalar)
{
  for(ull i = 0; i < matrow; ++i)
    for(ull j = 0; j < matcol; ++j)
      matrixgrid[i][j] *= scalar;
  return *this;
}

template <class U>
void matrix<U>::resize(const ull& m, const ull& n)
{
  matrow = m;
  matcol = n;
  matrixgrid.resize(m);
  for(ull i = 0; i < m; ++i)
    matrixgrid[i].resize(n);
}

template <class U>
void matrix<U>::resize(const ull& n)
{
  matrow = n;
  matcol = n;
  matrixgrid.resize(n);
  for(ull i = 0; i < n; ++i)
    matrixgrid[i].resize(n);
}

template <class U>
ull matrix<U>::row_size() const
{
  return matrow;
}

template <class U>
ull matrix<U>::col_size() const
{
  return matcol;
}


/*
The Class Square Matrix is defined below
*/

template <class U>
sqr_matrix<U>::sqr_matrix():matrix<U>()
{

}

template <class U>
sqr_matrix<U>::sqr_matrix(const ull& n):matrix<U>(n,n)
{

}

template <class U>
sqr_matrix<U>::sqr_matrix(const std::vector <std::vector <U> >& arg1)
{
  for(ull i = 0; i < arg1.size(); ++i)
    if(arg1[i].size() != arg1.size())
      throw INVALID_MATRIX_DIMENSION();
  matrix<U>::matrow = arg1.size();
  matrix<U>::matcol = 0;
  if( matrix<U>::matrow > 0)
    matrix<U>::matcol = matrix<U>::matrixgrid[0].size();
  matrix<U>::matrixgrid = arg1;
}

template <class U>
sqr_matrix<U>::sqr_matrix(const sqr_matrix& arg1):matrix<U>(arg1)
{

}

template <class U>
sqr_matrix<U>& sqr_matrix<U>::operator =(const matrix<U>& arg1)
{
  if(arg1.row_size() != arg1.col_size())
    throw INVALID_MATRIX_BASE_ASSIGNMENT();
  matrix<U>::matrow = arg1.row_size();
  matrix<U>::matcol = arg1.col_size();
  matrix<U>::matrixgrid.resize(matrix<U>::matrow);
  for(ull i = 0; i < matrix<U>::matrow; ++i)
  {
    matrix<U>::matrixgrid[i].resize(matrix<U>::matcol);
    for(ull j = 0; j < matrix<U>::matcol; ++j)
      matrix<U>::matrixgrid[i][j] = arg1(i,j);
  }
  return *this;
}

template <class U>
sqr_matrix<U> sqr_matrix<U>::iden(const ull& n)
{
  sqr_matrix tempmat(n);
  for(ull i = 0; i < n; ++i)
    for(ull j = 0; j < n; ++j)
      tempmat(i,j) = ((i == j)? U(1) : U(0));
  return tempmat;
}

template <class U>
sqr_matrix<U> sqr_matrix<U>::exp(const ull& index) const
{
  sqr_matrix answer(iden(matrix<U>::matrow));
  sqr_matrix basemat(*this);
  ull tempindex = index;
  while(tempindex > 0)
  {
    if((tempindex&1) == 1)
      answer = answer*basemat;
    basemat = basemat*basemat;
    tempindex >>= 1;
  }
  return answer;
}

template <class U>
void sqr_matrix<U>::resize(const ull& n)
{
  matrix<U>::resize(n);
}

template <class U>
ull sqr_matrix<U>::size() const
{
  return matrix<U>::matrow;
}

template <class U>
void sqr_matrix<U>::LUP(sqr_matrix& Lower, sqr_matrix& Upper, std::vector <ull>& Perm, const U tol) const
{
  ull n = matrix<U>::matrow;
  Lower.resize(n);
  Upper.resize(n);
  Perm.resize(n);
  sqr_matrix<U> tempmat = *this;

  for(ull i = 0; i < n; ++i)
    Perm[i] = i;

  for(ull k = 0; k < n; ++k)
  {
    U p = 0;
    ull k1 = k;
    for(ull i = k; i < n; ++i)
    {
      if((p*p) < (tempmat(i,k)*tempmat(i,k)))
      {
        p = tempmat(i,k);
        k1 = i;
      }
    }
    if(p*p <= tol)
      throw SINGULAR_MATRIX_ERROR();
    ull temp = Perm[k];
    Perm[k] = Perm[k1];
    Perm[k1] = temp;
    for(ull i = 0; i < n; ++i)
    {
      p = tempmat(k,i);
      tempmat(k,i) = tempmat(k1,i);
      tempmat(k1,i) = p;
    }
    for(ull i = k + 1; i < n; ++i)
    {
      for(ull j = k + 1; j < n; ++j)
        tempmat(i,j) -= (tempmat(i,k)*tempmat(k,j))/tempmat(k,k);
      tempmat(i,k) /= tempmat(k,k);
    }
  }


  for(ull i = 0; i < n; ++i)
  {
    for(ull j = 0; j < n; ++j)
    {
      if(i > j)
      {
        Lower(i,j) = tempmat(i,j);
        Upper(i,j) = 0;
      }
      else if(i == j)
      {
        Lower(i,j) = 1;
        Upper(i,j) = tempmat(i,j);
      }
      else
      {
        Lower(i,j) = 0;
        Upper(i,j) = tempmat(i,j);
      }
    }
  }
}

template <class U>
std::vector <U> sqr_matrix<U>::LUP_solve(sqr_matrix& Lower, sqr_matrix& Upper, std::vector <ull>& Perm, std::vector <U>& b) const
{
  ull n = matrix<U>::matrow;
  std::vector <U> answer(n),temp(n);
  for(ull i = 0; i < n; ++i)
  {
    temp[i] = b[Perm[i]];
    for(ull j = 0; j < i; ++j)
      temp[i] -= Lower(i,j)*temp[j];
  }
  for(ull i = 0; i < n; ++i)
  {
    answer[n - i - 1] = temp[n - i - 1];
    for(ull j = n - i; j < n; ++j)
      answer[n - i - 1] -= Upper(n - i - 1,j)*answer[j];
    answer[n - i - 1] /= Upper(n - i - 1,n - i - 1);
  }
  temp.clear();
  return answer;
}

template <class U>
std::vector <U> sqr_matrix<U>::eq_solve(std::vector <U>& b, const U tol) const
{
  sqr_matrix <U> Lower,Upper;
  std::vector <ull> Perm;
  LUP(Lower,Upper,Perm,tol);
  return LUP_solve(Lower,Upper,Perm,b);
}

template <class U>
U sqr_matrix<U>::det () const
{
  sqr_matrix <U> Lower,Upper;
  std::vector <ull> Perm;
  U answer = (U)1;
  ull n = matrix<U>::matrow;
  try
  {
    LUP(Lower,Upper,Perm);
  for(ull i = 0; i < n; ++i)
    answer *= Upper(i,i);
  }
  catch(SINGULAR_MATRIX_ERROR& excep)
  {
    answer = (U)0;
  }
  return answer;
}

template <class U>
sqr_matrix<U> sqr_matrix<U>::inv() const
{
  sqr_matrix <U> Lower,Upper;
  std::vector <ull> Perm;
  LUP(Lower,Upper,Perm);
  ull n = matrix<U>::matrow;
  sqr_matrix <U> inverse(n);
  std::vector <U> answer(n);
  for(ull i = 0; i < n; ++i)
  {
    for(ull j = 0; j < n; ++j)
      answer[j] = ((i == j)? (U)1: (U)0);
    answer = LUP_solve(Lower,Upper,Perm,answer);
    for(ull j = 0; j < n; ++j)
      inverse(j,i) = answer[j];
  }
  return inverse;
}

/*
* The class COO is defined below
*/

template <class U>
COO<U>::COO():matrow(0),matcol(0)
{

}

template <class U>
COO<U>::COO(const ull& n):matrow(n),matcol(n)
{

}

template <class U>
COO<U>::COO(const ull& m, const ull& n):matrow(m),matcol(n)
{

}

template <class U>
COO<U>::COO(const COO& arg1)
{
  matrow = arg1.matrow;
  matcol = arg1.matcol;
  Coordlist = arg1.Coordlist;
}

template <class U>
COO<U>::COO(const matrix<U>& arg1)
{
  matrow = arg1.row_size();
  matcol = arg1.col_size();
  for(ull i = 0; i < matrow; ++i)
  {
    for(ull j = 0; j < matcol; ++j)
    {
      U temp = arg1(i,j);
      if(temp != (U)0)
        Coordlist[std::make_pair(i,j)] = temp;
    }
  }
}

template <class U>
inline
U& COO<U>::operator ()(const ull& r,const ull& c)
{
  if(r >= matrow || c >= matcol)
    throw INVALID_MATRIX_SUBSCRIPT();
  return Coordlist[std::make_pair(r,c)];
}

template <class U>
inline
const U COO<U>::operator ()(const ull& r, const ull& c) const
{
  if(r >= matrow || c >= matcol)
    throw INVALID_MATRIX_SUBSCRIPT();
  typename std::map < std::pair <ull, ull> , U, LeftBottom >::const_iterator it = Coordlist.find(std::make_pair(r,c));
  const U val = ((it == Coordlist.end()) ? (U)0 : it->second);
  return val;
}


template <class U>
COO<U> COO<U>::operator +(const COO& arg1) const
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    COO tempmat(matrow,matcol);
    for(auto y: arg1.Coordlist)
    {
      U temp = this->operator()((y.first).first,(y.first).second) + y.second;
      if(temp != (U)0)
        tempmat.Coordlist[y.first] = temp;
    }
    return tempmat;
  }
  else
    throw INVALID_MATRIX_ADDITION();
}

template <class U>
COO<U>& COO<U>::operator +=(const COO& arg1)
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    for(auto y: arg1.Coordlist)
      Coordlist[y.first] += y.second;
    return *this;
  }
  else
    throw INVALID_MATRIX_ADDITION();
}

template <class U>
COO<U> COO<U>::operator -(const COO& arg1) const
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    COO tempmat(matrow,matcol);
    for(auto y: arg1.Coordlist)
    {
      U temp = this->operator()((y.first).first,(y.first).second) - y.second;
      if(temp != (U)0)
        tempmat.Coordlist[y.first] = temp;
    }
    return tempmat;
  }
  else
    throw INVALID_MATRIX_SUBTRACTION();
}

template <class U>
COO<U>& COO<U>::operator -=(const COO& arg1)
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    for(auto y: arg1.Coordlist)
      Coordlist[y.first] -= y.second;
    return *this;
  }
  else
    throw INVALID_MATRIX_SUBTRACTION();
}

template <class U>
COO<U> COO<U>::operator *(const COO& arg1) const
{
  if(matcol == arg1.matrow)
  {
    COO tempmat(matrow,arg1.matcol);
    for(auto y: arg1.Coordlist)
    {
      ull rowno = (y.first).first, colno = (y.first).second;
      for(ull i = 0; i < matrow; ++i)
      {
        U temp = this->operator()(i,rowno);
        if(temp != (U)0)
          tempmat.Coordlist[std::make_pair(i,colno)] += temp*(y.second);
      }
    }
    return tempmat;
  }
  else
    throw INVALID_MATRIX_MULTIPLICATION();
}

template <class U>
COO<U>& COO<U>::operator *=(const COO& arg1)
{
  return (*this = (*this) * arg1);
}

template <class U>
COO<U> COO<U>::operator -() const
{
  COO<U> tempmat(matrow,matcol);
  for(auto y: Coordlist)
    tempmat.Coordlist[y.first] = -y.second;
  return tempmat;
}

template <class U>
COO<U> COO<U>::operator *(const U& scalar) const
{
  COO tempmat(matrow,matcol);
  if(scalar != (U)0)
  {
    for(auto it = Coordlist.begin();it != Coordlist.end(); ++it)
      tempmat.Coordlist[it->first] = (it->second)*scalar;
  }
  return tempmat;
}

template <class U>
COO<U>& COO<U>::operator *=(const U& scalar)
{
  if(scalar == (U)0)
    Coordlist.clear();
  else
  {
    for(auto it = Coordlist.begin();it != Coordlist.end(); ++it)
      (it->second) *= scalar;
  }
  return *this;
}


template <class U>
COO<U> COO<U>::iden(const ull& n)
{
  COO<U> tempmat(n,n);
  for(ull i = 0; i < n; ++i)
    tempmat(i,i) = (U)1;
  return tempmat;
}

template <class U>
COO<U> COO<U>::transpose() const
{
  COO tempmat(matcol,matrow);
  for(auto y: Coordlist)
    tempmat.Coordlist[std::make_pair((y.first).second,(y.first).first)] = y.second;
  return tempmat;
}

template <class U>
void COO<U>::resize(const ull& m, const ull& n)
{
  matrow = m;
  matcol = n;
  std::map < std::pair <ull, ull> , U, LeftBottom > temp;
  for(auto y: Coordlist)
  {
    if((y.first).first < matrow && (y.first).second < matcol)
      temp[y.first] = y.second;
  }
  Coordlist = temp;
  temp.clear();
}

template <class U>
void COO<U>::resize(const ull& n)
{
  matrow = n;
  matcol = n;
  std::map < std::pair <ull, ull> , U, LeftBottom > temp;
  for(auto y: Coordlist)
  {
    if((y.first).first < matrow && (y.first).second < matcol)
      temp[y.first] = y.second;
  }
  Coordlist = temp;
  temp.clear();
}

template <class U>
ull COO<U>::row_size() const
{
  return matrow;
}

template <class U>
ull COO<U>::col_size() const
{
  return matcol;
}

/*
* The class CSR is defined below
*/

template <class U>
CSR<U>::CSR():matrow(0),matcol(0)
{

}

template <class U>
CSR<U>::CSR(const ull& n):_IEntry(n+1,0),matrow(n),matcol(n)
{

}

template <class U>
CSR<U>::CSR(const ull& m, const ull& n):_IEntry(m+1,0),matrow(m),matcol(n)
{

}

template <class U>
CSR<U>::CSR(const CSR& arg1)
{
  matrow = arg1.matrow;
  matcol = arg1.matcol;
  _Entry = arg1._Entry;
  _IEntry = arg1._IEntry;
  _JEntry = arg1._JEntry;
}

template <class U>
CSR<U>::CSR(const matrix<U>& arg1)
{
  matrow = arg1.row_size();
  matcol = arg1.col_size();
  _IEntry.resize(matrow + 1);
  _IEntry[0] = 0;
  for(ull i = 0; i < matrow; ++i)
  {
    _IEntry[i + 1] = _IEntry[i];
    for(ull j = 0; j < matcol; ++j)
    {
      U temp = arg1(i,j);
      if(temp != (U)0)
      {
        _Entry.push_back(temp);
        _JEntry.push_back(j);
        _IEntry[i + 1]++;
      }
    }
  }
}

template <class U>
void CSR<U>::insert(const ull& r, const ull& c, const U& val)
{
  if(r >= matrow || c >= matcol)
    throw INVALID_MATRIX_SUBSCRIPT();
  ull left = _IEntry[r], right = _IEntry[r+1], pos = _IEntry[r];
  if(right > left)
  {
    pos = lower_bound(_JEntry.begin() + left, _JEntry.begin() + right, c) - _JEntry.begin();
    if(pos < right && _JEntry[pos] == c)
      _Entry[pos] = val;
  }
  _Entry.emplace(_Entry.begin() + pos, val);
  _JEntry.emplace(_JEntry.begin() + pos, c);
  for(ull i = r + 1; i < matrow + 1; ++i)
    _IEntry[i]++;
}


template <class U>
const U CSR<U>::operator ()(const ull& r, const ull& c) const
{
  if(r >= matrow || c >= matcol)
    throw INVALID_MATRIX_SUBSCRIPT();
  ull left = _IEntry[r], right = _IEntry[r+1];
  if(right > left)
  {
    ull pos = lower_bound(_JEntry.begin() + left, _JEntry.begin() + right, c) - _JEntry.begin();
    if(pos < right && _JEntry[pos] == c)
      return _Entry[pos];
  }
  return (U)0;
}


template <class U>
CSR<U> CSR<U>::operator +(const CSR& arg1) const
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    CSR tempmat(matrow,matcol);
    for(ull i = 0; i < matrow; ++i)
    {
      ull beg = _IEntry[i], end = _IEntry[i + 1], beg1 = arg1._IEntry[i], end1 = arg1._IEntry[i + 1];

      while( beg < end && beg1 < end1)
      {
        U temp(_Entry[beg]);
        ull pos(_JEntry[beg]);
        if(_JEntry[beg] < arg1._JEntry[beg1])
          beg++;
        else if(_JEntry[beg] > arg1._JEntry[beg1])
        {
          temp = arg1._Entry[beg1];
          pos = arg1._JEntry[beg1];
          beg1++;
        }
        else
        {
          temp += arg1._Entry[beg1];
          beg++;
          beg1++;
        }
        if(temp != (U)0)
        {
          tempmat._Entry.push_back(temp);
          tempmat._JEntry.push_back(pos);
          tempmat._IEntry[i + 1]++;
        }
      }

      while(beg < end)
      {
        if(_Entry[beg] != (U)0)
        {
          tempmat._Entry.push_back(_Entry[beg]);
          tempmat._JEntry.push_back(_JEntry[beg]);
          tempmat._IEntry[i + 1]++;
        }
        beg++;
      }

      while(beg1 < end1)
      {
        if(arg1._Entry[beg1] != (U)0)
        {
          tempmat._Entry.push_back(arg1._Entry[beg1]);
          tempmat._JEntry.push_back(arg1._JEntry[beg1]);
          tempmat._IEntry[i + 1]++;
        }
        beg1++;
      }

      tempmat._IEntry[i + 1] += tempmat._IEntry[i];
    }
    return tempmat;
  }
  else
    throw INVALID_MATRIX_ADDITION();
}


template <class U>
CSR<U>& CSR<U>::operator +=(const CSR& arg1)
{
  return (*this = *this + arg1);
}

template <class U>
CSR<U> CSR<U>::operator -(const CSR& arg1) const
{
  if(matrow == arg1.matrow && matrow == arg1.matcol)
  {
    CSR tempmat(matrow,matcol);
    for(ull i = 0; i < matrow; ++i)
    {
      ull beg = _IEntry[i], end = _IEntry[i + 1], beg1 = arg1._IEntry[i], end1 = arg1._IEntry[i + 1];

      while( beg < end && beg1 < end1)
      {
        U temp(_Entry[beg]);
        ull pos(_JEntry[beg]);
        if(_JEntry[beg] < arg1._JEntry[beg1])
          beg++;
        else if(_JEntry[beg] > arg1._JEntry[beg1])
        {
          temp = -arg1._Entry[beg1];
          pos = arg1._JEntry[beg1];
          beg1++;
        }
        else
        {
          temp -= arg1._Entry[beg1];
          beg++;
          beg1++;
        }
        if(temp != (U)0)
        {
          tempmat._Entry.push_back(temp);
          tempmat._JEntry.push_back(pos);
          tempmat._IEntry[i + 1]++;
        }
      }

      while(beg < end)
      {
        if(_Entry[beg] != (U)0)
        {
          tempmat._Entry.push_back(_Entry[beg]);
          tempmat._JEntry.push_back(_JEntry[beg]);
          tempmat._IEntry[i + 1]++;
        }
        beg++;
      }

      while(beg1 < end1)
      {
        if(arg1._Entry[beg1] != (U)0)
        {
          tempmat._Entry.push_back(-arg1._Entry[beg1]);
          tempmat._JEntry.push_back(arg1._JEntry[beg1]);
          tempmat._IEntry[i + 1]++;
        }
        beg1++;
      }

      tempmat._IEntry[i + 1] += tempmat._IEntry[i];
    }
    return tempmat;
  }
  else
    throw INVALID_MATRIX_SUBTRACTION();
}


template <class U>
CSR<U>& CSR<U>::operator -=(const CSR& arg1)
{
  return (*this = *this - arg1);
}

template <class U>
CSR<U> CSR<U>::operator -() const
{
  CSR<U> tempmat(matrow,matcol);
  tempmat._IEntry = _IEntry;
  tempmat._JEntry = _JEntry;
  for(ull i = 0; i < _Entry.size(); ++i)
    tempmat._Entry.push_back(-_Entry[i]);
  return tempmat;
}


template <class U>
CSR<U> CSR<U>::operator *(const CSR& arg1) const
{
  if(matcol == arg1.matrow)
  {
    CSR tempmat(matrow,arg1.matcol);
    CSR arg2 = arg1.transpose();
    for(ull i = 0; i < matrow; ++i)
    {
      tempmat._IEntry[i + 1] = tempmat._IEntry[i];
      for(ull j = 0; j < arg2.matrow; ++j)
      {
        ull ii = _IEntry[i], jj = arg2._IEntry[j];
        U temp = (U)0;
        while(ii < _IEntry[i + 1] && jj < arg2._IEntry[j + 1])
        {
          if(_JEntry[ii] == arg2._JEntry[jj])
          {
            temp += _Entry[ii]*(arg2._Entry[jj]);
            ii++;
            jj++;
          }
          else if(_JEntry[ii] < arg2._JEntry[jj])
            ii++;
          else
            jj++;
        }
        if(temp != (U)0)
        {
          tempmat._Entry.push_back(temp);
          tempmat._JEntry.push_back(j);
          tempmat._IEntry[i + 1]++;
        }
      }
    }
    return tempmat;
  }
  else
    throw INVALID_MATRIX_MULTIPLICATION();
}

template <class U>
CSR<U>& CSR<U>::operator *=(const CSR& arg1)
{
  return (*this = (*this)*arg1);
}

template <class U>
CSR<U> CSR<U>::operator *(const U& scalar) const
{
  CSR tempmat(matrow,matcol);
  if(scalar != (U)0)
  {
    tempmat = *this;
    for(ull i = 0; i < _Entry.size(); ++i)
       tempmat._Entry[i] *= scalar;
  }
  return tempmat;
}

template <class U>
CSR<U>& CSR<U>::operator *=(const U& scalar)
{
  if(scalar == (U)0)
  {
    _Entry.clear();
    _JEntry.clear();
    for(ull i = 0; i < matrow + 1; ++i)
      _IEntry[i] = 0;
  }
  else
  {
    for(ull i = 0; i < _Entry[i].begin(); ++i)
      _Entry[i] *= scalar;
  }
  return *this;
}

template <class U>
CSR<U> CSR<U>::iden(const ull& n)
{
  CSR<U> tempmat(n,n);
  for(ull i = 0; i < n; ++i)
  {
    tempmat._IEntry[i + 1] = i + 1;
    tempmat._JEntry.push_back(i);
    tempmat._Entry.push_back((U)1);
  }
  return tempmat;
}


template <class U>
CSR<U> CSR<U>::transpose() const
{
  CSR tempmat(matcol,matrow);
  tempmat._JEntry.resize(_JEntry.size());
  tempmat._Entry.resize(_Entry.size());
  for(ull i: _JEntry)
    tempmat._IEntry[i + 1]++;
  std::vector <ull> temp = tempmat._IEntry;
  for(ull i = 0; i < tempmat.matrow; ++i)
    tempmat._IEntry[i + 1] += tempmat._IEntry[i];
  for(ull i = 0; i < matrow; ++i)
  {
    for(ull j = _IEntry[i]; j < _IEntry[i + 1]; ++j)
    {
      ull k = _JEntry[j] + 1;
      tempmat._JEntry[tempmat._IEntry[k] - temp[k]] = i;
      tempmat._Entry[tempmat._IEntry[k] - temp[k]] = _Entry[j];
      temp[k]--;
    }
  }
  temp.clear();
  return tempmat;
}

template <class U>
void CSR<U>::resize(const ull& m, const ull& n)
{
  _IEntry.resize(m + 1,_IEntry.back());
  std::vector <ull> IA(m + 1, 0), JA;
  std::vector <U> A;
  for(ull i = 0; i < m; ++i)
  {
    IA[i + 1] += IA[i];
    for(ull j = _IEntry[i]; j < _IEntry[i + 1]; ++j)
    {
      if(_JEntry[j] < n)
      {
        JA.push_back(_JEntry[j]);
        A.push_back(_Entry[j]);
        IA[i + 1]++;
      }
    }
  }
  _JEntry = JA;
  _IEntry = IA;
  _Entry = A;
  JA.clear();
  IA.clear();
  A.clear();
  matrow = m;
  matcol = n;
}


template <class U>
void CSR<U>::resize(const ull& n)
{
  _IEntry.resize(n + 1,_IEntry.back());
  std::vector <ull> IA(n + 1, 0), JA;
  std::vector <U> A;
  for(ull i = 0; i < n; ++i)
  {
    IA[i + 1] += IA[i];
    for(ull j = _IEntry[i]; j < _IEntry[i + 1]; ++j)
    {
      if(_JEntry[j] < n)
      {
        JA.push_back(_JEntry[j]);
        A.push_back(_Entry[j]);
        IA[i + 1]++;
      }
    }
  }
  _JEntry = JA;
  _IEntry = IA;
  _Entry = A;
  JA.clear();
  IA.clear();
  A.clear();
  matrow = n;
  matcol = n;
}


template <class U>
ull CSR<U>::row_size() const
{
  return matrow;
}

template <class U>
ull CSR<U>::col_size() const
{
  return matcol;
}

#endif
