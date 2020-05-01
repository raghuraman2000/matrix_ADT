#include "matrices.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <cmath>
typedef long long ll;
typedef unsigned long long ull;
typedef long double lld;

template <class V>
std::istream& operator >> (std::istream& is, matrix<V>& arg1)
{
  for(ull i=0;i<arg1.row_size();++i)
    for(ull j=0;j<arg1.col_size();++j)
      is >> arg1(i,j);
  return is;
}

template <class V>
std::ostream& operator << (std::ostream& os, const matrix<V>& arg1)
{
  os << '\n';
  for(ull i=0;i<arg1.row_size();++i)
  {
    for(ull j=0;j<arg1.col_size();++j)
      os << std::left << std::setw(15) <<  arg1(i,j);
    os << '\n' << '\n';
  }
  return os;
}

template <class V>
std::istream& operator >> (std::istream& is, COO<V>& arg1)
{
  for(ull i=0;i<arg1.row_size();++i)
    for(ull j=0;j<arg1.col_size();++j)
      is >> arg1(i,j);
  return is;
}

template <class V>
std::ostream& operator << (std::ostream& os, const COO<V>& arg1)
{
  os << '\n';
  for(ull i=0;i<arg1.row_size();++i)
  {
    for(ull j=0;j<arg1.col_size();++j)
      os << std::left << std::setw(15) <<  arg1(i,j);
    os << '\n' << '\n';
  }
  return os;
}

template <class V>
std::ostream& operator << (std::ostream& os, const CSR<V>& arg1)
{
  os << '\n';
  for(ull i=0;i<arg1.row_size();++i)
  {
    for(ull j=0;j<arg1.col_size();++j)
      os << std::left << std::setw(15) <<  arg1(i,j);
    os << '\n' << '\n';
  }
  return os;
}

int main()
{
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(NULL);
  std::cout.tie(NULL);

  matrix <ll> a(2,2);
  std::cout << "\nEnter a 2 by 2 matrix:\n\n";
  std::cin >> a;
  std::cout << a;
  matrix <ll> b(2,1);
  std::cout << "\nEnter a 2 by 1 matrix:\n\n";
  std::cin >> b;
  std::cout << b;
  matrix <ll> c;
  c = a * b;
  std::cout << "\nThis is the product of the top two matrices:\n\n";
  std::cout << c;
  c *= 5;
  std::cout << "\nMultiplying every entry by 5:\n\n";
  std::cout << c;
  matrix <ll> d = c.transpose()*6;
  d(0,1) += 7;
  std::cout << "\nTaking 6 times the transpose and adding 7 to the (0,1) entry:\n\n";
  std::cout << d;
  sqr_matrix <lld> f(2),g;
  std::cout << "\nEnter a square matrix of size 2:\n\n";
  std::cin >> f;
  std::cout << f;
  g = f.exp(10);
  std::cout << "\nThis is the above matrix raised to the power 10:\n\n";
  std::cout << g;
  std::vector <lld> h(2);
  h[0] = h[1] = 1;
  h = g.eq_solve(h); // gx = h
  std::cout << "\nSuppose that the above matrix is the coefficient matrix of a system of two linear equations.The RHS is say 1,1. As follows:\n\n";
  std::cout << "89x + 55y = 1\n\n55x + 34y = 1";
  std::cout << "\nThen, the solution is as follows:\n\n";
  std::cout << h[0] << ' ' << h[1] << '\n';
  std::cout << "\nFinding the determinant of the previous square matrix:\n\n";
  std::cout << g.det() << '\n';
  std::cout << "\nFinding the inverse of the previous square matrix:\n\n";
  std::cout << g.inv();
  std::cout << "\nMultiplying the matrix by its negative:\n\n";
  g *= -g;
  std::cout << g;
  sqr_matrix <lld> p(3);
  std::cout << "\nEnter a square matrix of size 3:\n\n";
  std::cin >> p;
  std::cout << p;
  std::cout << "\nIt's inverse:\n\n";
  std::cout << p.inv();
  matrix <lld> q, r;
  p.QR(q,r);
  std::cout << "\nIt's decomposition into an orthogonal matrix and an upper triangular matrix:\n\n";
  std::cout << q << r;
  matrix <lld> s(3);
  std::cout << "\nEnter a 3 by 3 matrix:\n\n";
  std::cin >> s;
  std:: cout << s;
  std::cout << "\nFinding it's column rank:\n\n";
  std::cout << s.col_rank() << '\n';
  std::cout << "\nEnter 4 points in the 2 dimensional Cartesian plane:\n\n";
  std::vector < std::pair<lld,lld> > func(4);
  for(ll i=0;i<4;++i)
  {
    std::cin >> func[i].first >> func[i].second;
    std::cout << func[i].first << ' ' << func[i].second << '\n';
  }
  std::vector <lld> fit = matrix<lld>::Poly_fit(func,2);
  std::cout << "\nFinding the best quadratic fit :\n\n";
  for(ull i = 0; i < 3; ++i)
    std::cout << fit[i] << " x^" << i << ((i == 2)? "\n" : " + ");

  COO <ll> t;
  std::cout << "\nTrying to make a 3 by 3 COO matrix and changing (1,1) and (0,0) to be 5:\n\n";
  t.resize(3,3);
  std::cout << t.row_size() << ' ' << t.col_size() << '\n';
  t(1,1) = 5;
  t(0,0) = t(1,1);
  std::cout << t;
  COO <ll> u(3,3);
  std::cout << "\nEnter a 3 by 3 COO matrix:\n\n";
  std::cin >> u;
  std::cout << u;
  t *= 2;
  t += u;
  std::cout << (u - t);
  std::cout << (-t*3).transpose();
  COO <ll> v(3,4);
  v(0,0) = v(1,1) = v(2,2) = 3;
  v(0,1) = 2;
  v(1,3) = 7;
  std::cout << "\nOutputing a 3 by 4 matrix and the product of it's transpose and itself:\n\n";
  std::cout << v << v.transpose()*v << '\n';
  COO<ll> aa(a), bb(b);
  std::cout << "\nChecking the multiplication of two COO matrices equal to the first two matrices in the demo:\n\n";
  std::cout << aa * bb;

  CSR <ll> ww(3);
  std::cout << "\nTrying to create a 3 by 3 CSR matrix:\n\n";
  std::cout << ww.row_size() << ' ' << ww.col_size() << '\n';
  ww.insert(0,1,7);
  ww.insert(0,0,ww(0,1));
  std::cout << ww;
  CSR <ll> xx(3,3);
  xx.insert(0,2,8);
  xx.insert(2,1,6);
  xx.insert(2,0,5);
  std::cout << "\nAdding another 3 by 3 CSR matrix to the above matrix:\n\n";
  std::cout << ww + xx;
  ww -= xx;
  std::cout << "\nOutputing a matrix and it's transpose:\n\n";
  std::cout << ww << ww.transpose();
  ww.resize(1,2);
  ww.resize(3);
  std::cout << "\nResizing the above matrix to a 1 by 2 and then a 3 by 3:\n\n";
  std::cout << ww;
  std::cout << "\nOutputing -9 times the 3 by 3 identity:\n\n";
  std::cout << CSR<ll>::iden(3) * -9;
  CSR <lld> qq(q), rr(r);
  qq *= -rr;
  std::cout << "\nChecking the multiplication of two COO matrices equal to the QR decomposition in the demo with the second one being negative:\n\n";
  std::cout << qq;
  return 0;
}
