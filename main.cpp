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
  std::cin >> a;
  std::cout << a;
  matrix <ll> b(2,1);
  std::cin >> b;
  std::cout << b;
  matrix <ll> c;
  c = a * b;
  std::cout << c;
  c *= 5;
  std::cout << c;
  matrix <ll> d = c.transpose()*6;
  d(0,1) += 7;
  std::cout << d;
  sqr_matrix <lld> f(2),g;
  std::cin >> f;
  g = f.exp(10);
  std::cout << g;
  std::vector <lld> h(2);
  h[0] = h[1] = 1;
  h = g.eq_solve(h);
  std::cout << h[0] << ' ' << h[1] << '\n';
  std::cout << g.det() << '\n';
  std::cout << g.inv();
  g *= -g;
  std::cout << g;
  sqr_matrix <lld> p(3);
  std::cin >> p;
  std::cout << p.inv();
  matrix <lld> q, r;
  p.QR(q,r);
  std::cout << q << r;
  matrix <lld> s(3);
  std::cin >> s;
  std::cout << s.col_rank() << '\n';
  std::vector < std::pair<lld,lld> > func(4);
  for(ll i=0;i<4;++i)
    std::cin >> func[i].first >> func[i].second;
  std::vector <lld> fit = matrix<lld>::Poly_fit(func,2);
  for(auto y: fit)
    std::cout << y << '\n';
  std::cout << '\n';

  COO <ll> t;
  t.resize(3,3);
  std::cout << t.row_size() << ' ' << t.col_size() << '\n';
  t(1,1) = 5;
  t(0,0) = t(1,1);
  std::cout << t;
  COO <ll> u(3,3);
  std::cin >> u;
  t *= 2;
  t += u;
  std::cout << (u - t);
  std::cout << (-t*3).transpose();
  COO <ll> v(3,4);
  v(0,0) = v(1,1) = v(2,2) = 3;
  v(0,1) = 2;
  v(1,3) = 7;
  std::cout << v << v.transpose()*v << '\n';
  COO<ll> aa(a), bb(b);
  std::cout << aa * bb;

  CSR <ll> ww(3);
  std::cout << ww.row_size() << ' ' << ww.col_size() << '\n';
  ww.insert(0,1,7);
  ww.insert(0,0,ww(0,1));
  std::cout << ww;
  CSR <ll> xx(3,3);
  xx.insert(0,2,8);
  xx.insert(2,1,6);
  xx.insert(2,0,5);
  std::cout << ww + xx;
  ww -= xx;
  std::cout << ww << ww.transpose();
  ww.resize(1,2);
  ww.resize(3);
  std::cout << ww;
  std::cout << CSR<ll>::iden(3) * -9;
  CSR <lld> qq(q), rr(r);
  qq *= -rr;
  std::cout << qq;
  return 0;
}
