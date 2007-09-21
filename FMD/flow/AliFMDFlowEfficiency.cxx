/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
#include <flow/AliFMDFlowEfficiency.h>
#include <cmath>
#include <iostream>
template <typename T>
T sign(const T& x, const T& y) 
{
  return (y >= 0 ? fabs(x) : -fabs(x));
}

//____________________________________________________________________
Double_t
AliFMDFlowEfficiency::LnGamma(Double_t z)
{
  if (z <= 0) return 0;
  Double_t c[] = { 2.5066282746310005, 76.18009172947146, 
		   -86.50532032941677, 24.01409824083091,  
		   -1.231739572450155, 0.1208650973866179e-2, 
		   -0.5395239384953e-5};
  Double_t x = z;
  Double_t y = x;
  Double_t t = x + 5.5;
  t        = (x+.5 * log(t) - t);
  Double_t s = 1.000000000190015;
  for (UInt_t i = 1; i < 7; i++) { 
    y  += 1;
    s  += c[i] / y;
  }
  Double_t v = t + log(c[0]*s/x);
  return v;
}

//____________________________________________________________________
Double_t
AliFMDFlowEfficiency::BetaCf(Double_t x, Double_t a, Double_t b)
{
  Int_t    itmax = 500;   // Maximum # of iterations
  Double_t eps   = 3e-14; // Precision. 
  Double_t fpmin = 1e-30;
  
  Double_t qab = a + b;
  Double_t qap = a + 1;
  Double_t qam = a - 1;
  Double_t c   = 1;
  Double_t d   = 1 - qab * x  / qap;
  if (fabs(d) < fpmin) d = fpmin;
  d          = 1/d;
  Double_t h   = d;
  Int_t m;
  for (m = 1; m <= itmax; m++) { 
    Int_t    m2 =  m * 2;
    Double_t aa =  m * (b - m) * x / ((qam + m2) * (a + m2));
    d         =  1 + aa * d;
    if (fabs(d) < fpmin) d = fpmin;
    c         =  1 + aa / c;
    if (fabs(c) < fpmin) d = fpmin;
    d         =  1 / d;
    h         *= d * c;
    aa        =  -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d         = 1 + aa * d;
    if (fabs(d) < fpmin) d = fpmin;
    c         =  1 + aa / c;
    if (fabs(c) < fpmin) d = fpmin;
    d         =  1 / d;
    Double_t dd = d * c;
    h         *= dd;
    if (fabs(dd - 1) <= eps) break;
  }
  if (m > itmax) { 
    std::cerr << "a or b too big, or max # iterations too small"
	      << std::endl;
  }
  return h;
}

  

//____________________________________________________________________
Double_t
AliFMDFlowEfficiency::IBetaI(Double_t a, Double_t b, Double_t x)
{
  Double_t bt = 0;
  if (x < 0 || x > 1) { 
    std::cerr << "x not in [0,1]" << std::endl;
    return 0;
  }
  if (x == 0 || x == 1) bt = 0;
  else 
    bt = (exp(LnGamma(a + b) - LnGamma(a) - LnGamma(b) 
	      + a * log(x) + b * log(1 - x)));
  if (x < (a + 1) / (a + b + 2)) 
    return bt * BetaCf(x, a, b) / a;
  return 1 - bt * BetaCf(1-x, b, a)/ b;
}

//____________________________________________________________________
Double_t
AliFMDFlowEfficiency::BetaAB(Double_t a, Double_t b, Int_t k, Int_t n)
{
  if (a == b) return 0;
  Int_t c1 = k + 1;
  Int_t c2 = n - k + 1;
  return IBetaI(c1, c2, b) - IBetaI(c1, c2, a);
}

//____________________________________________________________________
Double_t
AliFMDFlowEfficiency::SearchUpper(Double_t low, Int_t k, Int_t n, Double_t c)
{
  Double_t integral = BetaAB(low, 1, k, n);
  if (integral == c) return 1; // lucky -- this is the solution
  if (integral <  c) return -1;// no solution exists
  Double_t tooHigh = 1;         // upper edge estimate
  Double_t tooLow  = low;
  Double_t test    = 0;

  // Use a bracket-and-bisect search.  Loop 20 times, to end up with a
  // root guaranteed accurate to better than 0.1%.  
  for (UInt_t i = 0; i < 20; i++) { 
    test     = 0.5 * (tooLow + tooHigh);
    integral = BetaAB(low, test, k, n);
    if (integral > c) tooHigh = test;
    else              tooLow  = test;
  }
  return test;
}
//____________________________________________________________________
Double_t
AliFMDFlowEfficiency::SearchLower(Double_t high, Int_t k, Int_t n, Double_t c)
{
  Double_t integral = BetaAB(0, high, k, n);
  if (integral == c) return 0; // lucky -- this is the solution
  if (integral <  c) return -1;// no solution exists
  Double_t tooLow  = 0;         // lower edge estimate
  Double_t tooHigh = high;
  Double_t test     = 0;

  // Use a bracket-and-bisect search.  Loop 20 times, to end up with a
  // root guaranteed accurate to better than 0.1%.  
  for (UInt_t i = 0; i < 20; i++) { 
    test     = 0.5 * (tooHigh + tooLow);
    integral = BetaAB(test, high, k, n);
    if (integral > c) tooLow  = test;
    else              tooHigh = test;
  }
  return test;
}

//____________________________________________________________________
Double_t 
AliFMDFlowEfficiency::Brent(Double_t ax, Double_t bx, Double_t cx, 
			    Double_t& xmin, Int_t k, Int_t n, Double_t c, 
			    Double_t tol)
{
  const Int_t    kItMax = 100;
  const Double_t kGold  = 0.3819660;
  const Double_t kEps   = 1e-10;
  
  Int_t iter;
  Double_t a  = (ax < cx ? ax : cx);
  Double_t b  = (ax > cx ? ax : cx);
  Double_t d  = 0;
  Double_t e  = 0;
  Double_t x  = bx;
  Double_t w  = bx;
  Double_t v  = bx;
  Double_t fw = Length(x, k, n, c);
  Double_t fv = fw;
  Double_t fx = fv;
  
  for (iter = 1; iter <= kItMax; iter++) { 
    Double_t xm   = .5 * (a + b);
    Double_t tol1 = tol * fabs(x) + kEps;
    Double_t tol2 = 2 * tol1;
    if (fabs(x - xm) < (tol2 - .5 * (b-a))) {
      xmin = x;
      return fx;
    }
    if (fabs(e) > tol1) { 
      Double_t r = (x - w) * (fx - fv);
      Double_t q = (x - v) * (fx - fw);
      Double_t p = (x - v) * q - (x - w) * r;
      q        = 2 * (q - r);
      if (q > 0) p *= -1;
      q        = fabs(q);
      Double_t t = e;
      e        = d;
      if (fabs(p) >= fabs(.5 * q * t) || 
	  p <= q * (a - x)            || 
	  p >= q * (b - x)) {
	e = (x > xm ? a - x : b - x);
	d = kGold * e;
      }
      else {
	d        = p / q;
	Double_t u = x + d;
	if (u - a < tol2 || b - u < tol2) 
	  d = sign(tol1, xm - x);
      }
    }
    else { 
      e = (x >= xm ? a - x : b - x);
      d = kGold * e;
    }
    Double_t u  = (fabs(d) >= tol1 ? x + d : x + sign(tol1, d));
    Double_t fu = Length(u, k, n, c);
    if (fu <= fx) { 
      if (u >= x) a = x;
      else        b = x;
      v  = w;
      w  = x;
      x  = u;
      fv = fw;
      fw = fx;
      fx = fu;
    }
    else { 
      if (u < x) a = u;
      else       b = u;
      if (fu <= fw || w == x) { 
	v  = w;
	w  = u;
	fv = fw;
	fw = fu;
      }
      else if (fu <= fv || v == x || v == w) { 
	v  = u;
	fv = fu;
      }
    }
  }
  std::cerr << "Brett: too many iterations" << std::endl;
  xmin = x;
  return fx;
}
  
//____________________________________________________________________
Double_t 
AliFMDFlowEfficiency::Length(Double_t l, Int_t k, Int_t n, Double_t c)
{
  Double_t h = SearchUpper(l, k, n, c);
  if (h == -1) return 2;
  return (h-l);
}

//____________________________________________________________________
Double_t 
AliFMDFlowEfficiency::Interval(Int_t k, Int_t n, 
			       Double_t& low, Double_t& high, 
			       Double_t c)
{
  if (n == 0) { 
    low  = 0;
    high = 1;
    return .5;
  }
  
  Double_t e = Double_t(k) / n;
  Double_t l, h;
  
  if (k == 0) { 
    l = 0;
    h = SearchUpper(l, k, n, c);
  }
  else if (k == n) { 
    h = 1;
    l = SearchLower(h, k, n, c);
  }
  else { 
    Brent(0, .5, 1, l, n, k, 1e-9, c);
    h = l + Length(l, n, k, c);
  }
  low  = l;
  high = h;
  return e;
}

//____________________________________________________________________
//
// EOF
//

    
  
  
