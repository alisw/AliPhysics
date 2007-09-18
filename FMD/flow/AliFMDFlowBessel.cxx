/** @file 
    @brief Implementation of Bessel functions */
#include "flow/AliFMDFlowBessel.h"
#include <cmath>
#include <iostream>

/** Namespace for utility functions */
namespace 
{ 
  //__________________________________________________________________
  /** Utility function to compute 
      @f[ 
      envj_n(x) = \frac12 \log_{10}(6.28 n) - n \log_{10}(1.36 x / n)
      @f]
      @param n Order 
      @param x Argument 
      @return @f$ envj_n(x)@f$ */
  Double_t Envj(UInt_t n, Double_t x) { 
    return .5 * log10(6.28 * n) - n * log10(1.36 * x / n);
  }
  //__________________________________________________________________
  /** Utility function to do loop in Msta1 and Msta2 */ 
  UInt_t MstaLoop(Int_t n0, Double_t a0, Float_t mp) {
    Double_t f0 = Envj(n0, a0) - mp;
    Int_t    n1 = n0 + 5;
    Double_t f1 = Envj(n1, a0) - mp;
    Int_t    nn = 0;
    for (UInt_t i = 0; i < 20; i++) { 
      nn = n1 - (n1 - n0) / (1. - f0 / f1);
      if (fabs(nn - n1) < 1) break;
      n0 = n1;
      f0 = f1;
      n1 = nn;
      f1 = Envj(nn, a0) - mp;
    }
    return nn;
  }
  //__________________________________________________________________
  /** Determine the starting point for backward recurrence such that
      the magnitude of @f$J_n(x)@f$ at that point is about @f$
      10^{-mp}@f$  
      @param x   Argument to @f$ J_n(x)@f$ 
      @param mp  Value of magnitude @f$ mp@f$ 
      @return The starting point */
  UInt_t Msta1(Double_t x, UInt_t mp) { 
    Double_t a0 = fabs(x);
    Int_t    n0 = Int_t(1.1 * a0) + 1;
    return MstaLoop(n0, a0, mp);
  }
  //__________________________________________________________________
  /** Determine the starting point for backward recurrence such that
      all @f$ J_n(x)@f$ has @a mp significant digits.
      @param x   Argument to @f$ J_n(x)@f$ 
      @param n   Order of @f$ J_n@f$
      @param mp  Number of significant digits 
      @reutnr starting point */
  UInt_t Msta2(Double_t x, Int_t n, Int_t mp) { 
    Double_t a0  = fabs(x);
    Double_t hmp = 0.5 * mp;
    Double_t ejn = Envj(n, a0);
    Double_t obj = 0;
    Int_t    n0  = 0;
    if (ejn <= hmp) {
      obj = mp;
      n0  = Int_t(1.1 * a0) + 1; //Bug for x<0.1-vl, 2-8.2002
    }
    else { 
      obj = hmp + ejn;
      n0  = n;
    }
    return MstaLoop(n0, a0, obj) + 10;
  }
}
//____________________________________________________________________
void 
AliFMDFlowBessel::I01(Double_t x, Double_t& bi0, Double_t& di0, Double_t& bi1, Double_t& di1)
{
  Double_t       x2 =  x * x;
  if (x == 0) { 
    bi0 = 1;
    bi1 = 0;
    di0 = 0;
    di1 = .5;
    return;
  }
    
  if (x <= 18) { 
    bi0      = 1;
    Double_t r = 1;
    for (UInt_t k = 1; k <= 50; k++) { 
      r   =  .25 * r * x2 / (k * k);
      bi0 += r;
      if (fabs(r / bi0) < 1e-15) break;
    }
    bi1 = 1;
    r   = 1;
    for (UInt_t k = 1; k <= 50; k++) { 
      r   =  .25 * r * x2 / (k * (k + 1));
      bi1 += r;
    }
    bi1 *= .5 * x;
  }
  else { 
    const Double_t a[] = { 0.125,            0.0703125, 
			 0.0732421875,     0.11215209960938, 
			 0.22710800170898, 0.57250142097473, 
			 1.7277275025845,  6.0740420012735, 
			 24.380529699556,  110.01714026925,
			 551.33589612202,  3038.0905109224 };
    const Double_t b[] = { -0.375, 	     -0.1171875, 
			 -0.1025390625,    -0.14419555664063,
			 -0.2775764465332, -0.67659258842468,
			 -1.9935317337513, -6.8839142681099,
			 -27.248827311269, -121.59789187654,
			 -603.84407670507, -3302.2722944809 };
    UInt_t k0 = 12;
    if (x >= 35) k0 = 9;
    if (x >= 50) k0 = 7;
      
    Double_t ca = exp(x) / sqrt(2 * M_PI * x);
    Double_t xr = 1. / x;
    bi0       = 1;
      
    for (UInt_t k = 1; k <= k0; k++) {
      bi0 += a[k-1] * pow(xr, Int_t(k));
      bi1 += b[k-1] * pow(xr, Int_t(k));
    }
      
    bi0 *= ca;
    bi1 *= ca;
  }
  di0 = bi1;
  di1 = bi0 - bi1 / x;
}
//____________________________________________________________________
UInt_t 
AliFMDFlowBessel::Ihalf(Int_t n, Double_t x, Double_t* bi, Double_t* di) 
{
  typedef Double_t (*fun_t)(Double_t);
  Int_t   p = (n > 0 ? 1     : -1);
  fun_t s = (n > 0 ? &sinh : &cosh);
  fun_t c = (n > 0 ? &cosh : &sinh);
  
  // Temporary buffer 
  Double_t tbi[p*n/2+2];
  Double_t tdi[p*n/2+2];
  
  // Calculate I(-1/2) for n>0 or I(1/2) for n<0 
  Double_t bim = sqrt(2) * (*c)(x) / sqrt(M_PI * x);

  // We calculate one more than n, that is we calculate 
  //   I(1/2), I(3/2),...,I(n/2+1) 
  // because we need that for the negative orders 
  for (Int_t i = 1; i <= p * n + 2; i += 2) { 
    switch (i) { 
    case 1: 
      // Explicit formula for 1/2
      tbi[0] = sqrt(2 * x / M_PI) * (*s)(x) / x;
      break;
    case 3:
      // Explicit formula for 3/2
      tbi[1] = sqrt(2 * x / M_PI) * ((*c)(x) / x - (*s)(x) / (x * x));
      break;
    default:
      // Recursive formula for n/2 in terms of (n/2-2) and (n/2-1) 
      tbi[i/2] = tbi[i/2-2] - (i-2) * tbi[i/2-1] / x;
      break;
    }
  }
  // We calculate the first one dI(1/2) here, so that we can use it in
  // the loop. 
  tdi[0] = bim - tbi[0] / (2 * x);

  // Now, calculate dI(3/2), ..., dI(n/2)
  for (Int_t i = 3; i <= p * n; i += 2) {
    tdi[i/2] = tbi[i/2-p*1] - p * i * tbi[i/2] / (2 * x);
  }

  // Finally, store the output 
  for (Int_t i = 0; i < p * n / 2 + 1; i++) { 
    UInt_t j = (p < 0 ? p * n / 2 - i : i);
    bi[i]    = tbi[j];
    di[i]    = tdi[j];
  }
  return n / 2;
}
		
//____________________________________________________________________
UInt_t 
AliFMDFlowBessel::Iwhole(UInt_t in, Double_t x, Double_t* bi, Double_t* di) 
{
  UInt_t mn = in;
  if (x < 1e-100) { 
    for (UInt_t k = 0; k <= in; k++) bi[k] = di[k]  = 0;
    bi[0] = 1;
    di[1] = .5;
    return 1;
  }
  Double_t bi0, di0, bi1, di1;
  I01(x, bi0, di0, bi1, di1);
  bi[0] = bi0; di[0] = di0; bi[1] = bi1; di[1] = di1; 

  if (in <= 1) return 2;
    
  if (x > 40 and Int_t(in) < Int_t(.25 * x)) { 
    Double_t h0 = bi0;
    Double_t h1 = bi1;
    for (UInt_t k = 2; k <= in; k++) { 
      bi[k] = -2 * (k - 1) / x * h1 + h0;
      h0    = h1;
      h1    = bi[k];
    }
  }
  else { 
    UInt_t      m  = Msta1(x, 100);
    if (m < in) mn = m;
    else        m  = Msta2(x, in, 15);
      
    Double_t f0 = 0;
    Double_t f1 = 1e-100;
    Double_t f  = 0;
    for (Int_t k = m; k >= 0; k--) { 
      f  = 2 * (k + 1) * f1 / x + f0;
      if (k <= Int_t(mn)) bi[k] = f;
      f0 = f1;
      f1 = f;
    }
    Double_t s0 = bi0 / f;
    for (UInt_t k = 0; k <= mn; k++) 
      bi[k] *= s0;
  }
  for (UInt_t k = 2; k <= mn; k++) 
    di[k] =  bi[k - 1] - k / x * bi[k];
  return mn;
}

//____________________________________________________________________
UInt_t 
AliFMDFlowBessel::Inu(Double_t n1, Double_t n2, Double_t x, Double_t* bi, Double_t* di)
{
  UInt_t  in1 = UInt_t(fabs(n1));
  UInt_t  in2 = UInt_t(fabs(n2));
  UInt_t  nt  = UInt_t(n2 - n1 + 1);
  if (Int_t(2 * fabs(n1)) % 2 == 1) { // Half-integer argument 
    if (Int_t(2 * fabs(n2)) % 2 != 1) { 
      std::cout << "If n1 is half-integer (" << n1 << ") then n2 "
		<< "must also be half integer" << std::endl;
      return 0;
    }
    // std::cout << "Half-integers " << n1 << "," << n2 << std::endl;

    Int_t s1    = (n1 < 0 ? -1 : 1);
    Int_t s2    = (n2 < 0 ? -1 : 1);
    Int_t l1    = (s1 < 0 ? (2*in1+1)/2 : 0);
    Int_t l2    = (s1 > 0 ? in1 : 0);
    Double_t tbi[nt+1], tdi[nt+1];
    if (s1 < 0) { 
      Ihalf(2 * n1, x, tbi, tdi);
      for (Int_t i = 0; i < l1 && i < Int_t(nt); i++) { 
	bi[i] = tbi[i];
	di[i] = tdi[i];
      }
    }
    if (s2 > 0) { 
      // std::cout << "Evaluating Ihalf(" << 2 * n2 << "," << x << ",...)";
      Ihalf(2 * n2, x, tbi, tdi);
      for (Int_t i = l1; i <= 2 * Int_t(in2) && i < Int_t(nt); i++) { 
	UInt_t j = i + l2 - l1;
	bi[i] = tbi[j];
	di[i] = tdi[j];
      }
    }
    return nt;
  }
  if (Int_t(n1) != n1 || Int_t(n2) != n2) { 
    std::cerr << "n1 (" << n1 << "/" << in1 << ") and "
	      << "n2 (" << n2 << "/" << in2 << ") must be "
	      << "half-integer or integer" << std::endl;
    return 0; // Not integer!
  }

  UInt_t  n   = UInt_t(in1 > in2 ? in1 : in2);
  Double_t  tbi[n+1];
  Double_t  tdi[n+1];
  UInt_t  r   = Iwhole(n, x, tbi, tdi);
  if (r < n) 
    std::cerr << "Only got " << r << "/" << n << " values" 
	      << std::endl;
  for (UInt_t i = 0; i < nt; i++) { 
    UInt_t j = (i + n1 < 0 ? -n1-i : n1+i);
    bi[i]    = tbi[j];
    di[i]    = tdi[j];
  }
  return nt;
}

    
//____________________________________________________________________
Double_t 
AliFMDFlowBessel::I(Double_t n, Double_t x) 
{ 
  Double_t i[2], di[2];
  Inu(n, n, x, i, di);
  return i[0];
}
      
//____________________________________________________________________
Double_t 
AliFMDFlowBessel::DiffI(Double_t n, Double_t x) 
{ 
  Double_t i[2], di[2];
  Inu(n, n, x, i, di);
  return di[0];
}

//____________________________________________________________________
//
// EOF
//
