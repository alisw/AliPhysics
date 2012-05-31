/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.2  2007/10/12 13:36:27  cvetan
Coding convention fixes from Stefan

Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 

////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorFFT class
////
//// Wrapper class to perform Fast Fourier Transformations.
//// The code is based on the Gnu Scientific Library. 
//// See documentation of gsl for further details.
//// 
//// Author: Stefan Kniege, IKF, Frankfurt
////       
////
/////////////////////////////////////////////////////////////////////////


#include "AliTPCMonitorFFT.h"


ClassImp(AliTPCMonitorFFT)
  
//__________________________________________________________________
  AliTPCMonitorFFT::AliTPCMonitorFFT()
{
  // Constuctor
}
 
//__________________________________________________________________
AliTPCMonitorFFT::~AliTPCMonitorFFT()
{
  // Destructor
}


//__________________________________________________________________
Int_t AliTPCMonitorFFT::ComplexRadix2ForwardWrap(Double_t* data, Int_t stride, size_t n)
{

  // Wrapper function for forward fourier transformation
  Int_t direction = 1;
  Int_t ret = ComplexRadix2TransformWrap(data,stride,n,direction);
  return ret ;
} 

//__________________________________________________________________
Int_t AliTPCMonitorFFT::ComplexRadix2BackwardWrap(Double_t* data,Int_t stride, size_t n)
{
  // Wrapper function for backward  fourier transformation
  Int_t direction = -1;
  Int_t ret = ComplexRadix2TransformWrap(data,stride,n,direction);
  return ret ;
} 

//__________________________________________________________________
Int_t AliTPCMonitorFFT::ComplexRadix2InverseWrap(Double_t* data,Int_t stride, size_t n)
{
  // Wrapper function for inverse fourier transformation
  Int_t direction = -1;
  Int_t ret = ComplexRadix2TransformWrap(data,stride,n,direction);
  
  if (ret)
    {
      return ret;
    }
  
  const Double_t knorm = 1.0 / n;
  size_t i;
  for (i = 0; i < n; i++)
    {
      REAL(data,stride,i) *= knorm;
      IMAG(data,stride,i) *= knorm;
    }
  
  return ret;
} 

//__________________________________________________________________
Int_t AliTPCMonitorFFT::ComplexRadix2TransformWrap(Double_t* data,   Int_t stride, size_t n, Int_t sign)
{
  // Wrapper function for fourier transformation depending on sign  ( forward 1;backward -1 )
  Int_t result ;
  size_t dual;
  size_t bit; 
  size_t logn = 0;
  Int_t status;
  
  if (n == 1) /* identity operation */
    {
      return 0 ;
    }
  /* make sure that n is a power of 2 */
  
  result = FFTBinaryLogn(n) ;
  
  if (result == -1) 
    {
      cout << " not a power of 2 " << endl;
    } 
  else 
    {
      logn = result ;
    }

  /* apply fft recursion */

  dual = n / 2;

  for (bit = 0; bit < logn; bit++)
    {
      Double_t wreal = 1.0;
      Double_t wimag = 0.0;

      const Double_t ktheta = 2.0 * ((Int_t) sign) * M_PI / ((Double_t) (2 * dual));

      const Double_t ks = sin (ktheta);
      const Double_t kt = sin (ktheta / 2.0);
      const Double_t ks2 = 2.0 * kt * kt;

      size_t a, b;

      for (b = 0; b < dual; b++)
	{
	  for (a = 0; a < n; a+= 2 * dual)
	    {
	      const size_t ki = b + a;
	      const size_t kj = b + a + dual;
              
	      const Double_t kt1real = REAL(data,stride,ki) + REAL(data,stride,kj);
	      const Double_t kt1imag = IMAG(data,stride,ki) + IMAG(data,stride,kj);
	      const Double_t kt2real = REAL(data,stride,ki) - REAL(data,stride,kj);
	      const Double_t kt2imag = IMAG(data,stride,ki) - IMAG(data,stride,kj);

	      REAL(data,stride,ki) = kt1real;
	      IMAG(data,stride,ki) = kt1imag;
	      REAL(data,stride,kj) = wreal*kt2real - wimag * kt2imag;
	      IMAG(data,stride,kj) = wreal*kt2imag + wimag * kt2real;
	    }

	  /* trignometric recurrence for w-> exp(i ktheta) w */

	  {
	    const Double_t ktmpreal = wreal - ks * wimag - ks2 * wreal;
	    const Double_t ktmpimag = wimag + ks * wreal - ks2 * wimag;
	    wreal = ktmpreal;
	    wimag = ktmpimag;
	  }
	}
      dual /= 2;
    }

  /* bit reverse the ordering of output data for decimation in
     frequency algorithm */
  
  //status = FUNCTION(fft_complex,bitreverse_order)(data, stride, n, logn) ;
  status  = ComplexBitReverseOrderWrap(data, stride, n, logn) ;
  return 0;
}

//__________________________________________________________________
Int_t AliTPCMonitorFFT::ComplexBitReverseOrderWrap(Double_t* data, Int_t stride, size_t n, Int_t /*logn*/) const
{
  // Wrapper function from gnu scientific library
  /* This is the Goldrader bit-reversal algorithm */

  size_t i;
  size_t j = 0;
  
//   logn = 0 ; /* not needed for this algorithm */
  
  for (i = 0; i < n - 1; i++)
    {
      size_t k = n / 2 ;
      
      if (i < j)
        {
          const Double_t ktmpreal = REAL(data,stride,i);
          const Double_t ktmpimag = IMAG(data,stride,i);
          REAL(data,stride,i) = REAL(data,stride,j);
          IMAG(data,stride,i) = IMAG(data,stride,j);
          REAL(data,stride,j) = ktmpreal;
          IMAG(data,stride,j) = ktmpimag;
        }
      
      while (k <= j) 
        {
          j = j - k ;
          k = k / 2 ;
        }
      
      j += k ;
    }

  return 0;

}

//__________________________________________________________________
Int_t AliTPCMonitorFFT::FFTBinaryLogn(size_t n) const
{

  // Return log on base 2
  size_t ntest ;
  size_t binarylogn = 0 ;
  size_t k = 1;
  
  while (k < n)
    {
      k *= 2;
      binarylogn++;
    }
  
  ntest = (1 << binarylogn) ;
  
  if (n != ntest )       
    {
      return -1 ; /* n is not a power of 2 */
    } 
  
  return binarylogn;
}

