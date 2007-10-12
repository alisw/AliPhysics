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
Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 

////////////////////////////////////////////////////////////////////////
//
// AliTPCMonitorFFT class
//
// Wrapper class to perform Fast Fourier Transformations.
// The code is based on the Gnu Scientific Library. 
// See documentation of gsl for further details.
// 
// Author: Stefan Kniege, IKF, Frankfurt
//       
//
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
  
  const Double_t norm = 1.0 / n;
  size_t i;
  for (i = 0; i < n; i++)
    {
      REAL(data,stride,i) *= norm;
      IMAG(data,stride,i) *= norm;
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
      Double_t w_real = 1.0;
      Double_t w_imag = 0.0;

      const Double_t theta = 2.0 * ((Int_t) sign) * M_PI / ((Double_t) (2 * dual));

      const Double_t s = sin (theta);
      const Double_t t = sin (theta / 2.0);
      const Double_t s2 = 2.0 * t * t;

      size_t a, b;

      for (b = 0; b < dual; b++)
	{
	  for (a = 0; a < n; a+= 2 * dual)
	    {
	      const size_t i = b + a;
	      const size_t j = b + a + dual;
              
	      const Double_t t1_real = REAL(data,stride,i) + REAL(data,stride,j);
	      const Double_t t1_imag = IMAG(data,stride,i) + IMAG(data,stride,j);
	      const Double_t t2_real = REAL(data,stride,i) - REAL(data,stride,j);
	      const Double_t t2_imag = IMAG(data,stride,i) - IMAG(data,stride,j);

	      REAL(data,stride,i) = t1_real;
	      IMAG(data,stride,i) = t1_imag;
	      REAL(data,stride,j) = w_real*t2_real - w_imag * t2_imag;
	      IMAG(data,stride,j) = w_real*t2_imag + w_imag * t2_real;
	    }

	  /* trignometric recurrence for w-> exp(i theta) w */

	  {
	    const Double_t tmp_real = w_real - s * w_imag - s2 * w_real;
	    const Double_t tmp_imag = w_imag + s * w_real - s2 * w_imag;
	    w_real = tmp_real;
	    w_imag = tmp_imag;
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
Int_t AliTPCMonitorFFT::ComplexBitReverseOrderWrap(Double_t* data, Int_t stride, size_t n, Int_t logn)
{
  // Wrapper function from gnu scientific library
  /* This is the Goldrader bit-reversal algorithm */

  size_t i;
  size_t j = 0;
  
  logn = 0 ; /* not needed for this algorithm */
  
  for (i = 0; i < n - 1; i++)
    {
      size_t k = n / 2 ;
      
      if (i < j)
        {
          const Double_t tmp_real = REAL(data,stride,i);
          const Double_t tmp_imag = IMAG(data,stride,i);
          REAL(data,stride,i) = REAL(data,stride,j);
          IMAG(data,stride,i) = IMAG(data,stride,j);
          REAL(data,stride,j) = tmp_real;
          IMAG(data,stride,j) = tmp_imag;
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
Int_t AliTPCMonitorFFT::FFTBinaryLogn(size_t n)
{

  // Return log on base 2
  size_t ntest ;
  size_t binary_logn = 0 ;
  size_t k = 1;
  
  while (k < n)
    {
      k *= 2;
      binary_logn++;
    }
  
  ntest = (1 << binary_logn) ;
  
  if (n != ntest )       
    {
      return -1 ; /* n is not a power of 2 */
    } 
  
  return binary_logn;
}

