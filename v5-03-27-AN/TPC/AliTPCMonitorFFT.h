#ifndef ALITPCMONITORFFT_H
#define ALITPCMONITORFFT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorFFT class
////
//// Wrapper class to perform Fast Fourier Transformations
//// Code based on Gnu Scientific Library 
//// See documentation of gsl for further details
//// 
//// Author: Stefan Kniege, IKF, Frankfurt
////       
////
/////////////////////////////////////////////////////////////////////////



#include <iostream> 
#include "TNamed.h"
#include <math.h>
#define  REAL(a,stride,i) ((a)[2*(stride)*(i)])
#define  IMAG(a,stride,i) ((a)[2*(stride)*(i)+1])
using namespace std; 

class AliTPCMonitorFFT : public TNamed { 
 public:
    
    AliTPCMonitorFFT();
    ~AliTPCMonitorFFT(); 
    
    Int_t     ComplexRadix2ForwardWrap(  Double_t* data, Int_t stride, size_t n ); 
    Int_t     ComplexRadix2BackwardWrap( Double_t* data, Int_t stride, size_t n );  
    Int_t     ComplexRadix2InverseWrap(  Double_t* data, Int_t stride, size_t n ); 
    Int_t     ComplexRadix2TransformWrap(Double_t* data, Int_t stride, size_t n, Int_t sign ); 
    Int_t     ComplexBitReverseOrderWrap(Double_t* data, Int_t stride, size_t n, Int_t logn) const ; 
    Int_t     FFTBinaryLogn(size_t n) const ;

 private:
    
    ClassDef(AliTPCMonitorFFT,1);
};
#endif
 
