#ifndef AliAODRedCov_H
#define AliAODRedCov_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Reduced Cov Matrix
//     Author: fca
//-------------------------------------------------------------------------

#include <Rtypes.h>
#include <TMath.h>

template <Int_t N> class AliAODRedCov {


   //
   //  Class containing reduced cov matrix, see example here for a track
   //
   //       X          Y          Z         Px        Py        Pz
   //
   // X  fDiag[ 0]  
   //
   // Y  fOdia[ 0]  fDiag[ 1]
   //
   // Z  fOdia[ 1]  fOdia[ 2]  fDiag[ 2]
   //
   // Px fOdia[ 3]  fOdia[ 4]  fOdia[ 5]  fDiag[ 3]
   //
   // Py fOdia[ 6]  fOdia[ 7]  fOdia[ 8]  fOdia[ 9]  fDiag[ 4]
   //
   // Pz fOdia[10]  fOdia[11]  fOdia[12]  fOdia[13]  fOdia[14]  fDiag[ 5]
   //

 public:
  AliAODRedCov() {
    for(Int_t i=0; i<N; i++)         fDiag[i]       = 0.;
    for(Int_t i=0; i<N*(N-1)/2; i++) fODia[i]       = 0.;
  }
  virtual ~AliAODRedCov() {}
  template <class T> void GetCovMatrix(T *cmat) const;
  template <class T> void SetCovMatrix(T *cmat);
  
 private:
  Double32_t   fDiag[N];         // Diagonal elements
  Double32_t   fODia[N*(N-1)/2]; // [-1, 1,8] 8 bit precision for off diagonal elements
  
  ClassDef(AliAODRedCov,1)

 };

//Cint craps out here, we protect this part
#if !defined(__CINT__) && !defined(__MAKECINT__)

//#define DEBUG

//______________________________________________________________________________
template <Int_t N> template <class T> inline void AliAODRedCov<N>::GetCovMatrix(T *cmat) const
{
  //
  // Returns the external cov matrix
  //

  for(Int_t i=0; i<N; ++i) {
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
      cmat[i*(i+1)/2+j] = (fDiag[j] >= 0. && fDiag[i] >= 0.) ? fODia[(i-1)*i/2+j]*fDiag[j]*fDiag[i]: -999.;
#ifdef DEBUG
      printf("cmat[%2d] = fODia[%2d]*fDiag[%2d]*fDiag[%2d] = %f\n",
	     i*(i+1)/2+j,(i-1)*i/2+j,j,i,cmat[i*(i+1)/2+j]);
#endif
    }

    // Diagonal elements
    cmat[i*(i+1)/2+i] = (fDiag[i] >= 0.) ? fDiag[i]*fDiag[i] : -999.;
#ifdef DEBUG
    printf("cmat[%2d] = fDiag[%2d]*fDiag[%2d] = %f\n",
	   i*(i+1)/2+i,i,i,cmat[i*(i+1)/2+i]);
#endif
  }
}


//______________________________________________________________________________
template <Int_t N> template <class T> inline void AliAODRedCov<N>::SetCovMatrix(T *cmat)
{
  //
  // Sets the external cov matrix
  //

  if(cmat) {
    
#ifdef DEBUG
    for (Int_t i=0; i<(N*(N+1))/2; i++) {
      printf("cmat[%d] = %f\n", i, cmat[i]);
    }
#endif
    
    // Diagonal elements first
    for(Int_t i=0; i<N; ++i) {
      fDiag[i] = (cmat[i*(i+1)/2+i] >= 0.) ? TMath::Sqrt(cmat[i*(i+1)/2+i]) : -999.;
#ifdef DEBUG
	printf("fDiag[%2d] = TMath::Sqrt(cmat[%2d]) = %f\n",
	       i,i*(i+1)/2+i, fDiag[i]);
#endif
  }
  
  // ... then the ones off diagonal
  for(Int_t i=0; i<N; ++i) 
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
      fODia[(i-1)*i/2+j] = (fDiag[i] > 0. && fDiag[j] > 0.) ? cmat[i*(i+1)/2+j]/(fDiag[j]*fDiag[i]) : 0.;
      // check for division by zero (due to diagonal element of 0) and for fDiag != -999. (due to negative input diagonal element).
      if (fODia[(i-1)*i/2+j]>1.) { // check upper boundary
#ifdef DEBUG
	printf("out of bounds: %f\n", fODia[(i-1)*i/2+j]);
#endif
	fODia[(i-1)*i/2+j] = 1.;
      }
      if (fODia[(i-1)*i/2+j]<-1.) { // check lower boundary
#ifdef DEBUG
	printf("out of bounds: %f\n", fODia[(i-1)*i/2+j]);
#endif
	fODia[(i-1)*i/2+j] = -1.; 
      }
#ifdef DEBUG
	printf("fODia[%2d] = cmat[%2d]/(fDiag[%2d]*fDiag[%2d]) = %f\n",
	       (i-1)*i/2+j,i*(i+1)/2+j,j,i,fODia[(i-1)*i/2+j]); 
#endif
    }
  } else {
    for(Int_t i=0; i< N; ++i) fDiag[i]=-999.;
    for(Int_t i=0; i< N*(N-1)/2; ++i) fODia[i]=0.;
  }

  return;
}

#undef DEBUG

#endif
#endif
