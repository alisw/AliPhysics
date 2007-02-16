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
  AliAODRedCov() {}
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
#ifdef DEBUG
      printf("cmat[%2d] = fODia[%2d]*fDiag[%2d]*fDiag[%2d];\n",
	     i*(i+1)/2+j,(i-1)*i/2+j,j,i);
#endif
      cmat[i*(i+1)/2+j] = fODia[(i-1)*i/2+j]*fDiag[j]*fDiag[i];}

    // Diagonal elements
#ifdef DEBUG
    printf("cmat[%2d] = fDiag[%2d]*fDiag[%2d];\n",
	   i*(i+1)/2+i,i,i);
#endif
    cmat[i*(i+1)/2+i] = fDiag[i]*fDiag[i];
  }
}


//______________________________________________________________________________
template <Int_t N> template <class T> inline void AliAODRedCov<N>::SetCovMatrix(T *cmat)
{
  //
  // Sets the external cov matrix
  //

  if(cmat) {
    // Diagonal elements first
    for(Int_t i=0; i<N; ++i) {
#ifdef DEBUG
      printf("fDiag[%2d] = TMath::Sqrt(cmat[%2d]);\n",
	     i,i*(i+1)/2+i);
#endif
      fDiag[i] = TMath::Sqrt(cmat[i*(i+1)/2+i]);}

  // ... then the ones off diagonal
  for(Int_t i=0; i<N; ++i) 
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
#ifdef DEBUG
      printf("fODia[%2d] = cmat[%2d]/(fDiag[%2d]*fDiag[%2d]);\n",
	     (i-1)*i/2+j,i*(i+1)/2+j,j,i);
#endif
      fODia[(i-1)*i/2+j] = cmat[i*(i+1)/2+j]/(fDiag[j]*fDiag[i]);
    }
  } else {
    for(Int_t i=0; i< N; ++i) fDiag[i]=-999.;
    for(Int_t i=0; i< N*(N-1)/2; ++i) fODia[i]=0.;
  }
}

#undef DEBUG

#endif
#endif
