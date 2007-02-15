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


#endif
