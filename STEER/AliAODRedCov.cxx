/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

#include <TMath.h>

#include "AliAODRedCov.h"

templateClassImp(AliAODRedCov)

//______________________________________________________________________________
template <Int_t N> template <class T> void AliAODRedCov<N>::GetCovMatrix(T *cmat) const
{
  //
  // Returns the external cov matrix
  //

  for(Int_t i=0; i<N; ++i) {
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
      printf("cmat[%2d] = fODia[%2d]*fDiag[%2d]*fDiag[%2d];\n",
	     i*(i+1)/2+j,(i-1)*i/2+j,j,i);
      cmat[i*(i+1)/2+j] = fODia[(i-1)*i/2+j]*fDiag[j]*fDiag[i];}

    // Diagonal elements
    printf("cmat[%2d] = fDiag[%2d]*fDiag[%2d];\n",
	   i*(i+1)/2+i,i,i);
    cmat[i*(i+1)/2+i] = fDiag[i]*fDiag[i];
  }
}


//______________________________________________________________________________
template <Int_t N> template <class T> void AliAODRedCov<N>::SetCovMatrix(T *cmat)
{
  //
  // Sets the external cov matrix
  //

  if(cmat) {
    // Diagonal elements first
    for(Int_t i=0; i<N; ++i) {
      printf("fDiag[%2d] = TMath::Sqrt(cmat[%2d]);\n",
	     i,i*(i+1)/2+i);
      fDiag[i] = TMath::Sqrt(cmat[i*(i+1)/2+i]);}

  // ... then the ones off diagonal
  for(Int_t i=0; i<N; ++i) 
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
      printf("fODia[%2d] = cmat[%2d]/(fDiag[%2d]*fDiag[%2d]);\n",
	     (i-1)*i/2+j,i*(i+1)/2+j,j,i);
      fODia[(i-1)*i/2+j] = cmat[i*(i+1)/2+j]/(fDiag[j]*fDiag[i]);
    }
  } else {
    for(Int_t i=0; i< N; ++i) fDiag[i]=-999.;
    for(Int_t i=0; i< N*(N-1)/2; ++i) fODia[i]=0.;
  }
}
