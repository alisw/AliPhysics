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
*/

//_________________________________________________________________________
// Class for PMD reconstructed space points 
// usually coming from the clusterisation algorithms
// run on the digits
//
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TObjArray.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPMDRecPoint.h"
#include "AliGeometry.h"
#include "AliDigitNew.h"

ClassImp(AliPMDRecPoint)


//____________________________________________________________________________
//____________________________________________________________________________
void AliPMDRecPoint::AddDigit(AliDigitNew & digit)
{
  // adds a digit to the digits list
  // and accumulates the total amplitude and the multiplicity 
  
  
  if ( fMulDigit >= fMaxDigit ) { // increase the size of the list 
    int * tempo = new ( int[fMaxDigit*=2] ) ; 
    
    Int_t index ; 
    
    for ( index = 0 ; index < fMulDigit ; index++ )
      tempo[index] = fDigitsList[index] ; 
    
    delete fDigitsList ; 
    fDigitsList = tempo ; 
  }
  
  fDigitsList[fMulDigit] = digit.GetIndexInList()  ; 
  fMulDigit++ ; 
  fAmp += digit.GetAmp() ; 
}

//____________________________________________________________________________
void AliPMDRecPoint::Copy(AliPMDRecPoint& recp) const
{
  //
  // Copy *this onto pts
  //
  // Copy all first
  if(this != &recp) {
    ((TObject*) this)->Copy((TObject&)recp);
    recp.fAmp = fAmp;
    recp.fGeom = fGeom;
    recp.fIndexInList = fIndexInList;
    recp.fLocPos = fLocPos;
    recp.fLocPosM = new TMatrix(*fLocPosM);
    recp.fMaxDigit = fMaxDigit;
    recp.fMulDigit = fMulDigit;
    recp.fMaxTrack = fMaxTrack;
    recp.fMulTrack = fMulTrack;
    
    // Duplicate pointed objects
    recp.fDigitsList = new Int_t[fMulDigit];
    memcpy(recp.fDigitsList,fDigitsList,fMulDigit*sizeof(Int_t));
    recp.fTracksList = new Int_t[fMulTrack];
    memcpy(recp.fTracksList,fTracksList,fMulTrack*sizeof(Int_t));
  }
}

//____________________________________________________________________________
void AliPMDRecPoint::GetCovarianceMatrix(TMatrix & mat)
{
  // returns the covariant matrix for the local position
  
  mat = *fLocPosM ; 

}

//____________________________________________________________________________
void AliPMDRecPoint::GetLocalPosition(TVector3 & pos) const
{
  // returns the position of the cluster in the local reference system of the sub-detector

  pos = fLocPos;

 
}

//____________________________________________________________________________
AliPMDRecPoint & AliPMDRecPoint::operator= (const AliPMDRecPoint &recp)
{
  recp.Copy(*this);
  return (*this);
}


//____________________________________________________________________________
void AliPMDRecPoint::GetGlobalPosition(TVector3 & gpos, TMatrix & gmat) const
{
  // returns the position of the cluster in the global reference system of ALICE
  // and the uncertainty on this position
  

  fGeom->GetGlobal(this, gpos, gmat) ;
 
}









