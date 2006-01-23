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

/* $Id$ */

//-*-C++-*-
//_________________________________________________________________________
// Base Class for reconstructed space points 
// usually coming from the clusterisation algorithms
// run on the digits
//
//*-- Author : Yves Schutz  SUBATECH 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliRecPoint.h"
#include "AliGeometry.h"
#include "AliDigitNew.h"

ClassImp(AliRecPoint)


//_______________________________________________________________________
AliRecPoint::AliRecPoint():
  fAmp(0),
  fGeom(0),
  fIndexInList(-1), // to be set when the point is already stored
  fLocPos(0,0,0),
  fLocPosM(0),
  fMaxDigit(100),
  fMulDigit(0),
  fMaxTrack(5),
  fMulTrack(0),
  fDigitsList(0),
  fTracksList(0)
{
  //
  // default ctor  
  //
}

//_______________________________________________________________________
AliRecPoint::AliRecPoint(const char * ):
  fAmp(0),
  fGeom(0),
  fIndexInList(-1), // to be set when the point is already stored
  fLocPos(0,0,0),
  fLocPosM(new TMatrixF(3,3)),
  fMaxDigit(100),
  fMulDigit(0),
  fMaxTrack(5),
  fMulTrack(0),
  fDigitsList(new int[fMaxDigit]),
  fTracksList(new int[fMaxTrack])
{
  //
  // Standard ctor  
  //
}

//_______________________________________________________________________
AliRecPoint::AliRecPoint(const AliRecPoint& recp):
  TObject(recp),
  fAmp(0),
  fGeom(0),
  fIndexInList(-1), // to be set when the point is already stored
  fLocPos(0,0,0),
  fLocPosM(0),
  fMaxDigit(100),
  fMulDigit(0),
  fMaxTrack(5),
  fMulTrack(0),
  fDigitsList(0),
  fTracksList(0)
{
  //
  // Copy constructor
  //
  recp.Copy(*this);
}

//_______________________________________________________________________
AliRecPoint::~AliRecPoint()
{
  // dtor
  
  delete fLocPosM ; 
  delete [] fDigitsList ; 
  delete [] fTracksList ;  
  
}
  
//_______________________________________________________________________
void AliRecPoint::AddDigit(AliDigitNew & digit)
{
  // adds a digit to the digits list
  // and accumulates the total amplitude and the multiplicity 
  
  
  if ( fMulDigit >= fMaxDigit ) { // increase the size of the list 
    int * tempo = new int[fMaxDigit*2]; 
    
    Int_t index ; 
    
    for ( index = 0 ; index < fMulDigit ; index++ )
      tempo[index] = fDigitsList[index] ; 
    
    delete [] fDigitsList ; 
    fDigitsList = tempo ; 
  }
  
  fDigitsList[fMulDigit] = digit.GetIndexInList()  ; 
  fMulDigit++ ; 
  fAmp += digit.GetAmp() ; 
}

//_______________________________________________________________________
// void AliRecPoint::AddTrack(AliTrack & track)
// {
//   // adds a digit to the digits list
//   // and accumulates the total amplitude and the multiplicity 


//   if ( fMulTrack >= fMaxTrack ) { // increase the size of the list 
//     int * tempo = new int[fMaxTrack*=2] ; 
//     Int_t index ; 
//     for ( index = 0 ; index < fMulTrack ; index++ )
//       tempo[index] = fTracksList[index] ; 
//     delete fTracksList ; 
//     fTracksList = tempo ; 
//   }

//   fTracksList[fMulTrack++]=  (int) &Track  ; 
// }

//_______________________________________________________________________
void AliRecPoint::Copy(TObject& recp) const
{
  //
  // Copy *this onto pts
  //
  // Copy all first
  if((TObject*)this != &recp) {
    ((TObject*) this)->Copy(recp);
    (dynamic_cast<AliRecPoint&>(recp)).fAmp = fAmp;
    (dynamic_cast<AliRecPoint&>(recp)).fGeom = fGeom;
    (dynamic_cast<AliRecPoint&>(recp)).fIndexInList = fIndexInList;
    (dynamic_cast<AliRecPoint&>(recp)).fLocPos = fLocPos;
    (dynamic_cast<AliRecPoint&>(recp)).fLocPosM = new TMatrixF(*fLocPosM);
    (dynamic_cast<AliRecPoint&>(recp)).fMaxDigit = fMaxDigit;
    (dynamic_cast<AliRecPoint&>(recp)).fMulDigit = fMulDigit;
    (dynamic_cast<AliRecPoint&>(recp)).fMaxTrack = fMaxTrack;
    (dynamic_cast<AliRecPoint&>(recp)).fMulTrack = fMulTrack;
    
    // Duplicate pointed objects
    (dynamic_cast<AliRecPoint&>(recp)).fDigitsList = new Int_t[fMulDigit];
    memcpy((dynamic_cast<AliRecPoint&>(recp)).fDigitsList,fDigitsList,fMulDigit*sizeof(Int_t));
    (dynamic_cast<AliRecPoint&>(recp)).fTracksList = new Int_t[fMulTrack];
    memcpy((dynamic_cast<AliRecPoint&>(recp)).fTracksList,fTracksList,fMulTrack*sizeof(Int_t));
  }
}

//_______________________________________________________________________
void AliRecPoint::GetCovarianceMatrix(TMatrixF & mat) const
{
  // returns the covariant matrix for the local position
  
  mat = *fLocPosM ; 

}

//____________________________________________________________________________
void AliRecPoint::GetLocalPosition(TVector3 & pos) const
{
  // returns the position of the cluster in the local reference system of the sub-detector

  pos = fLocPos;

 
}

//____________________________________________________________________________
void AliRecPoint::GetGlobalPosition(TVector3 & gpos, TMatrixF & gmat) const
{
  // returns the position of the cluster in the global reference system of ALICE
  // and the uncertainty on this position
  

  fGeom->GetGlobal(this, gpos, gmat) ;
 
}


