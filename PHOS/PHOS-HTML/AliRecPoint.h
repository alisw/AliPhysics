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
Revision 1.1  1999/12/17 09:01:14  fca
Y.Schutz new classes for reconstruction

*/

//-*-C++-*-
//_________________________________________________________________________
// Base Class of Cluster (empty cxx needed by Root)
//*-- Author : Yves Schutz  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TObjArray.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliRecPoint.h"

ClassImp(AliRecPoint)


//____________________________________________________________________________
AliRecPoint::AliRecPoint()
{
  // ctor  
  fAmp = 0.0 ; 
  
  fLocPos.SetXYZ(0., 0., 0.) ;
  fLocPosM = new TMatrix(3,3) ;
  fMaxDigit = 100 ; 
  fMulDigit = 0 ; 
  fDigitsList = new int[fMaxDigit]; ; 
  fMaxTrack = 5 ; 
  fMulTrack = 0 ; 
  fTracksList = new int[fMaxTrack]; ; 
}

//____________________________________________________________________________
AliRecPoint::~AliRecPoint()
{
  // dtor
  
  delete fLocPosM ; 
  if ( fDigitsList )     delete fDigitsList ; 
  if ( fTracksList )     delete fTracksList ;  
  
}
  
//____________________________________________________________________________
void AliRecPoint::AddDigit(AliDigitNew & digit)
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
  
  fDigitsList[fMulDigit++]=  (int) &digit  ; 
  fAmp += digit.GetAmp() ; 
}

//____________________________________________________________________________
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

//____________________________________________________________________________
void AliRecPoint::GetCovarianceMatrix(TMatrix & mat)
{
  // returns the covariant matrix for the local position
  
  mat = *fLocPosM ; 

}

//____________________________________________________________________________
void AliRecPoint::GetLocalPosition(TVector3 & pos)
{
  // returns the position of the cluster in the local reference system of the sub-detector

  pos = fLocPos;

 
}

//____________________________________________________________________________
void AliRecPoint::GetGlobalPosition(TVector3 & gpos, TMatrix & gmat)
{
  // returns the position of the cluster in the global reference system of ALICE
  // and the uncertainty on this position
  

  fGeom->GetGlobal(this, gpos, gmat) ;
 
}

//______________________________________________________________________________
void AliRecPoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliRecPoint.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fAmp;
      R__b.ReadArray(fDigitsList);
      R__b >> fGeom;
      fLocPos.Streamer(R__b);
      R__b >> fLocPosM;
      R__b >> fMulDigit;
      R__b >> fMulTrack;
      R__b.ReadArray(fTracksList);
   } else {
      R__b.WriteVersion(AliRecPoint::IsA());
      TObject::Streamer(R__b);
      R__b << fAmp;
      R__b.WriteArray(fDigitsList, fMaxDigit);
      R__b << fGeom;
      fLocPos.Streamer(R__b);
      R__b << fLocPosM;
      R__b << fMulDigit;
      R__b << fMulTrack;
      R__b.WriteArray(fTracksList, fMaxTrack);
   }
}
