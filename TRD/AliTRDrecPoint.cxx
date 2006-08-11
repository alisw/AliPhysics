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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD reconstructed point                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRun.h"

#include "AliTRDgeometry.h"
#include "AliTRDrecPoint.h"

ClassImp(AliTRDrecPoint)

//_____________________________________________________________________________
AliTRDrecPoint::AliTRDrecPoint()
  :AliRecPoint()
  ,fDetector(0)
  ,fTimeBin(0)
  ,fUsed(0)
  ,fY(0)
  ,fZ(0)
  ,fSigmaY2(0)
  ,fSigmaZ2(0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 3; i++) {
    fTracks[i] = 0;
  }
  fGeom = AliTRDgeometry::GetGeometry();

}

//_____________________________________________________________________________
AliTRDrecPoint::AliTRDrecPoint(const char * opt)
  :AliRecPoint(opt)
  ,fDetector(0)
  ,fTimeBin(0)
  ,fUsed(0)
  ,fY(0)
  ,fZ(0)
  ,fSigmaY2(0)
  ,fSigmaZ2(0)
{
  //
  // Standard constructor
  //

  for (Int_t i = 0; i < 3; i++) {
    fTracks[i] = 0;
  }
  fGeom = AliTRDgeometry::GetGeometry();

}

//_____________________________________________________________________________
AliTRDrecPoint::~AliTRDrecPoint()
{
  //
  // AliTRDrecPoint destructor
  //

}

//_____________________________________________________________________________
void AliTRDrecPoint::AddDigit(Int_t digit)
{
  //
  // Adds the index of a digit to the digits list
  //

  // First resize the list 
  // (no clusters with more than 3 digits for the TRD
  if ((fMulDigit == 0) && 
      (fMaxDigit >= 5)) {
    fMaxDigit = 5;
    delete fDigitsList;
    fDigitsList = new int[fMaxDigit];
  }

  // Increase the size of the list if necessary
  if (fMulDigit >= fMaxDigit) { 
    fMaxDigit *= 2;
    Int_t *tempo = new Int_t[fMaxDigit]; 
    Int_t  index; 
    for (index = 0; index < fMulDigit; index++) {
      tempo[index] = fDigitsList[index]; 
    }
    delete fDigitsList; 
    fDigitsList = tempo; 
  }
  
  fDigitsList[fMulDigit++] = digit;

}

//_____________________________________________________________________________
void AliTRDrecPoint::SetLocalPosition(TVector3 & /*pos*/)
{
  //
  // Sets the position of the point in the local coordinate system
  // (row,col,time) and calculates the error matrix in the same
  // system.
  //

  AliFatal("Not implemented");

}

//_____________________________________________________________________________
void AliTRDrecPoint::SetTrackingYZ(Float_t /*sigmaY*/, Float_t /*sigmaZ*/)
{
  //
  // Sets the position of the point in the local coordinate system
  // of tracking sector
  //
  
  AliFatal("Not implemented");

}                                    

//_____________________________________________________________________________
void AliTRDrecPoint::AddTrackIndex(Int_t *track)
{
  //
  // Adds track index. Currently assumed that track is an array of
  // size 9, and up to 3 track indexes are stored in fTracks[3].
  // Indexes are sorted according to:
  //  1) index of max number of appearances is stored first
  //  2) if two or more indexes appear equal number of times, the lowest
  //     ones are stored first;
  //

  const Int_t kSize = 9;
  Int_t entries[kSize][2];

  Int_t i = 0;
  Int_t j = 0;
  Int_t k = 0;
  Int_t index;

  Bool_t indexAdded;

  for (i = 0; i < kSize; i++) {
    entries[i][0] = -1;
    entries[i][1] =  0;
  }

  for (k = 0; k < kSize; k++) {
    index      = track[k];
    indexAdded = kFALSE; 
    j = 0;
    if (index >= 0) {
      while ((!indexAdded) && (j < kSize)) {
        if ((entries[j][0] == index) || 
            (entries[j][1] ==     0)) {
          entries[j][0] = index;
          entries[j][1] = entries[j][1]+1;
          indexAdded    = kTRUE;
        }
        j++;
      }
    }
  }

  // Sort by number of appearances and index value
  Int_t swap = 1;
  Int_t tmp0;
  Int_t tmp1;
  while (swap > 0) {
    swap = 0;
    for (i = 0; i < (kSize - 1); i++) {
      if ((entries[i][0]   >= 0) && 
          (entries[i+1][0] >= 0)) {
        if ((entries[i][1] < entries[i+1][1]) ||
            ((entries[i][1] == entries[i+1][1]) &&
             (entries[i][0]  > entries[i+1][0]))) {
               tmp0            = entries[i][0];
               tmp1            = entries[i][1];
               entries[i][0]   = entries[i+1][0];
               entries[i][1]   = entries[i+1][1];
               entries[i+1][0] = tmp0;
               entries[i+1][1] = tmp1;
               swap++;
        }
      }
    }
  }

  // Set track indexes
  for(i = 0; i < 3; i++) {
    fTracks[i] = entries[i][0];
  }

  return;

}                    







