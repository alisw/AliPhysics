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
//  TRD cluster                                                              //
//                                                                           //
/////////////////////////////////////////////////////////////////////////////// 

#include "AliTRDcluster.h"
#include "AliTRDrecPoint.h"

ClassImp(AliTRDcluster)

//___________________________________________________________________________
AliTRDcluster::AliTRDcluster() 
  :AliCluster() 
  ,fDetector(0)
  ,fX(0)
  ,fTimeBin(0)
  ,fQ(0)
  ,fNPads(0)
  ,fCenter(0)
{ 
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 7; i++) {
    fSignals[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(const AliTRDrecPoint &p)
  :AliCluster()
  ,fDetector(p.GetDetector())
  ,fX(0)
  ,fTimeBin(p.GetLocalTimeBin())
  ,fQ(p.GetEnergy())
  ,fNPads(0)
  ,fCenter(0)
{
  //
  // Constructor from AliTRDrecPoint
  //

  fTracks[0] = p.GetTrackIndex(0);
  fTracks[1] = p.GetTrackIndex(1);
  fTracks[2] = p.GetTrackIndex(2);
  fY         = p.GetY();
  fZ         = p.GetZ();

  //fSigmaY2   = p.GetSigmaY2();
  //fSigmaZ2   = p.GetSigmaZ2();  
  // Why is this ????
  fSigmaY2   = 0.2;
  fSigmaZ2   = 5.0;  

}

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(const AliTRDcluster &c)
  :AliCluster()
  ,fDetector(c.fDetector)
  ,fX(c.fX)
  ,fTimeBin(c.fTimeBin)
  ,fQ(c.fQ)
  ,fNPads(c.fNPads)
  ,fCenter(c.fCenter)
{
  //
  // Copy constructor 
  //

  fTracks[0] = c.GetLabel(0);
  fTracks[1] = c.GetLabel(1);
  fTracks[2] = c.GetLabel(2);

  fY         = c.GetY();
  fZ         = c.GetZ();
  fSigmaY2   = c.GetSigmaY2();
  fSigmaZ2   = c.GetSigmaZ2();  

  for (Int_t i = 0; i < 7; i++) {
    fSignals[i] = c.fSignals[i];
  }

}

//_____________________________________________________________________________
void AliTRDcluster::AddTrackIndex(Int_t *track)
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
  Int_t  entries[kSize][2];

  Int_t  i = 0;
  Int_t  j = 0;
  Int_t  k = 0;
  Int_t  index;
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
          entries[j][1] = entries[j][1] + 1;
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
             (entries[i][0] >  entries[i+1][0]))) {
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
  for (i = 0; i < 3; i++) {
    SetLabel(entries[i][0],i);
  }

  return;

}          

//_____________________________________________________________________________
void AliTRDcluster::SetSignals(Short_t *signals)
{
  //
  // Write signals in the cluster
  //

  for (Int_t i = 0; i < 7; i++) {
    fSignals[i] = signals[i];
  }

}

//_____________________________________________________________________________
Float_t AliTRDcluster::GetSumS() const
{
  //
  // Returns the total charge from a not unfolded cluster
  //

  Float_t sum = 0.0;
  for (Int_t i = 0; i < 7; i++) {
    sum += fSignals[i];
  }

  return sum;

}
