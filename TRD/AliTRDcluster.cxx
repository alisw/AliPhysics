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

#include "AliLog.h"
#include "AliTRDcluster.h"

ClassImp(AliTRDcluster)

//___________________________________________________________________________
AliTRDcluster::AliTRDcluster() 
  :AliCluster() 
  ,fPadCol(0)
  ,fPadRow(0)
  ,fPadTime(0)
  ,fLocalTimeBin(0)
  ,fNPads(0)
  ,fClusterMasking(0)
  ,fDetector(0)
  ,fQ(0)
  ,fCenter(0)
{ 
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 7; i++) {
    fSignals[i] = 0;
  }

}

//___________________________________________________________________________
AliTRDcluster::AliTRDcluster(Int_t det, Float_t q
                           , Float_t *pos, Float_t *sig
                           , Int_t *tracks, Char_t npads, Short_t *signals
                           , UChar_t col, UChar_t row, UChar_t time
                           , Char_t timebin, Float_t center, UShort_t volid)
  :AliCluster(volid,pos[0],pos[1],pos[2],sig[0],sig[1],0.0,0x0) 
  ,fPadCol(col)
  ,fPadRow(row)
  ,fPadTime(time)
  ,fLocalTimeBin(timebin)
  ,fNPads(npads)
  ,fClusterMasking(0)
  ,fDetector(det)
  ,fQ(q)
  ,fCenter(center)
{ 
  //
  // Constructor
  //

  for (Int_t i = 0; i < 7; i++) {
    fSignals[i] = signals[i];
  }

  if (tracks) {
    AddTrackIndex(tracks);
  }

}

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(const AliTRDcluster &c)
  :AliCluster(c)
  ,fPadCol(c.fPadCol)
  ,fPadRow(c.fPadRow)
  ,fPadTime(c.fPadTime)
  ,fLocalTimeBin(c.fLocalTimeBin)
  ,fNPads(c.fNPads)
  ,fClusterMasking(c.fClusterMasking)
  ,fDetector(c.fDetector)
  ,fQ(c.fQ)
  ,fCenter(c.fCenter)
{
  //
  // Copy constructor 
  //

  SetBit(kInChamber, c.IsInChamber());
  SetLabel(c.GetLabel(0),0);
  SetLabel(c.GetLabel(1),1);
  SetLabel(c.GetLabel(2),2);

  SetY(c.GetY());
  SetZ(c.GetZ());
  SetSigmaY2(c.GetSigmaY2());
  SetSigmaZ2(c.GetSigmaZ2());  

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
void AliTRDcluster::Clear(Option_t *)
{
  //
  // Reset all member to the default value
  //
  fPadCol=0;
  fPadRow=0;
  fPadTime=0;
  fLocalTimeBin=0;
  fNPads=0;
  fClusterMasking=0;
  fDetector=0;
  for (Int_t i=0; i < 7; i++) fSignals[i]=0;
  fQ = 0;
  fCenter = 0;
  for (Int_t i = 0; i < 3; i++) SetLabel(0,i);
  SetX(0);
  SetY(0);
  SetZ(0);
  SetSigmaY2(0);
  SetSigmaZ2(0);
  SetVolumeId(0);
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


//_____________________________________________________________________________
Bool_t AliTRDcluster::IsEqual(const TObject *o) const
{
  //
  // Compare relevant information of this cluster with another one
  //
  
  const AliTRDcluster *inCluster = dynamic_cast<const AliTRDcluster*>(o);
  if (!o || !inCluster) return kFALSE;

  if ( AliCluster::GetX() != inCluster->GetX() ) return kFALSE;
  if ( AliCluster::GetY() != inCluster->GetY() ) return kFALSE;
  if ( AliCluster::GetZ() != inCluster->GetZ() ) return kFALSE;
  if ( fQ != inCluster->fQ ) return kFALSE;
  if ( fDetector != inCluster->fDetector ) return kFALSE;
  if ( fPadCol != inCluster->fPadCol ) return kFALSE;
  if ( fPadRow != inCluster->fPadRow ) return kFALSE;
  if ( fPadTime != inCluster->fPadTime ) return kFALSE;
  if ( fClusterMasking != inCluster->fClusterMasking ) return kFALSE;
  if ( IsInChamber() != inCluster->IsInChamber() ) return kFALSE;
  if ( IsShared() != inCluster->IsShared() ) return kFALSE;
  if ( IsUsed() != inCluster->IsUsed() ) return kFALSE;
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliTRDcluster::Print(Option_t *o) const
{
  AliInfo(Form("Det[%3d] LTrC[%7.2f %7.2f %7.2f] Q[%f] Stat[in(%c) use(%c) sh(%c)]", 
    fDetector, GetX(), GetY(), GetZ(), fQ, 
    IsInChamber() ? 'y' : 'n', IsUsed() ? 'y' : 'n', IsShared() ? 'y' : 'n'));

  if(strcmp(o, "a")!=0) return;
  AliInfo(Form("LChC[c(%3d) r(%2d) t(%2d)] t-t0[%2d] Npad[%d] cen[%5.3f] mask[%d]", fPadCol, fPadRow, fPadTime, fLocalTimeBin, fNPads, fCenter, fClusterMasking)); 
  AliInfo(Form("Signals[%3d %3d %3d %3d %3d %3d %3d]", fSignals[0], fSignals[1], fSignals[2], fSignals[3], fSignals[4], fSignals[5], fSignals[6]));
}


//_____________________________________________________________________________
void AliTRDcluster::SetPadMaskedPosition(UChar_t position)
{
  //
  // store the pad corruption position code
  // 
  // Code: 1 = left cluster
  //       2 = middle cluster;
  //       4 = right cluster
  //
  for(Int_t ipos = 0; ipos < 3; ipos++)
    if(TESTBIT(position, ipos))
      SETBIT(fClusterMasking, ipos);
}

//_____________________________________________________________________________
void AliTRDcluster::SetPadMaskedStatus(UChar_t status)
{
  //
  // store the status of the corrupted pad
  //
  // Code: 2 = noisy
  //       4 = Bridged Left
  //       8 = Bridged Right
  //      32 = Not Connected
  for(Int_t ipos = 0; ipos < 5; ipos++)
    if(TESTBIT(status, ipos))
      SETBIT(fClusterMasking, ipos + 3);
}
