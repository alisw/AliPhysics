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

  AliTRDcluster::AliTRDcluster() : AliCluster() { 
  //
  // default constructor
  //
  fQ=0; 
  fTimeBin=0; 
  fDetector=0; 
  fNPads=0; 
  for (Int_t i = 0;i<7; i++) fSignals[i]=0;
}
//_____________________________________________________________________________
  AliTRDcluster::AliTRDcluster(const AliTRDrecPoint &p):AliCluster()
{
  //
  // Constructor from AliTRDrecPoint
  //

  fDetector   = p.GetDetector();
  fTimeBin    = p.GetLocalTimeBin();

  fTracks[0]  = p.GetTrackIndex(0);
  fTracks[1]  = p.GetTrackIndex(1);
  fTracks[2]  = p.GetTrackIndex(2);

  fQ          = p.GetEnergy();

  fY          = p.GetY();
  fZ          = p.GetZ();
  fSigmaY2    = p.GetSigmaY2();
  fSigmaZ2    = p.GetSigmaZ2();  

  fSigmaY2    = 0.2;
  fSigmaZ2    = 5.;  
  fNPads      =0;
  fCenter     = 0;
}

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(const AliTRDcluster &c):AliCluster()
{
  //
  // Copy constructor 
  //

  fTracks[0]  = c.GetLabel(0);
  fTracks[1]  = c.GetLabel(1);
  fTracks[2]  = c.GetLabel(2);

  fY          = c.GetY();
  fZ          = c.GetZ();
  fSigmaY2    = c.GetSigmaY2();
  fSigmaZ2    = c.GetSigmaZ2();  

  fDetector   = c.GetDetector();
  fTimeBin    = c.GetLocalTimeBin();
  fQ          = c.GetQ();
  fNPads      = c.fNPads;
  fCenter     = c.fCenter;
  for (Int_t i=0;i<7;i++) fSignals[i] = c.fSignals[i];
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

  Int_t entries[kSize][2], i, j, index;

  Bool_t indexAdded;

  for (i=0; i<kSize; i++) {
    entries[i][0]=-1;
    entries[i][1]=0;
  }                                 

  for (Int_t k=0; k<kSize; k++) {
    index=track[k];
    indexAdded=kFALSE; 
    j=0;
    if (index >= 0) {
      while ( (!indexAdded) && ( j < kSize ) ) {
        if ((entries[j][0]==index) || (entries[j][1]==0)) {
          entries[j][0]=index;
          entries[j][1]=entries[j][1]+1;
          indexAdded=kTRUE;
        }
        j++;
      }
    }
  }             

  // sort by number of appearances and index value
  Int_t swap=1, tmp0, tmp1;
  while ( swap > 0) {
    swap=0;
    for(i=0; i<(kSize-1); i++) {
      if ((entries[i][0] >= 0) && (entries[i+1][0] >= 0)) {
        if ((entries[i][1] < entries[i+1][1]) ||
            ((entries[i][1] == entries[i+1][1]) &&
             (entries[i][0] > entries[i+1][0]))) {
               tmp0=entries[i][0];
               tmp1=entries[i][1];
               entries[i][0]=entries[i+1][0];
               entries[i][1]=entries[i+1][1];
               entries[i+1][0]=tmp0;
               entries[i+1][1]=tmp1;
               swap++;
        }
      }
    }
  }               

  // set track indexes
  for(i=0; i<3; i++) SetLabel(entries[i][0],i);

  return;

}          

void AliTRDcluster::SetSignals(Short_t*signals){
  //
  // write signals in the cluster
  //
  for (Int_t i = 0;i<7;i++) fSignals[i]=signals[i];
}

Float_t AliTRDcluster::GetSumS() const
{
  //
  // return total charge in non unfolded cluster
  //
  Float_t sum=0;
  for (Int_t i = 0;i<7;i++) sum+=fSignals[i];
  return sum;
}
Float_t AliTRDcluster::GetCenterS() const
{
  //
  //
  //
  Float_t sum=0;
  Float_t sum2=0;
  for (Int_t i = 0;i<7;i++) {    
    sum+=fSignals[i];
    sum2+=i*fSignals[i];
  }
  if (sum>0) return sum2/sum-2;
  return 0;

}
