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

//-------------------------------------------------------------------------
//               Implementation of the AliESDfriendTrack class
//  This class keeps complementary to the AliESDtrack information 
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------
#include "AliTrackPointArray.h"
#include "AliESDfriendTrack.h"
#include "AliESD.h"

ClassImp(AliESDfriendTrack)

  AliESDfriendTrack::AliESDfriendTrack(): TObject(), f1P(0), fPoints(0) {
  //
  // Default constructor
  //
  Int_t i;
  for (i=0; i<AliESDtrack::kMaxITScluster; i++) fITSindex[i]=-2;
  for (i=0; i<AliESDtrack::kMaxTPCcluster; i++) fTPCindex[i]=-2;
  for (i=0; i<AliESDtrack::kMaxTRDcluster; i++) fTRDindex[i]=-2;
}

AliESDfriendTrack::AliESDfriendTrack(const AliESDfriendTrack &t): 
TObject(t),
f1P(t.f1P),
fPoints(0)
{
  //
  // Copy constructor
  //
  Int_t i;
  for (i=0; i<AliESDtrack::kMaxITScluster; i++) fITSindex[i]=t.fITSindex[i];
  for (i=0; i<AliESDtrack::kMaxTPCcluster; i++) fTPCindex[i]=t.fTPCindex[i];
  for (i=0; i<AliESDtrack::kMaxTRDcluster; i++) fTRDindex[i]=t.fTRDindex[i];
  if (t.fPoints) fPoints=new AliTrackPointArray(*t.fPoints);
}

AliESDfriendTrack::AliESDfriendTrack(const AliESDtrack &t): 
TObject(t),
f1P(t.Get1P()),
fPoints(0) 
{
  //
  // Extracts the complementary info from the ESD track
  //
  t.GetITSclusters(fITSindex); 
  t.GetTPCclusters(fTPCindex); 
  t.GetTRDclusters(fTRDindex); 
  const AliTrackPointArray *points=t.GetTrackPointArray();
  if (points) fPoints=new AliTrackPointArray(*points);
}

AliESDfriendTrack::~AliESDfriendTrack() {delete fPoints;}
