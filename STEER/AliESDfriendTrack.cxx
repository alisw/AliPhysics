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
#include "TObjArray.h"
#include "AliKalmanTrack.h"

ClassImp(AliESDfriendTrack)

AliESDfriendTrack::AliESDfriendTrack(): 
TObject(), 
f1P(0), 
fPoints(0),
fCalibContainer(0),
fITStrack(0),
fTRDtrack(0)
{
  //
  // Default constructor
  //
  Int_t i;
  for (i=0; i<kMaxITScluster; i++) fITSindex[i]=-2;
  for (i=0; i<kMaxTPCcluster; i++) fTPCindex[i]=-2;
  for (i=0; i<kMaxTRDcluster; i++) fTRDindex[i]=-2;
}

AliESDfriendTrack::AliESDfriendTrack(const AliESDfriendTrack &t): 
TObject(t),
f1P(t.f1P),
fPoints(0),
fCalibContainer(0),
fITStrack(0),
fTRDtrack(0)
{
  //
  // Copy constructor
  //
  Int_t i;
  for (i=0; i<kMaxITScluster; i++) fITSindex[i]=t.fITSindex[i];
  for (i=0; i<kMaxTPCcluster; i++) fTPCindex[i]=t.fTPCindex[i];
  for (i=0; i<kMaxTRDcluster; i++) fTRDindex[i]=t.fTRDindex[i];
  if (t.fPoints) fPoints=new AliTrackPointArray(*t.fPoints);
  if (t.fCalibContainer) {
     fCalibContainer = new TObjArray(5);
     Int_t no=t.fCalibContainer->GetEntriesFast();
     for (i=0; i<no; i++) {
       TObject *o=t.fCalibContainer->At(i);
       fCalibContainer->AddLast(o->Clone());
     }  
  }
}

AliESDfriendTrack::~AliESDfriendTrack() {
  //
  // Simple destructor
  //
   delete fPoints;
   if (fCalibContainer) fCalibContainer->Delete();
   delete fCalibContainer;
   delete fITStrack;
   delete fTRDtrack;
}


void AliESDfriendTrack::AddCalibObject(TObject * calibObject){
  //
  // add calibration object to array -
  // track is owner of the objects in the container 
  //
  if (!fCalibContainer) fCalibContainer = new TObjArray(5);
  fCalibContainer->AddLast(calibObject);
}

TObject * AliESDfriendTrack::GetCalibObject(Int_t index){
  //
  //
  //
  if (!fCalibContainer) return 0;
  if (index>=fCalibContainer->GetEntriesFast()) return 0;
  return fCalibContainer->At(index);
}
