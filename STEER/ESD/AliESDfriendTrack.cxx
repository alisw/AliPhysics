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
#include "TClonesArray.h"
#include "AliKalmanTrack.h"

ClassImp(AliESDfriendTrack)

AliESDfriendTrack::AliESDfriendTrack(): 
TObject(), 
f1P(0), 
fnMaxITScluster(0),
fnMaxTPCcluster(0),
fnMaxTRDcluster(0),
fITSindex(0x0),
fTPCindex(0x0),
fTRDindex(0x0),
fPoints(0),
fCalibContainer(0),
fITStrack(0),
fTRDtrack(0),
fTPCOut(0),
fITSOut(0),
fTRDIn(0)
{
  //
  // Default constructor
  //
	//  Int_t i;
  //  fITSindex = new Int_t[fnMaxITScluster];
  //fTPCindex = new Int_t[fnMaxTPCcluster];
  //fTRDindex = new Int_t[fnMaxTRDcluster];
  //for (i=0; i<kMaxITScluster; i++) fITSindex[i]=-2;
  //for (i=0; i<kMaxTPCcluster; i++) fTPCindex[i]=-2;
  //for (i=0; i<kMaxTRDcluster; i++) fTRDindex[i]=-2;
  
  //fHmpPhotClus->SetOwner(kTRUE); 
  
}

AliESDfriendTrack::AliESDfriendTrack(const AliESDfriendTrack &t): 
TObject(t),
f1P(t.f1P),
fnMaxITScluster(t.fnMaxITScluster),
fnMaxTPCcluster(t.fnMaxTPCcluster),
fnMaxTRDcluster(t.fnMaxTRDcluster),
fITSindex(0x0),
fTPCindex(0x0),
fTRDindex(0x0),
fPoints(0),
fCalibContainer(0),
fITStrack(0),
fTRDtrack(0),
fTPCOut(0),
fITSOut(0),
fTRDIn(0)
{
  //
  // Copy constructor
  //
  AliDebug(2,"Calling copy constructor");

  Int_t i;
  if (fnMaxITScluster != 0){
	  fITSindex = new Int_t[fnMaxITScluster];
	  for (i=0; i<fnMaxITScluster; i++) fITSindex[i]=t.fITSindex[i];
  }
  if (fnMaxTPCcluster != 0){
	  fTPCindex = new Int_t[fnMaxTPCcluster];
	  for (i=0; i<fnMaxTPCcluster; i++) fTPCindex[i]=t.fTPCindex[i];
  }
  if (fnMaxTRDcluster != 0){
	  fTRDindex = new Int_t[fnMaxTRDcluster];
	  for (i=0; i<fnMaxTRDcluster; i++) fTRDindex[i]=t.fTRDindex[i]; 
  }
  AliDebug(2,Form("fnMaxITScluster = %d",fnMaxITScluster));
  AliDebug(2,Form("fnMaxTPCcluster = %d",fnMaxTPCcluster));
  AliDebug(2,Form("fnMaxTRDcluster = %d",fnMaxTRDcluster));
  if (t.fPoints) fPoints=new AliTrackPointArray(*t.fPoints);
  if (t.fCalibContainer) {
     fCalibContainer = new TObjArray(5);
     Int_t no=t.fCalibContainer->GetEntriesFast();
     for (i=0; i<no; i++) {
       TObject *o=t.fCalibContainer->At(i);
       if (o) fCalibContainer->AddLast(o->Clone());
     }  
  }

  if (t.fTPCOut) fTPCOut = new AliExternalTrackParam(*(t.fTPCOut));
  if (t.fITSOut) fITSOut = new AliExternalTrackParam(*(t.fITSOut));
  if (t.fTRDIn)  fTRDIn = new AliExternalTrackParam(*(t.fTRDIn));
  
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
   delete fTPCOut;
   delete fITSOut;
   delete fTRDIn;
   delete[] fITSindex;
   delete[] fTPCindex;
   delete[] fTRDindex;
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


void AliESDfriendTrack::SetTPCOut(const AliExternalTrackParam &param) {
  // 
  // backup TPC out track
  //
  delete fTPCOut;
  fTPCOut=new AliExternalTrackParam(param);
} 
void AliESDfriendTrack::SetITSOut(const AliExternalTrackParam &param) {
  //
  // backup ITS out track
  //
  delete fITSOut;
  fITSOut=new AliExternalTrackParam(param);
} 
void AliESDfriendTrack::SetTRDIn(const AliExternalTrackParam  &param)  {
  //
  // backup TRD in track
  //
  delete fTRDIn;
  fTRDIn=new AliExternalTrackParam(param);
} 

void AliESDfriendTrack::SetITSIndices(Int_t* indices, Int_t n){

	//
	// setting fITSindex
	// instantiating the pointer if still NULL
	//

	fnMaxITScluster = n;
	AliDebug(2,Form("fnMaxITScluster = %d",fnMaxITScluster));
	if (fITSindex == 0x0){
		fITSindex = new Int_t[fnMaxITScluster];
	}
	for (Int_t i = 0; i < fnMaxITScluster; i++){
		fITSindex[i] = indices[i];
	}
}

void AliESDfriendTrack::SetTPCIndices(Int_t* indices, Int_t n){

	//
	// setting fTPCindex
	// instantiating the pointer if still NULL
	//

	fnMaxTPCcluster = n;
	AliDebug(2,Form("fnMaxTPCcluster = %d",fnMaxTPCcluster));
	if (fTPCindex == 0x0){
		fTPCindex = new Int_t[fnMaxTPCcluster];
	}
	for (Int_t i = 0; i < fnMaxTPCcluster; i++){
		fTPCindex[i] = indices[i];
	}
}

void AliESDfriendTrack::SetTRDIndices(Int_t* indices, Int_t n){

	//
	// setting fTRDindex
	// instantiating the pointer if still NULL
	//

	fnMaxTRDcluster = n;
	AliDebug(2,Form("fnMaxTRDcluster = %d",fnMaxTRDcluster));
	if (fTRDindex == 0x0){
		fTRDindex = new Int_t[fnMaxTRDcluster];
	}
	for (Int_t i = 0; i < fnMaxTRDcluster; i++){
		fTRDindex[i] = indices[i];
	}
}

