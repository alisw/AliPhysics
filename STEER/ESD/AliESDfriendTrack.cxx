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
#include "AliPoolsSet.h"

ClassImp(AliESDfriendTrack)


AliPoolsSet* AliESDfriendTrack::fgPools = 0;

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
  AliPoolN* poolN = fgPools ? fgPools->GetPoolN() : 0;
  Int_t i;
  if (fnMaxITScluster != 0){
    fITSindex = poolN ? poolN->BookI(fnMaxITScluster,&fITSindex) : new Int_t[fnMaxITScluster];
    for (i=fnMaxITScluster; i--;) fITSindex[i]=t.fITSindex[i];
  }
  if (fnMaxTPCcluster != 0){
    fTPCindex = poolN ? poolN->BookI(fnMaxTPCcluster,&fTPCindex) : new Int_t[fnMaxTPCcluster];
    for (i=fnMaxTPCcluster; i--;) fTPCindex[i]=t.fTPCindex[i];
  }
  if (fnMaxTRDcluster != 0){
    fTRDindex = poolN ? poolN->BookI(fnMaxTRDcluster,&fTRDindex) : new Int_t[fnMaxTRDcluster];
    for (i=fnMaxTRDcluster; i--;) fTRDindex[i]=t.fTRDindex[i]; 
  }
  AliDebug(2,Form("fnMaxITScluster = %d",fnMaxITScluster));
  AliDebug(2,Form("fnMaxTPCcluster = %d",fnMaxTPCcluster));
  AliDebug(2,Form("fnMaxTRDcluster = %d",fnMaxTRDcluster));
  //
  if (t.fCalibContainer) { // RS!!!
     fCalibContainer = new TObjArray(5);
     fCalibContainer->SetOwner(kFALSE); //RS !!! objects are owned by pools
     Int_t no=t.fCalibContainer->GetEntriesFast();
     for (i=0; i<no; i++) {
       TObject *o=t.fCalibContainer->At(i);
       if (o) {
	 TObject* ocl = o->Clone();
	 fCalibContainer->AddLast(ocl);
	 // RS: Special treatment for CalibContainer objects: we don't know at which pool (if any) the 
	 // object was created, so we have to create the clone outside of any pool
	 if (ocl->InheritsFrom("AliExternalTrackParam")) ((AliExternalTrackParam*)ocl)->SetPoolID(-1);
	 else                                                                     ocl->SetUniqueID(0);
       }
     }  
  }
  //
  AliClonesPool* poolETP = fgPools ? fgPools->GetPoolExtTrPar() : 0;
  AliClonesPool* poolTPA = fgPools ? fgPools->GetPoolTrPoints() : 0;
  if (t.fPoints) {
    if (poolTPA) {fPoints = new(poolTPA->NextFreeSlot()) AliTrackPointArray(*t.fPoints); poolTPA->RegisterClone(fPoints);}
    else fPoints = new AliTrackPointArray(*t.fPoints);
  }
  if (poolETP) {
    if (t.fTPCOut) {fTPCOut = new(poolETP->NextFreeSlot()) AliExternalTrackParam(*(t.fTPCOut)); poolETP->RegisterClone(fTPCOut);}
    if (t.fITSOut) {fITSOut = new(poolETP->NextFreeSlot()) AliExternalTrackParam(*(t.fITSOut)); poolETP->RegisterClone(fITSOut);}
    if (t.fTRDIn)  {fTRDIn  = new(poolETP->NextFreeSlot()) AliExternalTrackParam(*(t.fTRDIn));  poolETP->RegisterClone(fTRDIn);}
  }
  else {
    if (t.fTPCOut) fTPCOut = new AliExternalTrackParam(*(t.fTPCOut));
    if (t.fITSOut) fITSOut = new AliExternalTrackParam(*(t.fITSOut));
    if (t.fTRDIn)  fTRDIn  = new AliExternalTrackParam(*(t.fTRDIn));
  }
}

//__________________________________________________________________
AliESDfriendTrack::AliESDfriendTrack(AliESDfriendTrack *t, Bool_t detach): 
TObject(*t),
f1P(t->f1P),
fnMaxITScluster(t->fnMaxITScluster),
fnMaxTPCcluster(t->fnMaxTPCcluster),
fnMaxTRDcluster(t->fnMaxTRDcluster),
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
  // Semi-shallow copy constructor: pointers of dynamic content is simply transfered, if detach is true the pointers of source are set to 0
  //
  AliDebug(2,"Calling semi-shallow copy constructor");
  if (fnMaxITScluster != 0) fITSindex = t->fITSindex;
  if (fnMaxTPCcluster != 0) fTPCindex = t->fTPCindex;
  if (fnMaxTRDcluster != 0) fTRDindex = t->fTRDindex;
  //
  AliDebug(2,Form("fnMaxITScluster = %d",fnMaxITScluster));
  AliDebug(2,Form("fnMaxTPCcluster = %d",fnMaxTPCcluster));
  AliDebug(2,Form("fnMaxTRDcluster = %d",fnMaxTRDcluster));
  //
  if (t->fCalibContainer) fCalibContainer = t->fCalibContainer;
  //
  if (t->fPoints) fPoints = t->fPoints;
  //
  if (t->fTPCOut) fTPCOut = t->fTPCOut;
  if (t->fITSOut) fITSOut = t->fITSOut;
  if (t->fTRDIn)  fTRDIn  = t->fTRDIn;
  //
  if (detach) {
    t->fITSindex = t->fTPCindex = t->fTRDindex = 0;
    t->fnMaxITScluster = t->fnMaxTPCcluster = t->fnMaxTRDcluster = 0;
    t->fCalibContainer = 0;
    t->fPoints = 0;
    t->fTPCOut = t->fITSOut = t->fTRDIn = 0;
  }
  //  
}


AliESDfriendTrack::~AliESDfriendTrack() 
{
  //
  // Simple destructor
  //
  Clear();
}

void AliESDfriendTrack::Clear(Option_t*) 
{
  // clear dynamical part
  if (fCalibContainer) {
    int n = fCalibContainer->GetEntriesFast();
    for (int i=0;i<=n;i++) {
      TObject* obj = fCalibContainer->RemoveAt(i);
      if (!obj) continue;
      // was it created in the pool?
      if (fgPools) {
	int poolID = (obj->InheritsFrom("AliExternalTrackParam")) ? ((AliExternalTrackParam*)obj)->GetPoolID() : obj->GetUniqueID()-1;
	if (poolID<0) delete obj; // not in the pool. those which are in the pool will be deleted on pool reset      
      }
      else delete obj;
    }
    delete fCalibContainer;
    fCalibContainer = 0;
  }
  //
  AliClonesPool* poolETP = fgPools ? fgPools->GetPoolExtTrPar() : 0;
  AliClonesPool* poolTPA = fgPools ? fgPools->GetPoolTrPoints() : 0;
  AliClonesPool* poolITS = fgPools ? fgPools->GetPoolTrITS()    : 0;
  AliClonesPool* poolTRD = fgPools ? fgPools->GetPoolTrTRD()    : 0;
  AliPoolN     * poolN   = fgPools ? fgPools->GetPoolN()        : 0;
  //
  if (!poolTPA) delete fPoints;
  else if (!poolTPA->IsReset()) poolTPA->MarkSlotFree(fPoints);
  fPoints = 0;
  //
  if (!poolETP) { // during reco these objects are stored in AliReconstruction pools
     delete fTPCOut;
     delete fITSOut;
     delete fTRDIn;
  }
  else if (!poolETP->IsReset()) {
    poolETP->MarkSlotFree(fTPCOut);
    poolETP->MarkSlotFree(fITSOut);
    poolETP->MarkSlotFree(fTRDIn);
  }
  fTPCOut = fITSOut = fTRDIn = 0;
  //
  if (!poolITS) delete fITStrack;
  else if (!poolITS->IsReset()) poolITS->MarkSlotFree(fITStrack);
  fITStrack = 0;
  //
  if (!poolTRD) delete fTRDtrack;
  else if (!poolTRD->IsReset()) poolTRD->MarkSlotFree(fTRDtrack); // no TRD tracks?
  fTRDtrack = 0;
  //
  //
  if (!poolETP) {
    delete[] fITSindex;
    delete[] fTPCindex;
    delete[] fTRDindex;
  }
  else if (!poolN->IsReset()) {
    if (fITSindex) poolN->FreeSlot(&fITSindex);
    if (fTPCindex) poolN->FreeSlot(&fTPCindex);
    if (fTRDindex) poolN->FreeSlot(&fTRDindex);
  }     
  fITSindex = fTPCindex = fTRDindex = 0;
  //
  TObject::Clear();
}

void AliESDfriendTrack::AddCalibObject(TObject * calibObject)
{
  // add calibration object to array -
  // track is owner of the objects in the container 
  //
  if (!fCalibContainer) fCalibContainer = new TObjArray(5); //RS!!! container does not own its objects
  fCalibContainer->AddLast(calibObject);
}

TObject * AliESDfriendTrack::GetCalibObject(Int_t index)
{
  // get object
  if (!fCalibContainer) return 0;
  if (index>=fCalibContainer->GetEntriesFast()) return 0;
  return fCalibContainer->At(index);
}


void AliESDfriendTrack::SetTPCOut(const AliExternalTrackParam &param) 
{
  // backup TPC out track
  //
  if (fTPCOut) {*fTPCOut = param; return;}
  //
  AliClonesPool* poolETP = fgPools ? fgPools->GetPoolExtTrPar() : 0;
  if (poolETP) {fTPCOut = new(poolETP->NextFreeSlot()) AliExternalTrackParam(param); poolETP->RegisterClone(fTPCOut);}
  else          fTPCOut = new AliExternalTrackParam(param);
  //
} 

void AliESDfriendTrack::SetITSOut(const AliExternalTrackParam &param) 
{
  // backup ITS out track
  //
  if (fITSOut) {*fITSOut = param; return;}
  //
  AliClonesPool* poolETP = fgPools ? fgPools->GetPoolExtTrPar() : 0;
  if (poolETP) {fITSOut = new(poolETP->NextFreeSlot()) AliExternalTrackParam(param); poolETP->RegisterClone(fITSOut);}
  else          fITSOut = new AliExternalTrackParam(param);
  //
} 

void AliESDfriendTrack::SetTRDIn(const AliExternalTrackParam  &param)  
{
  // backup TRD in track
  //
  if (fTRDIn) {*fTRDIn = param; return;}
  //
  AliClonesPool* poolETP = fgPools ? fgPools->GetPoolExtTrPar() : 0;
  if (poolETP) {fTRDIn = new(poolETP->NextFreeSlot()) AliExternalTrackParam(param); poolETP->RegisterClone(fTRDIn);}
  else          fTRDIn = new AliExternalTrackParam(param);
  //  
} 

void AliESDfriendTrack::SetITSIndices(Int_t* indices, Int_t n)
{
  //
  // setting fITSindex
  // instantiating the pointer if still NULL
  //  
  if ( !(fITSindex && (fnMaxITScluster==n)) ) {
    fnMaxITScluster = n;
    AliPoolN* poolN = fgPools ? fgPools->GetPoolN() : 0;
    //
    if (fITSindex) poolN ? poolN->FreeSlot(&fITSindex) : delete[] fITSindex; fITSindex = 0;
    if (poolN) fITSindex = poolN->BookI(fnMaxITScluster,&fITSindex);
    else       fITSindex = new Int_t[fnMaxITScluster];
  }
  AliDebug(2,Form("fnMaxITScluster = %d",n));
  for (Int_t i=fnMaxITScluster; i--;) fITSindex[i] = indices[i];
}

void AliESDfriendTrack::SetTPCIndices(Int_t* indices, Int_t n)
{
  //
  // setting fTPCindex
  // instantiating the pointer if still NULL
  //
  if ( !(fTPCindex && (fnMaxTPCcluster==n)) ) {
    fnMaxTPCcluster = n;
    AliPoolN* poolN = fgPools ? fgPools->GetPoolN() : 0;
    //
    if (fTPCindex) poolN ? poolN->FreeSlot(&fTPCindex) : delete[] fTPCindex; fTPCindex = 0;
    if (poolN) fTPCindex = poolN->BookI(fnMaxTPCcluster,&fTPCindex);
    else       fTPCindex = new Int_t[fnMaxTPCcluster];
  }
  AliDebug(2,Form("fnMaxTPCcluster = %d",n));
  for (Int_t i=fnMaxTPCcluster; i--;) fTPCindex[i] = indices[i];  
  //
}

void AliESDfriendTrack::SetTRDIndices(Int_t* indices, Int_t n)
{
  // setting fTRDindex
  // instantiating the pointer if still NULL
  //
  if ( !(fTRDindex && (fnMaxTRDcluster==n)) ) {
    fnMaxTRDcluster = n;
    AliPoolN* poolN = fgPools ? fgPools->GetPoolN() : 0;
    //
    if (fTRDindex) poolN ? poolN->FreeSlot(&fTRDindex) : delete[] fTRDindex; fTRDindex = 0;
    if (poolN) fTRDindex = poolN->BookI(fnMaxTRDcluster,&fTRDindex);
    else       fTRDindex = new Int_t[fnMaxTRDcluster];
  }
  AliDebug(2,Form("fnMaxTRDcluster = %d",n));
  for (Int_t i=fnMaxTRDcluster; i--;) fTRDindex[i] = indices[i];  
  //
}

