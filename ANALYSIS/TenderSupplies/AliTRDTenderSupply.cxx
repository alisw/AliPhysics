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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD tender: reapply pid on the fly                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TList.h>
#include <TObjString.h>
#include <TTree.h>
#include <TString.h>

#include <AliCDBEntry.h>
#include <AliCDBId.h>
#include <AliCDBManager.h>
#include <AliTRDCalDet.h>

#include <AliLog.h>
#include <TTree.h>
#include <TChain.h>
#include <AliPID.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliESDpid.h>
#include <AliESDtrack.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliTRDPIDResponse.h>
#include <AliTender.h>

#include "AliTRDTenderSupply.h"

ClassImp(AliTRDTenderSupply)

//_____________________________________________________
AliTRDTenderSupply::AliTRDTenderSupply() :
  AliTenderSupply(),
  fESD(NULL),
  fESDpid(NULL),
  fChamberGainOld(NULL),
  fChamberGainNew(NULL),
  fChamberVdriftOld(NULL),
  fChamberVdriftNew(NULL),
  fFileNNref(""),
  fPIDmethod(kNNpid),
  fPthreshold(0.8),
  fNBadChambers(0),
  fGainCorrection(kTRUE)
{
  //
  // default ctor
  //
  memset(fBadChamberID, 0, sizeof(Int_t) * kNChambers);
}

//_____________________________________________________
AliTRDTenderSupply::AliTRDTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fESD(NULL),
  fESDpid(NULL),
  fChamberGainOld(NULL),
  fChamberGainNew(NULL),
  fChamberVdriftOld(NULL),
  fChamberVdriftNew(NULL),
  fFileNNref(""),
  fPIDmethod(kNNpid),
  fPthreshold(0.8),
  fNBadChambers(0),
  fGainCorrection(kTRUE)
{
  //
  // named ctor
  //
  memset(fBadChamberID, 0, sizeof(Int_t) * kNChambers);
}

//_____________________________________________________
AliTRDTenderSupply::~AliTRDTenderSupply()
{
  //
  // dtor
  //
}

//_____________________________________________________
void AliTRDTenderSupply::Init()
{
  //
  // Initialise TRD tender
  //

  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  
  // 1DLQ PID implemented in the AliESD object
  fESDpid=fTender->GetESDhandler()->GetESDpid();
  if (!fESDpid) {
    fESDpid=new AliESDpid;
    fTender->GetESDhandler()->SetESDpid(fESDpid);
  }
  // Set Normalisation Factors
  if(mgr->GetMCtruthEventHandler()){
    // Assume MC
    //fESDpid->GetTRDResponse().SetGainNormalisationFactor(1.);
    SwitchOffGainCorrection();
  }
  else{
    // Assume Data
    //if(fPIDmethod == kNNpid) fPidRecal->SetGainScaleFactor(1.14);
    //fESDpid->GetTRDResponse().SetGainNormalisationFactor(1.14);
    SwitchOnGainCorrection();
  }
}

//_____________________________________________________
void AliTRDTenderSupply::ProcessEvent()
{
  //
  // Reapply pid information
  //
  if (fTender->RunChanged()){
    AliDebug(0, Form("AliTPCTenderSupply::ProcessEvent - Run Changed (%d)\n",fTender->GetRun()));
    if (fGainCorrection) SetChamberGain();
  }

  fESD = fTender->GetEvent();
  if (!fESD) return;
  Int_t ntracks=fESD->GetNumberOfTracks();
   
  //
  // recalculate PID probabilities
  //
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliESDtrack *track=fESD->GetTrack(itrack);
    // Recalculate likelihoods
    if(!(track->GetStatus() & AliESDtrack::kTRDout)) continue;
    if(fGainCorrection) ApplyGainCorrection(track);
    switch(fPIDmethod){
      case kNNpid:
        break;
      case k1DLQpid:
        fESDpid->MakeTRDPID(fESD->GetTrack(itrack));
        break;
      default:
        AliError("PID Method not implemented (yet)");
    }
  }
}

//_____________________________________________________
void AliTRDTenderSupply::SetChamberGain(){
  //
  // Load Chamber Gain factors into the Tender supply
  //
  
  //find previous entry from the UserInfo
  //   TTree *tree=((TChain*)fTender->GetInputData(0))->GetTree();
  AliAnalysisManager*mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisTaskSE *task = (AliAnalysisTaskSE*)mgr->GetTasks()->First();
  TTree *tree=((TChain*)task->GetInputData(0))->GetTree();
  if (!tree) {
  AliError("Tree not found in ESDhandler");
    return;
  }
 	 
  TList *userInfo=(TList*)tree->GetUserInfo();
  if (!userInfo) {
    AliError("No UserInfo found in tree");
    return;
  }

  TList *cdbList=(TList*)userInfo->FindObject("cdbList");
  if (!cdbList) {
    AliError("No cdbList found in UserInfo");
    if (AliLog::GetGlobalLogLevel()>=AliLog::kError) userInfo->Print();
    return;
  }
 	
  TIter nextCDB(cdbList);
  TObjString *os=0x0;
  while ( (os=(TObjString*)nextCDB()) ){
    if(os->GetString().Contains("TRD/Calib/ChamberGainFactor")){
      // Get Old gain calibration
      AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
 	   
      AliCDBEntry *entry=fTender->GetCDBManager()->Get(id->GetPath(), id->GetFirstRun(), id->GetVersion());
      if (!entry) {
        AliError("No previous gain calibration entry found");
        return;
      }

      fChamberGainOld = dynamic_cast<AliTRDCalDet *>(entry->GetObject());
 	   
      AliDebug(1, Form("Used old Gain entry: %s\n",entry->GetId().ToString().Data()));
    } else if(os->GetString().Contains("TRD/Calib/ChamberVdrift")){
      // Get Old drift velocity calibration
      AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
 	   
      AliCDBEntry *entry=fTender->GetCDBManager()->Get(id->GetPath(), id->GetFirstRun(), id->GetVersion());
      if (!entry) {
        AliError("No previous drift velocity calibration entry found");
        return;
      }

      fChamberVdriftOld = dynamic_cast<AliTRDCalDet *>(entry->GetObject());
 	   
      AliDebug(1, Form("Used old Drift velocity entry: %s\n",entry->GetId().ToString().Data()));
    
    }
  }

  // Get Latest Gain Calib Object
  AliCDBEntry *entryNew=fTender->GetCDBManager()->Get("TRD/Calib/ChamberGainFactor",fTender->GetRun());
  if (entryNew) {
    AliDebug(1, Form("Used new Gain entry: %s\n",entryNew->GetId().ToString().Data()));
    fChamberGainNew = dynamic_cast<AliTRDCalDet *>(entryNew->GetObject());
  } else
    AliError("No new gain calibration entry found");
  
  // Also get the latest Drift Velocity calibration object
  entryNew=fTender->GetCDBManager()->Get("TRD/Calib/ChamberVdrift",fTender->GetRun());
  if (entryNew) {
    AliDebug(1, Form("Used new Drift velocity entry: %s\n",entryNew->GetId().ToString().Data()));
    fChamberVdriftNew = dynamic_cast<AliTRDCalDet *>(entryNew->GetObject());
  } else
    AliError("No new drift velocity calibration entry found");

  if(!fChamberGainNew || !fChamberVdriftNew) AliError("No recent calibration found");
}

//_____________________________________________________
void AliTRDTenderSupply::ApplyGainCorrection(AliESDtrack * track){
  //
  // Apply new gain factors to the track
  //
  if(!fChamberGainNew || !fChamberGainOld){
    AliError("Cannot apply gain correction.");
    return;
  }
 
  if(!(track->GetStatus() & AliESDtrack::kTRDout)) return;
  Double_t p = track->GetOuterParam() ? track->GetOuterParam()->P() : track->P();
  if(p < fPthreshold) return; // Apply low momentum cutoff

  Bool_t applyCorrectionVdrift = kFALSE;
  if(fChamberVdriftOld && fChamberVdriftNew) applyCorrectionVdrift = kTRUE;

  Int_t chamberID[kNPlanes]; 
  for(Int_t il = 0; il < kNPlanes; il++) chamberID[il] = -1;
  if(!GetTRDchamberID(track, chamberID)) return;
  Int_t nTrackletsPID = 0, nTracklets = track->GetTRDntracklets();
  for(Int_t iplane = 0; iplane < kNPlanes; iplane++){
    if(chamberID[iplane] < 0) continue;

    // Mask out bad chambers
    Bool_t isMasked = kFALSE;
    for(UInt_t icam = 0; icam < fNBadChambers; icam++)
      if(fBadChamberID[icam] == chamberID[iplane]){
        isMasked = kTRUE;
        break;
      }
    if(isMasked){
      for(Int_t islice = 0; islice < track->GetNumberOfTRDslices(); islice++){
        track->SetTRDslice(0, iplane, islice);
      }
      continue;
    }

    // Take old and new gain factor and make ratio
    Double_t facOld = fChamberGainOld->GetValue(chamberID[iplane]);
    Double_t facNew = fChamberGainNew->GetValue(chamberID[iplane]); 
    Double_t correction = facNew/facOld;
    if(applyCorrectionVdrift){
      // apply also correction for drift velocity calibration
      Double_t vDriftOld = fChamberVdriftOld->GetValue(chamberID[iplane]);
      Double_t vDriftNew = fChamberVdriftNew->GetValue(chamberID[iplane]);
      correction *= vDriftNew/vDriftOld;
    }
    AliDebug(2, Form("Applying correction factor %f\n", correction));
    Int_t nSlices = 0;
    for(Int_t islice = 0; islice < track->GetNumberOfTRDslices(); islice++){
      Double_t qslice = track->GetTRDslice(iplane, islice);
      if(qslice <= 0.) continue; 
      track->SetTRDslice(qslice / correction, iplane, islice);
      nSlices++;
    }
    if(nSlices) nTrackletsPID++;
  }
  // Use nTrackletsPID to indicate the number of tracklets from good
  // chambers so they are used for the PID
  track->SetTRDntracklets(nTrackletsPID | (nTracklets << 3));
}

//_____________________________________________________
Bool_t AliTRDTenderSupply::GetTRDchamberID(AliESDtrack * const track, Int_t *detectors) {
  //
  // Calculate TRD chamber ID
  //
  Double_t xLayer[kNPlanes] = {300.2, 312.8, 325.4, 338., 350.6, 363.2};
  Double_t etamin[kNStacks] = {0.536, 0.157, -0.145, -0.527,-0.851};
  Double_t etamax[kNStacks] = {0.851, 0.527, 0.145, -0.157,-0.536};
  for(Int_t ily = 0; ily < kNPlanes; ily++) detectors[ily] = -1;

  const AliExternalTrackParam *trueparam = NULL;
  if(track->GetTPCInnerParam()) trueparam = track->GetTPCInnerParam();
  else if(track->GetOuterParam()) trueparam = track->GetOuterParam();
  else if(track->GetInnerParam()) trueparam = track->GetInnerParam();
  if(!trueparam){
    AliDebug(2, "No Track Params");
    return kFALSE;
  }

  AliExternalTrackParam workparam(*trueparam); // Do calculation on working Copy
  Double_t pos[3];
  Int_t nDet = 0;
  for(Int_t ily = 0; ily < kNPlanes; ily++){
    if(!workparam.PropagateTo(xLayer[ily], fESD->GetMagneticField())) {
      AliDebug(2, "Propagation failed");
      break;
    }
    workparam.GetXYZ(pos);
    Double_t trackAlpha = TMath::ATan2(pos[1], pos[0]);
    if(trackAlpha < 0) trackAlpha = 2 * TMath::Pi() + trackAlpha;
    Double_t secAlpha = 2 * TMath::Pi() / 18.;
   
    Int_t sector = static_cast<Int_t>(trackAlpha/secAlpha);
    Double_t etaTrack = track->Eta();
    Int_t stack = -1;
    for(Int_t istack = 0; istack < 5; istack++){
      if(etaTrack >= etamin[istack] && etaTrack <= etamax[istack]){
        stack = istack;
        break;
      }
    }
    if(stack < 0) {
      AliDebug(2, "Dead Area");
      continue;
    }

    detectors[ily] = sector * kNStacks * kNPlanes + stack * kNPlanes + ily;
    nDet++;
  }
  return nDet ? kTRUE : kFALSE;
}

