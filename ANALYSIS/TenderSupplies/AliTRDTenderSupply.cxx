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
#include <TDirectory.h>
#include <TFile.h>
#include <TList.h>
#include <TObjString.h>
#include <TTree.h>
#include <TString.h>
#include <TVectorD.h>

#include <AliCDBEntry.h>
#include <AliCDBId.h>
#include <AliCDBManager.h>
#include <AliOADBContainer.h>
#include <AliTRDCalDet.h>
#include "AliTRDonlineTrackMatching.h"

#include <AliLog.h>
#include <TTree.h>
#include <TChain.h>
#include <AliGeomManager.h>
#include <AliPID.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliESDpid.h>
#include <AliESDtrack.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliTrackerBase.h>
#include <AliTRDPIDReference.h>
#include <AliTRDPIDResponse.h>
#include <AliTRDCalChamberStatus.h>
#include <AliTender.h>

#include "AliTRDTenderSupply.h"

ClassImp(AliTRDTenderSupply)

//_____________________________________________________
AliTRDTenderSupply::AliTRDTenderSupply() :
  AliTenderSupply(),
  fESD(NULL),
  fESDpid(NULL),
  fTrdOnlineTrackMatcher(NULL),
  fChamberGainOld(NULL),
  fChamberGainNew(NULL),
  fChamberVdriftOld(NULL),
  fChamberVdriftNew(NULL),
  fRunByRunCorrection(NULL),
  fPIDmethod(k1DLQpid),
  fNormalizationFactor(1.),
  fPthreshold(0.8),
  fNBadChambers(0),
  fGeoFile(NULL),
  fGainCorrection(kTRUE),
//  fLoadReferences(kFALSE),
//  fLoadReferencesFromCDB(kFALSE),
  fLoadDeadChambers(kFALSE),
  fHasReferences(kFALSE),
  fHasNewCalibration(kTRUE),
  fDebugMode(kFALSE),
  fRedoTrdMatching(kTRUE),
  fNameRunByRunCorrection(),
  fNormalizationFactorArray(NULL)
{
  //
  // default ctor
  //
  memset(fBadChamberID, 0, sizeof(Int_t) * kNChambers);
  memset(fSlicesForPID, 0, sizeof(UInt_t) * 2);
}

//_____________________________________________________
AliTRDTenderSupply::AliTRDTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fESD(NULL),
  fESDpid(NULL),
  fTrdOnlineTrackMatcher(NULL),
  fChamberGainOld(NULL),
  fChamberGainNew(NULL),
  fChamberVdriftOld(NULL),
  fChamberVdriftNew(NULL),
  fRunByRunCorrection(NULL),
  fPIDmethod(k1DLQpid),
  fNormalizationFactor(1.),
  fPthreshold(0.8),
  fNBadChambers(0),
  fGeoFile(NULL),
  fGainCorrection(kTRUE),
//  fLoadReferences(kFALSE),
//  fLoadReferencesFromCDB(kFALSE),
  fLoadDeadChambers(kFALSE),
  fHasReferences(kFALSE),
  fHasNewCalibration(kTRUE),
  fDebugMode(kFALSE),
  fRedoTrdMatching(kTRUE),
  fNameRunByRunCorrection(),
  fNormalizationFactorArray(NULL)
{
  //
  // named ctor
  //
  memset(fSlicesForPID, 0, sizeof(UInt_t) * 2);
  memset(fBadChamberID, 0, sizeof(Int_t) * kNChambers);
}

//_____________________________________________________
AliTRDTenderSupply::~AliTRDTenderSupply()
{
  //
  // dtor
  //
    if(fNormalizationFactorArray) delete fNormalizationFactorArray;
    delete fTrdOnlineTrackMatcher;
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
  // Load References
  //if(fLoadReferences && !fLoadReferencesFromCDB) LoadReferences();
  //fESDpid->GetTRDResponse().SetGainNormalisationFactor(fNormalizationFactor);
  //fESDpid->SetTRDslicesForPID(fSlicesForPID[0], fSlicesForPID[1]);

  if(fNameRunByRunCorrection.Length()) LoadRunByRunCorrection(fNameRunByRunCorrection.Data());
  fTrdOnlineTrackMatcher=new AliTRDonlineTrackMatching();
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
    //if(fLoadReferences && !fHasReferences) LoadReferences();
    if(fLoadDeadChambers) LoadDeadChambersFromCDB();
    // Load Geometry
    if(AliGeomManager::GetGeometry()){
      AliInfo("Geometry already loaded by other tenders");
    } else {
      if(fGeoFile) AliInfo(Form("Load geometry from file %s\n", fGeoFile));
      else AliInfo("Load Geometry from OCDB\n");
      AliGeomManager::LoadGeometry(fGeoFile);
    }
  }


  fESD = fTender->GetEvent();
  if (!fESD) return;
  if(fNormalizationFactorArray) fNormalizationFactor = GetNormalizationFactor(fESD->GetRunNumber());
  Int_t ntracks=fESD->GetNumberOfTracks();



  if (fRedoTrdMatching) {
      if (!fTrdOnlineTrackMatcher->ProcessEvent(fESD)) {
	  AliError("TRD online track matching failed!");
      } 
  }



  //
  // recalculate PID probabilities
  //
  Int_t detectors[kNPlanes];
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    for(Int_t idet = 0; idet < 5; idet++) detectors[idet] = -1;
    AliESDtrack *track=fESD->GetTrack(itrack);
    // Recalculate likelihoods
    if(!(track->GetStatus() & AliESDtrack::kTRDout)) continue;
    AliDebug(2, Form("TRD track found, gain correction: %s, Number of bad chambers: %d\n", fGainCorrection ? "Yes" : "No", fNBadChambers));
    if(GetTRDchamberID(track, detectors)){
      if(fGainCorrection && fHasNewCalibration) ApplyGainCorrection(track, detectors);
      if(fNBadChambers) MaskChambers(track, detectors);
    }
    if(fRunByRunCorrection) ApplyRunByRunCorrection(track);
    if(fNormalizationFactor != 1.){
      //printf("Gain Factor: %f\n", fNormalizationFactor);
      // Renormalize charge
      Double_t qslice = -1;
      for(Int_t ily = 0; ily < 6; ily++){
        for(Int_t is = 0; is < track->GetNumberOfTRDslices(); is++){
          qslice = track->GetTRDslice(ily, is);
          //printf("Doing layer %d slice %d, value %f\n", ily, is, qslice);
          if(qslice >0){
            qslice *= fNormalizationFactor;
            //printf("qslice new: %f\n", qslice);
            track->SetTRDslice(qslice, ily, is);
          }
        }
      }
    }
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
void AliTRDTenderSupply::LoadDeadChambersFromCDB(){
  //
  // Load Dead Chambers from the OCDB
  //
  AliDebug(1, "Loading Dead Chambers from the OCDB");
  AliCDBEntry *en = fTender->GetCDBManager()->Get("TRD/Calib/ChamberStatus",fTender->GetRun());
  if(!en){
   AliError("Dead Chambers not in OCDB");
   return;
  }
  en->GetId().Print();

  AliTRDCalChamberStatus* chamberStatus = 0;
  if(en){
    chamberStatus = (AliTRDCalChamberStatus*)en->GetObject();
    if(!chamberStatus) AliError("List with the dead chambers not found");
    for(Int_t ichamber = 0; ichamber < 540; ichamber++) {
      if(!chamberStatus->IsGood(ichamber)){
        //printf("Chamber not installed %d\n",ichamber);
        AddBadChamber(ichamber);
      }
    }
  }
}

/*
//_____________________________________________________
void AliTRDTenderSupply::LoadReferences(){
  //
  // Load Reference from the OCDB/OADB into the PID Response
  //
  if(fLoadReferencesFromCDB){
    AliDebug(1, "Loading Reference Distributions from the OCDB");
    AliCDBEntry *en = fTender->GetCDBManager()->Get("TRD/Calib/PIDLQ1D");
    if(!en){
      AliError("References for 1D Likelihood Method not in OCDB");
      return;
    }
    en->GetId().Print();
    TObjArray *arr = dynamic_cast<TObjArray *>(en->GetObject());
    if(!arr) AliError("List with the references not found");
  
    // Get new references
    TIter refs(arr);
    TObject *o = NULL;
    AliTRDPIDReference *ref = NULL;
    while((o = refs())){
      if(!TString(o->IsA()->GetName()).CompareTo("AliTRDPIDReference")){
        ref = dynamic_cast<AliTRDPIDReference *>(o);
        break;
      }
    }
    if(ref){
      fESDpid->GetTRDResponse().Load(ref);
      fHasReferences = kTRUE;
      AliDebug(1, "Reference distributions loaded into the PID Response");
    } else {
      AliError("References not found");
    }
  } else {
    // Backward compatibility mode
    AliInfo("Loading Reference Distributions from ROOT file");
    fESDpid->GetTRDResponse().Load("$TRAIN_ROOT/util/tender/LQ1dRef_v3.root");
    fHasReferences = kTRUE;
  }
}
*/

//_____________________________________________________
void AliTRDTenderSupply::SetChamberGain(){
  //
  // Load Chamber Gain factors into the Tender supply
  //
  
  //find previous entry from the UserInfo
  TTree *tree=((TChain*)fTender->GetInputData(0))->GetTree();
  if (!tree) {
    fHasNewCalibration = kFALSE;
    AliError("Tree not found in ESDhandler");
    return;
  }
 	 
  TList *userInfo=(TList*)tree->GetUserInfo();
  if (!userInfo) {
    fHasNewCalibration = kFALSE;
    AliError("No UserInfo found in tree");
    return;
  }

  TList *cdbList=(TList*)userInfo->FindObject("cdbList");
  if (!cdbList) {
    fHasNewCalibration = kFALSE;
    AliError("No cdbList found in UserInfo");
    if (AliLog::GetGlobalLogLevel()>=AliLog::kError) userInfo->Print();
    return;
  }
  fHasNewCalibration = kTRUE;
 	
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

  if(!fChamberGainNew || !fChamberVdriftNew){
    AliError("No recent calibration found");
    fHasNewCalibration = kFALSE;
  }
}

//_____________________________________________________
void AliTRDTenderSupply::LoadRunByRunCorrection(const char *filename){
  //
  // Define run by run gain correction for the charge
  // 

  TDirectory *bkp = gDirectory;
  TFile *in = TFile::Open(filename);
  bkp->cd();
  fRunByRunCorrection = dynamic_cast<AliOADBContainer *>(in->Get("TRDchargeCorrection"));
  delete in;
  if(fRunByRunCorrection )
    AliDebug(2, Form("OADB Container has %d runs\n", fRunByRunCorrection->GetNumberOfEntries()));
  /* Temporarily out due to a bug in AliOADBContainer
  fRunByRunCorrection = new AliOADBContainer("TRDchargeCorrection");
  Int_t status = fRunByRunCorrection->InitFromFile(filename, "TRDchargeCorrection");
  if(!status) AliDebug(1, Form("Run-dependend gain correction factors loaded from OADB file %s", filename));
  else{
    AliDebug(1, "Failed Loading Run-dependend gain correction factors");
    delete fRunByRunCorrection;
    fRunByRunCorrection = NULL;
  }
  */
}

//_____________________________________________________
Bool_t AliTRDTenderSupply::IsBadChamber(Int_t chamberID){
  //
  // Check if the chamber id is in the list of bad chambers
  //
  Bool_t isBad = kFALSE;
  for(UInt_t icam = 0; icam < fNBadChambers; icam++)
    if(fBadChamberID[icam] == chamberID){
	    isBad = kTRUE;
	    //printf("cross checking: %i \n",chamberID);
      break;
    }
  return isBad;
}

//_____________________________________________________
void AliTRDTenderSupply::ApplyGainCorrection(AliESDtrack * track, const Int_t * const chamberID){
  //
  // Apply new gain factors to the track
  //
  if(!fChamberGainNew || !fChamberGainOld){
    AliError("Cannot apply gain correction.");
    return;
  }
 
  if(!(track->GetStatus() & AliESDtrack::kTRDout)) return;
  Bool_t applyCorrectionVdrift = kFALSE;
  if(fChamberVdriftOld && fChamberVdriftNew) applyCorrectionVdrift = kTRUE;

  for(Int_t iplane = 0; iplane < kNPlanes; iplane++){
    if(chamberID[iplane] < 0) continue;
    if(IsBadChamber(chamberID[iplane])) continue; // Don't apply gain correction for chambers which are in the list of bad chambers

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
    for(Int_t islice = 0; islice < track->GetNumberOfTRDslices(); islice++){
      Double_t qslice = track->GetTRDslice(iplane, islice);
      if(qslice <= 0.) continue; 
      track->SetTRDslice(qslice / correction, iplane, islice);
    }
  }
}

//_____________________________________________________
void AliTRDTenderSupply::ApplyRunByRunCorrection(AliESDtrack *const track) {
  // 
  // Equalize charge distribution by applying run-by-run correction (multiplicative)
  //

  TVectorD *corrfactor = dynamic_cast<TVectorD *>(fRunByRunCorrection->GetObject(fTender->GetRun()));
  if(!corrfactor){ 
    // No correction available - simply return
    AliDebug(2, "Couldn't derive gain correction factor from OADB");
    return;
  }
  else AliDebug(2, Form("Gain factor from OADB %f", (*corrfactor)[0]));
  Double_t slice = 0;
  for(Int_t ily = 0; ily < kNPlanes; ily++){
    for(Int_t islice = 0; islice < track->GetNumberOfTRDslices(); islice++){
      slice = track->GetTRDslice(ily, islice);
      if(slice < 0.001) continue;   // do not modify slices which are 0 or negative
      slice *= (*corrfactor)[0];
      track->SetTRDslice(slice, ily, islice);
    }
  }
}

//_____________________________________________________
void AliTRDTenderSupply::SetNormalizationFactor(Double_t norm, Int_t runMin, Int_t runMax) {
  //
  // Set the normalisation factor for a given run range
  //
  if(!fNormalizationFactorArray)
    fNormalizationFactorArray = new TObjArray;
  TVectorD *entry = new TVectorD(3);
  TVectorD &myentry = *entry;
  myentry(0) = runMin;
  myentry(1) = runMax;
  myentry(2) = norm;
  fNormalizationFactorArray->Add(entry);
}

//_____________________________________________________
Double_t AliTRDTenderSupply::GetNormalizationFactor(Int_t runnumber){
  // 
  // Load the normalization factor
  //
  Double_t norm = 1.;
  if(fNormalizationFactorArray){
    TVectorD *entry;
    Int_t runMin, runMax;
    TIter entries(fNormalizationFactorArray);
    while((entry = dynamic_cast<TVectorD *>(entries()))){
      TVectorD &myentry = *entry;
      runMin = TMath::Nint(myentry(0));
      runMax = TMath::Nint(myentry(1));
      if(runnumber >= runMin && runnumber <= runMax) norm = myentry(2);
    }
  }
  AliDebug(1, Form("Gain normalization factor: %f\n", norm));
  return norm;
}

//_____________________________________________________
void AliTRDTenderSupply::MaskChambers(AliESDtrack *const track, const Int_t * const chamberID){
  //
  // Mask out chambers which are in the list of bad chambers
  // Set chamber signal to 0 and reduce the number of tracklets used for PID
  //
  AliDebug(2, "Masking bad chambers for TRD track");
  Int_t nTrackletsPID = 0, nslice = 0, nTracklets = track->GetTRDntracklets();
  Bool_t badChamber = kFALSE;
  //Int_t nbad = 0 ;      // Number of bad chambers which contain also a signal
  //Int_t nsliceBad = 0;  // Number of slices in tracklet in a bad chamber
  for(Int_t iplane = 0; iplane < kNPlanes; iplane++){
    badChamber = kFALSE;
    nslice = 0; //nsliceBad = 0;
    if(IsBadChamber(chamberID[iplane])) badChamber = kTRUE;
    for(Int_t islice = 0; islice < track->GetNumberOfTRDslices(); islice++){
      if(badChamber){
        //if(track->GetTRDslice(iplane, islice)) nsliceBad++;
        track->SetTRDslice(-1, iplane, islice);
      } else if(track->GetTRDslice(iplane, islice) > 0.001) nslice++;
    }
    //if(nsliceBad) nbad++;
    if(nslice > 0) nTrackletsPID++;
  }
  //if(nbad) track->SetTRDncls(track->GetTRDncls() - 20 * nbad);      // subtract mean number of clusters per tracklet for bad tracklets 
  // Use nTrackletsPID to indicate the number of tracklets from good
  // chambers so they are used for the PID
  track->SetTRDntracklets(nTrackletsPID | (nTracklets << 3));
}

//_____________________________________________________
Bool_t AliTRDTenderSupply::GetTRDchamberID(AliESDtrack * const track, Int_t *detectors) {
  //
  // Calculate TRD chamber ID
  //
  Double_t p = track->GetOuterParam() ? track->GetOuterParam()->P() : track->P();
  if(p < fPthreshold) return kFALSE; // Apply low momentum cutoff

  Double_t xLayer[kNPlanes] = {300.2, 312.8, 325.4, 338., 350.6, 363.2};
  Double_t etamin[kNStacks] = {0.536, 0.157, -0.145, -0.527,-0.851};
  Double_t etamax[kNStacks] = {0.851, 0.527, 0.145, -0.157,-0.536};
  //Double_t zboundary[kNPlanes] = {302., 317., 328., 343., 350., 350.};
  for(Int_t ily = 0; ily < kNPlanes; ily++) detectors[ily] = -1;

  const AliExternalTrackParam *trueparam = NULL;
  if(track->GetOuterParam()) trueparam = track->GetOuterParam();
  else if(track->GetTPCInnerParam()) trueparam = track->GetTPCInnerParam();
  else if(track->GetInnerParam()) trueparam = track->GetInnerParam();
  if(!trueparam){
    AliDebug(2, "No Track Params");
    return kFALSE;
  }

  AliExternalTrackParam workparam(*trueparam); // Do calculation on working Copy
  Double_t pos[3];
  Int_t nDet = 0;
  for(Int_t ily = 0; ily < kNPlanes; ily++){
    //if(TMath::Abs(workparam.GetZ()) > zboundary[ily]) break;
    //if(!AliTrackerBase::PropagateTrackToBxByBz(&workparam, xLayer[ily], 0.139, 100)){   // Assuming the pion mass
    if(!workparam.PropagateTo(xLayer[ily], fESD->GetMagneticField())) {
      AliDebug(2, "Propagation failed");
      break;
    }
    workparam.GetXYZ(pos);
    Double_t trackAlpha = TMath::ATan2(pos[1], pos[0]);
    if(trackAlpha < 0) trackAlpha = 2 * TMath::Pi() + trackAlpha;
    Double_t secAlpha = 2 * TMath::Pi() / 18.;
   
    Int_t sector = static_cast<Int_t>(trackAlpha/secAlpha);

    if(fDebugMode){
      // Compare to simple propagation without magnetic field
      AliExternalTrackParam workparam1(*trueparam); // Do calculation on working Copy
      Double_t pos1[3];
      //if(!workparam1.PropagateTo(xLayer[ily], fESD->GetMagneticField())) {
      if(!AliTrackerBase::PropagateTrackToBxByBz(&workparam1, xLayer[ily], 0.139, 100)){   // Assuming the pion mass
        AliDebug(2, "Propagation failed");
        break;
      }
      workparam1.GetXYZ(pos1);
      Double_t trackAlpha1 = TMath::ATan2(pos1[1], pos1[0]);
      if(trackAlpha1 < 0) trackAlpha1 = 2 * TMath::Pi() + trackAlpha1;
   
      Int_t sector1 = static_cast<Int_t>(trackAlpha1/secAlpha);
      AliDebug(2, Form("Alpha: Old %f, New %f, diff %f", trackAlpha, trackAlpha1, trackAlpha-trackAlpha1));
      AliDebug(2, Form("Sector: Old %d, New %d", sector, sector1));
    }

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

