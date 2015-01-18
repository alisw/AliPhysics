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
// TPC tender: reapply pid on the fly                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TGraphErrors.h>

#include <AliDCSSensor.h>
#include <AliGRPObject.h>
#include <AliESDpid.h>
#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrack.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliSplineFit.h>
#include <AliCDBId.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliCDBRunRange.h>
#include <AliTender.h>
#include <AliTPCcalibDButil.h>
#include <AliPID.h>

#include "AliTPCTenderSupply.h"

ClassImp(AliTPCTenderSupply)

AliTPCTenderSupply::AliTPCTenderSupply() :
AliTenderSupply(),
fESDpid(0x0),
fGainNew(0x0),
fGainOld(0x0),
fGainAttachment(0x0),
fIsMC(kFALSE),
fGainCorrection(kTRUE),
fAttachmentCorrection(kFALSE),
fPcorrection(kFALSE),
fMultiCorrection(kFALSE),
fArrPidResponseMaster(0x0),
fMultiCorrMean(0x0),
fMultiCorrSigma(0x0),
fSpecificStorages(0x0),
fDebugLevel(0),
fMip(50),
fGRP(0x0),
fBeamType("PP"),
fLHCperiod(),
fMCperiod(),
fRecoPass(0)
{
  //
  // default ctor
  //
}

//_____________________________________________________
AliTPCTenderSupply::AliTPCTenderSupply(const char *name, const AliTender *tender) :
AliTenderSupply(name,tender),
fESDpid(0x0),
fGainNew(0x0),
fGainOld(0x0),
fGainAttachment(0x0),
fIsMC(kFALSE),
fGainCorrection(kTRUE),
fAttachmentCorrection(kFALSE),
fPcorrection(kFALSE),
fMultiCorrection(kFALSE),
fArrPidResponseMaster(0x0),
fMultiCorrMean(0x0),
fMultiCorrSigma(0x0),
fSpecificStorages(0x0),
fDebugLevel(0),
fMip(50),
fGRP(0x0),
fBeamType("PP"),
fLHCperiod(),
fMCperiod(),
fRecoPass(0)
{
  //
  // named ctor
  //
}

//_____________________________________________________
void AliTPCTenderSupply::Init()
{
  //
  // Initialise TPC tender
  //
  
  AliLog::SetClassDebugLevel("AliTPCTenderSupply",10);
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  //
  // Setup PID object
  //
  fIsMC=mgr->GetMCtruthEventHandler();
  
  // Check if another detector already created the esd pid object
  // if not we create it and set it to the ESD input handler
  fESDpid=fTender->GetESDhandler()->GetESDpid();
  if (!fESDpid) {
    fESDpid=new AliESDpid(fIsMC);
    fTender->GetESDhandler()->SetESDpid(fESDpid);
  } 
  
  //setup specific storages
  if (fSpecificStorages){
    TNamed *storage;
    TIter nextStorage(fSpecificStorages);;
    while ( (storage=(TNamed*)nextStorage()) ){
      fTender->GetCDBManager()->SetSpecificStorage(storage->GetName(),storage->GetTitle());
      AliInfo(Form("Setting specific storage: %s (%s)",storage->GetName(), storage->GetTitle()));
    }
  }

  if (fIsMC){
    //force no gain and P correction in MC
    fGainCorrection=kFALSE;
    fPcorrection=kFALSE;
    fAttachmentCorrection=kFALSE;
  }
}

//_____________________________________________________
void AliTPCTenderSupply::ProcessEvent()
{
  //
  // Reapply pid information
  //
  
  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
  
  //load gain correction if run has changed
  if (fTender->RunChanged()){
    SetBeamType();
    SetRecoInfo();
    if ( fBeamType == "PBPB" ) fMultiCorrection=kTRUE;

    if (fDebugLevel>0) AliInfo(Form("Run Changed (%d)",fTender->GetRun()));
    SetParametrisation();
    if (fGainCorrection) SetSplines();
  }
  
  //
  // get gain correction factor
  //
  Double_t corrFactor = GetGainCorrection();
  Double_t corrAttachSlope = 0;
  Double_t corrGainMultiplicityPbPb=1;
  if (fAttachmentCorrection && fGainAttachment) corrAttachSlope = fGainAttachment->Eval(event->GetTimeStamp());
  if (fMultiCorrection&&fMultiCorrMean) corrGainMultiplicityPbPb = fMultiCorrMean->Eval(GetTPCMultiplicityBin());
  
  //
  // - correct TPC signals
  // - recalculate PID probabilities for TPC
  // - correct TPC signal multiplicity dependence

  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliESDtrack *track=event->GetTrack(itrack);
    const AliExternalTrackParam *inner=track->GetInnerParam();
    
    // skip tracks without TPC information
    if (!inner) continue;

    //calculate total gain correction factor given by
    // o gain calibration factor
    // o attachment correction
    // o multiplicity correction in PbPb
    Float_t meanDrift= 250. - 0.5*TMath::Abs(2*inner->GetZ() + (247-83)*inner->GetTgl());
    Double_t corrGainTotal=corrFactor*(1 + corrAttachSlope*180.)/(1 + corrAttachSlope*meanDrift)/corrGainMultiplicityPbPb;

    // apply gain correction
    track->SetTPCsignal(track->GetTPCsignal()*corrGainTotal ,track->GetTPCsignalSigma(), track->GetTPCsignalN());

    // recalculate pid probabilities
    fESDpid->MakeTPCPID(track);
  }
}

//_____________________________________________________
Double_t AliTPCTenderSupply::GetTPCMultiplicityBin()
{
  //
  // Get TPC multiplicity in bins of 150
  //

  AliESDEvent *event=fTender->GetEvent();
  const AliESDVertex* vertexTPC = event->GetPrimaryVertexTPC();
  Double_t tpcMulti=0.;
  if(vertexTPC){
    Double_t vertexContribTPC=vertexTPC->GetNContributors();
    tpcMulti=vertexContribTPC/150.;
    if (tpcMulti>20.) tpcMulti=20.;
  }
  return tpcMulti;
}

//_____________________________________________________
Double_t AliTPCTenderSupply::GetMultiplicityCorrectionMean(Double_t tpcMulti)
{
  //
  // calculate correction factor for dEdx mean
  //
  Double_t meancorrection;
  if (fIsMC)
  {
    // MC data
    meancorrection=1.00054 + (0.00189566)*tpcMulti + (2.07777e-05)*tpcMulti*tpcMulti;
  }
  else
  {
    // real data
    meancorrection=0.999509 + (-0.00271488)*tpcMulti + (-2.98873e-06)*tpcMulti*tpcMulti;
  }
  
  return meancorrection;
  
}

//_____________________________________________________
Double_t AliTPCTenderSupply::GetMultiplicityCorrectionSigma(Double_t tpcMulti)
{
  //
  // calculate correction factor for dEdx sigma
  //
  Double_t sigmacorrection;
  if (fIsMC)
  {
    // MC data
    sigmacorrection=0.95972 + 0.0103721*tpcMulti;
  }
  else
  {
  // real data
    sigmacorrection=1.01817 + 0.0143673*tpcMulti;
  }
  return sigmacorrection;
  
}

//_____________________________________________________
void AliTPCTenderSupply::SetSplines()
{
  //
  // Get Gain splines from OCDB
  //
  
  AliInfo("Update Gain splines");
  
  //
  // Get GPR info for pressure correction
  //
  fPcorrection=kFALSE;
  
  AliCDBEntry *entryGRP=fTender->GetCDBManager()->Get("GRP/GRP/Data",fTender->GetRun());
  if (!entryGRP) {
    AliError("No new GRP entry found");
  } else {
    fGRP = (AliGRPObject*)entryGRP->GetObject();
  }
  if (fDebugLevel>1) AliInfo(Form("GRP entry used: %s",entryGRP->GetId().ToString().Data()));
  
  fGainNew=0x0;
  fGainOld=0x0;
  //
  //find previous entry from the UserInfo
  //
//   TTree *tree=((TChain*)fTender->GetInputData(0))->GetTree();
  AliAnalysisManager*mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisTaskSE *task = (AliAnalysisTaskSE*)mgr->GetTasks()->First();
  if (!task) return;
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
    if (!(os->GetString().Contains("TPC/Calib/TimeGain"))) continue;
    AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
    
    AliCDBEntry *entry=fTender->GetCDBManager()->Get(*id);
    if (!entry) {
      AliError("No previous gain calibration entry found");
      return;
    }
    
    if (fDebugLevel>1) AliInfo(Form("Used old Gain entry: %s",entry->GetId().ToString().Data()));
    
    TObjArray *arr=(TObjArray *)entry->GetObject();
    if (!arr) {
      AliError("Gain Splines array not found in calibration entry");
      return;
    }
    
    AliSplineFit *fit=(AliSplineFit*)arr->At(0);
    if (!fit) {
      AliError("Spline fit not found in array");
      return;
    }
    
    fGainOld = fit;
    delete id;
    break;
  }
  
  //
  //new gain correction
  //

  // This is in principle sill needed for the 2009 data and the 2010 b+c pass1 data
  //   however not supported any longer.
  // For the LHC10c pass2 there is an exception and we still load the splines, but a defined version
  // In case the attachment correction should be check again, this part of the code needs to be changed
  //   in order to load the gain entry again
  Bool_t special10cPass2=fLHCperiod=="LHC10C" && fRecoPass==2;
              
  AliCDBEntry *entryNew=0x0;
  if (special10cPass2) {
    entryNew=fTender->GetCDBManager()->Get("TPC/Calib/TimeGain",fTender->GetRun(),8);
  }
  if (!entryNew) {
    AliError("No new gain calibration entry found");
    return;
  }
  if (fDebugLevel>1) AliInfo(Form("Used new Gain entry: %s",entryNew->GetId().ToString().Data()));
  
  if (entryNew->GetId().GetLastRun()==AliCDBRunRange::Infinity()){
    if (fDebugLevel>0) AliInfo("Use P correction");
    fPcorrection=kTRUE;
  }
  
  TObjArray *arrSplines=(TObjArray *)entryNew->GetObject();
  if (!arrSplines) {
    AliError("Gain Splines array not found in new calibration entry");
    return;
  }
  
  fGainNew        = (AliSplineFit*)arrSplines->At(0);
  fGainAttachment = (TGraphErrors*)arrSplines->FindObject("TGRAPHERRORS_MEAN_ATTACHMENT_BEAM_ALL");
  
  if (!fGainNew) AliError("No recent spline fit object found");
}

//_____________________________________________________
Double_t AliTPCTenderSupply::GetGainCorrection()
{
  //
  // Calculate gain correction factor
  //
  AliESDEvent *event=fTender->GetEvent();
  UInt_t time=event->GetTimeStamp();
  
  Double_t gain=1;
  
  
  if (fGainOld){
    //TODO TODO TODO
    //first correction for the eval const problem
    // needs to be removed when the fix is in AliROOT and
    // the production was done with EvalGraphConst
    // This should be the case from V4-20-Rev11 on
    //TODO TODO TODO
    if ( fLHCperiod.Contains("LHC09") ||
         ((fLHCperiod=="LHC10B" || fLHCperiod=="LHC10C" || fLHCperiod=="LHC10D") && fRecoPass==2) ||
         (fLHCperiod=="LHC10E" && fRecoPass==1) ||
         (fLHCperiod=="LHC10H" && fRecoPass==1)
         ) {
      Double_t valDefault = fGainOld->Eval(time);
      Double_t valConst   = AliTPCcalibDButil::EvalGraphConst(fGainOld, time);
      gain = valDefault/valConst;
    }
    
    if (fGainNew){
      gain *= AliTPCcalibDButil::EvalGraphConst(fGainOld,time)/AliTPCcalibDButil::EvalGraphConst(fGainNew,time);
    }
  }
  
  //If there is only the default calibration, at least apply correction for pressure
  if (fPcorrection){
    if (fGRP) {
      Double_t pressure=fGRP->GetCavernAtmosPressure()->GetValue(time);
//       gain=fGainOld->Eval(time)/(7.03814-0.00459798*pressure)/49.53*fMip;
//       gain=fGainOld->Eval(time)/(7.03814-0.00459798*pressure)/51.51*fMip;
      gain=AliTPCcalibDButil::EvalGraphConst(fGainOld,time)/(7.03814-0.00459798*pressure)/51.51*fMip;
    }
  }
  
  return gain;
}

//_____________________________________________________
void AliTPCTenderSupply::SetBeamType()
{
  //
  // Set the beam type
  //

  fBeamType=fTender->GetEvent()->GetBeamType();
  if (fBeamType.IsNull()||fBeamType.Contains("No Beam")) fBeamType="p-p";
  fBeamType.ToUpper();
  fBeamType.ReplaceAll("-","");
}

//_____________________________________________________
void AliTPCTenderSupply::SetRecoInfo()
{
  //
  // Set reconstruction information
  //

  //reset information
  fRecoPass=0;
  fLHCperiod="";
  fMCperiod="";
  
  //Get the current file to check the reconstruction pass (UGLY, but not stored in ESD... )
  AliESDInputHandler *esdIH = dynamic_cast<AliESDInputHandler*> (fTender->GetESDhandler());
  if (!esdIH) return;
  TTree *tree= (TTree*)esdIH->GetTree();
  TFile *file= (TFile*)tree->GetCurrentFile();
  if (!file) {
    AliError("Current file not found, cannot set reconstruction information");
    return;
  }
  
  TString fileName(file->GetName());
  
  Int_t run=fTender->GetRun();
  
  
  TPRegexp reg(".*(LHC11[a-z]+[0-9]+[a-z]*)/.*");
  //find the period by run number (UGLY, but not stored in ESD... )
  if (run>=114737&&run<=117223)      { fLHCperiod="LHC10B"; fMCperiod="LHC10D1";  }
  else if (run>=118503&&run<=121040) { fLHCperiod="LHC10C"; fMCperiod="LHC10D1";  }
  else if (run>=122195&&run<=126437) { fLHCperiod="LHC10D"; fMCperiod="LHC10F6A"; }
  else if (run>=127719&&run<=130850) { fLHCperiod="LHC10E"; fMCperiod="LHC10F6A"; }
  else if (run>=133004&&run<=135029) { fLHCperiod="LHC10F"; fMCperiod="LHC10F6A"; }
  else if (run>=135654&&run<=136377) { fLHCperiod="LHC10G"; fMCperiod="LHC10F6A"; }
  else if (run>=136851&&run<=139517) { 
    fLHCperiod="LHC10H"; 
    fMCperiod="LHC10H8";  
    if (reg.MatchB(fileName)) fMCperiod="LHC11A10";
  }
  else if ( (run>=144871 && run <=146459 ) || ( run >=146686 && run<= 146860) ) {
   // low energy: 146686 - 146860
    fLHCperiod="LHC11A"; fMCperiod="LHC10F6A";
  }
  else if ( run>=148541  ){
    fLHCperiod="LHC11B"; fMCperiod="LHC10F6A";
  }
  
  //exception new pp MC productions from 2011
  if (fBeamType=="PP" && reg.MatchB(fileName)) fMCperiod="LHC11B2";
  
  //find pass from file name (UGLY, but not stored in ESD... )
  if (fileName.Contains("/pass1")) {
    fRecoPass=1;
  } else if (fileName.Contains("/pass2")) {
    fRecoPass=2;
  } else if (fileName.Contains("/pass3")) {
    fRecoPass=3;
  }
  
}

//_____________________________________________________
void AliTPCTenderSupply::SetParametrisation()
{
  //
  // Change PbPb multiplicity gain correction factor
  //
  
  if (fLHCperiod.IsNull()) {
    AliError("No period set, not changing parametrisation");
    return;
  }
  
  //Get CDB Entry with pid response parametrisations
  AliCDBEntry *pidCDB=fTender->GetCDBManager()->Get("TPC/Calib/PidResponse",fTender->GetRun());
  if (!fArrPidResponseMaster && pidCDB){
    fArrPidResponseMaster=dynamic_cast<TObjArray*>(pidCDB->GetObject());
    AliInfo(Form("Using pid response objects: %s",pidCDB->GetId().ToString().Data()));
  }

  if (!fArrPidResponseMaster){
    AliError("No valid PidResponse master found in OCDB");
    return;
  }
  //data type
  TString datatype="DATA";
  TString period=fLHCperiod;
  //in case of mc fRecoPass is per default 1
  if (fIsMC) {
    datatype="MC";
    fRecoPass=1;
    period=fMCperiod;
  }

  // Set PbPb correction
  fMultiCorrMean=(TF1*)fArrPidResponseMaster->FindObject(Form("TF1_%s_ALL_%s_PASS%d_%s_MEAN",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));

  if (fMultiCorrMean) AliInfo(Form("Setting multiplicity correction function: %s",fMultiCorrMean->GetName()));
  
}

//____________________________________________________________
void AliTPCTenderSupply::AddSpecificStorage(const char* cdbPath, const char* storage)
{
  //
  // Add a specific storage to be set up in Init
  //
  if (!fSpecificStorages) fSpecificStorages=new TObjArray;
  fSpecificStorages->Add(new TNamed(cdbPath, storage));
}

