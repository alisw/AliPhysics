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

#include <AliDCSSensor.h>
#include <AliGRPObject.h>
#include <AliESDpid.h>
#include <AliLog.h>
#include <AliESDEvent.h>
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
fGainCorrection(kTRUE),
fPcorrection(kFALSE),
fArrPidResponseMaster(0x0),
fDebugLevel(0),
fMip(50),
fGRP(0x0)
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
fGainCorrection(kTRUE),
fPcorrection(kFALSE),
fArrPidResponseMaster(0x0),
fDebugLevel(0),
fMip(50),
fGRP(0x0)
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
  
  // Check if another detector already created the esd pid object
  // if not we create it and set it to the ESD input handler
  fESDpid=fTender->GetESDhandler()->GetESDpid();
  if (!fESDpid) {
    fESDpid=new AliESDpid;
    fTender->GetESDhandler()->SetESDpid(fESDpid);
  }
  
  //
  //set bethe bloch parameters depending on whether we have MC or real data
  //
  // for the moment we set the values hardwired. In future they should be stored either in
  // the OCDB or an equivalent calibration data base
  //
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  
    // assume data if there is no mc handler
  if (!mgr->GetMCtruthEventHandler()){
    alephParameters[0] = 0.0283086/0.97;
    //alephParameters[0] = 0.0283086;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;
    //temporary solution
    //fESDpid->GetTPCResponse().SetMip(fMip);
    if (fDebugLevel>0) AliInfo(Form("Use Data parametrisation, Mip set to: %.3f\n",fMip));
    //fESDpid->GetTPCResponse().SetMip(49.2);
  } else {
    //force no gain and P correction in MC
    fGainCorrection=kFALSE;
    fPcorrection=kFALSE;
    if (fDebugLevel>0) AliInfo("Use MC parametrisation\n");
  }
  
  fESDpid->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);
  
  //set detector resolution parametrisation
  fESDpid->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
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
    if (fDebugLevel>0) AliInfo(Form("Run Changed (%d)\n",fTender->GetRun()));
    SetParametrisation();
    if (fGainCorrection) SetSplines();
  }
  
  //
  // get gain correction factor
  //
  Double_t corrFactor = 1;
  if (fGainCorrection) corrFactor=GetGainCorrection();
  
  //
  // - correct TPC signals
  // - recalculate PID probabilities for TPC
  //
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliESDtrack *track=event->GetTrack(itrack);
    if (fGainCorrection)
      track->SetTPCsignal(track->GetTPCsignal()*corrFactor,track->GetTPCsignalSigma(),track->GetTPCsignalN());
    fESDpid->MakeTPCPID(track);
  }
  
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
  if (fDebugLevel>1) AliInfo(Form("GRP entry used: %s\n",entryGRP->GetId().ToString().Data()));
  
  fGainNew=0x0;
  fGainOld=0x0;
  //
  //find previous entry from the UserInfo
  //
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
    if (!(os->GetString().Contains("TPC/Calib/TimeGain"))) continue;
    AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
    
    AliCDBEntry *entry=fTender->GetCDBManager()->Get(*id);
    if (!entry) {
      AliError("No previous gain calibration entry found");
      return;
    }
    
    if (fDebugLevel>1) AliInfo(Form("Used old Gain entry: %s\n",entry->GetId().ToString().Data()));
    
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
  AliCDBEntry *entryNew=fTender->GetCDBManager()->Get("TPC/Calib/TimeGain",fTender->GetRun());
  if (!entryNew) {
    AliError("No new gain calibration entry found");
    return;
  }
  if (fDebugLevel>1) AliInfo(Form("Used new Gain entry: %s\n",entryNew->GetId().ToString().Data()));
  
  if (entryNew->GetId().GetLastRun()==AliCDBRunRange::Infinity()){
    if (fDebugLevel>0) AliInfo("Use P correction\n");
    fPcorrection=kTRUE;
  }
  
  TObjArray *arrSplines=(TObjArray *)entryNew->GetObject();
  if (!arrSplines) {
    AliError("Gain Splines array not found in new calibration entry");
    return;
  }
  
  fGainNew = (AliSplineFit*)arrSplines->At(0);
  
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
    //first correction for the eval
    // needs to be removed when the fix is in AliROOT and
    // the production was done with EvalGraphConst
    //TODO TODO TODO
    if (fTender->GetRun()<=138154||fTender->GetRun()==138197){
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
void AliTPCTenderSupply::SetParametrisation()
{
  //
  // Change BB parametrisation for current run
  //

  //Get CDB Entry with pid response parametrisations
  AliCDBEntry *pidCDB=fTender->GetCDBManager()->Get("TPC/Calib/PidResponse",fTender->GetRun());
  if (pidCDB){
    fArrPidResponseMaster=(TObjArray*)pidCDB->GetObject();
  }

  //Get the current file to check the reconstruction pass (UGLY, but not stored in ESD... )
  AliESDInputHandler *esdIH = dynamic_cast<AliESDInputHandler*> (fTender->GetESDhandler());
  if (!esdIH) return;
  TTree *tree= (TTree*)esdIH->GetTree();
  TFile *file= (TFile*) tree->GetCurrentFile();
  if (!file) {
    AliError("File not found, not changing parametrisation");
    return;
  }

  Int_t run=fTender->GetRun();
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  Double_t alephParameters[5]={0,0,0,0,0};

  //find the period by run number (UGLY, but not stored in ESD... )
  TString period;
  if (run>=114737&&run<=117223) period="LHC10B";
  else if (run>=118503&&run<=121040) period="LHC10C";
  else if (run>=122195&&run<=126437) period="LHC10D";
  else if (run>=127719&&run<=130850) period="LHC10E";
  else if (run>=133004&&run<=135029) period="LHC10F";
  else if (run>=135654&&run<=136377) period="LHC10G";
  else if (run>=136851&&run<=139517) period="LHC10H";
  else if (run>=139699) period="LHC11A";
  
  //find pass from file name (UGLY, but not stored in ESD... )
  TString fileName(file->GetName());
  Int_t pass=0;
  if (fileName.Contains("/pass1")) {
    pass=1;
  } else if (fileName.Contains("/pass2")) {
    pass=2;
  }

  //beam type
  TString beamtype=fTender->GetEvent()->GetBeamType();
  if (beamtype.IsNull()||beamtype.Contains("No Beam")) beamtype="p-p";
  beamtype.ToUpper();
  beamtype.ReplaceAll("-","");
  
  //
  // Set default parametrisations for data and MC
  //
  Bool_t isMC=mgr->GetMCtruthEventHandler();
  if (isMC){
    //MC data
    alephParameters[0] = 2.15898e+00/50.;
    alephParameters[1] = 1.75295e+01;
    alephParameters[2] = 3.40030e-09;
    alephParameters[3] = 1.96178e+00;
    alephParameters[4] = 3.91720e+00;
    
  } else {
    //data parametrisation

    
    //use defaut data parametrisation in case no other will be selected
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;

    fESDpid->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
    
    if (pass==1){
      
    } else if (pass==2){
      //find period
      if (run>=114737&&run<=117223){
        //LHC10b
      } else if (run>=118503&&run<=121040) {
        //LHC10c
      } else if (run>=122195){
        //LHC10d +
        // last run in LHC10d: &&run<126437
        alephParameters[0] = 1.63246/50.;
        alephParameters[1] = 2.20028e+01;
        alephParameters[2] = TMath::Exp(-2.48879e+01);
        alephParameters[3] = 2.39804e+00;
        alephParameters[4] = 5.12090e+00;
      
        //
        fESDpid->GetTPCResponse().SetSigma(2.30176e-02, 5.60422e+02);
      }
    }
    
    if ( beamtype == "PBPB" ){
      AliInfo("BETHE-BLOCH parametrization for PbPb !!!!!!!!!!!!!!!!!!!!!!\n");

      alephParameters[0] = 1.25202/50.;   //was 1.79571/55.;
      alephParameters[1] = 2.74992e+01;   //was 22.0028;
      alephParameters[2] = TMath::Exp(-3.31517e+01);  //was1.55354e-11;
      alephParameters[3] = 2.46246;       //was 2.39804; 
      alephParameters[4] = 6.78938;       //was 5.1209;
    }
  }
  
  fESDpid->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);

  //data type
  TString datatype="DATA";
  //in case of mc pass is per default 1
  if (isMC) {
    datatype="MC";
    pass=1;
  }

  //
  //set the new parametrisation
  //

  if (fArrPidResponseMaster){
    TObject *grAll=0x0;
    //for MC don't use period information
    if (isMC) period="[A-Z0-9]*";
    //pattern for the default entry (valid for all particles)
    TPRegexp reg(Form("TSPLINE3_%s_([A-Z]*)_%s_PASS%d_%s_MEAN",datatype.Data(),period.Data(),pass,beamtype.Data()));
    
    //loop over entries and filter them
    for (Int_t iresp=0; iresp<fArrPidResponseMaster->GetEntriesFast();++iresp){
      TObject *responseFunction=fArrPidResponseMaster->At(iresp);
      TString responseName=responseFunction->GetName();
      
      if (!reg.MatchB(responseName)) continue;
      
      TObjArray *arr=reg.MatchS(responseName);
      TString particleName=arr->At(1)->GetName();
      delete arr;
      if (particleName.IsNull()) continue;
      if (particleName=="ALL") grAll=responseFunction;
      else {
        //find particle id
        for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
          TString particle=AliPID::ParticleName(ispec);
          particle.ToUpper();
          if ( particle == particleName ){
            //test if there is already a function set. If yes, cleanup
            TObject *old=const_cast<TObject*>(fESDpid->GetTPCResponse().GetResponseFunction((AliPID::EParticleType)ispec));
            if (old) delete old;
            fESDpid->GetTPCResponse().SetResponseFunction((AliPID::EParticleType)ispec,responseFunction);
            fESDpid->GetTPCResponse().SetUseDatabase(kTRUE);
            AliInfo(Form("Adding graph: %d - %s\n",ispec,responseFunction->GetName()));
            break;
          }
        }
      }
    }

    //set default response function to all particles which don't have a specific one
    if (grAll){
      for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
        if (!fESDpid->GetTPCResponse().GetResponseFunction((AliPID::EParticleType)ispec)){
          fESDpid->GetTPCResponse().SetResponseFunction((AliPID::EParticleType)ispec,grAll);
          AliInfo(Form("Adding graph: %d - %s\n",ispec,grAll->GetName()));
        }
      }
    }
  }
}

