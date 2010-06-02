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
#include <TObjString.h>
#include <TChain.h>

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
#include <AliTender.h>

#include "AliTPCTenderSupply.h"


AliTPCTenderSupply::AliTPCTenderSupply() :
  AliTenderSupply(),
  fESDpid(0x0),
  fGainNew(0x0),
  fGainOld(0x0),
  fGainCorrection(kTRUE),
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
    fESDpid->GetTPCResponse().SetMip(47.9);
    //fESDpid->GetTPCResponse().SetMip(49.2);
  } else {
    //force no gain correction in MC
    fGainCorrection=kFALSE;
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
    if (fGainCorrection) SetSplines();
  }
  
  //
  // get gain correction factor
  //
  Double_t corrFactor = GetGainCorrection();
  
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
  AliCDBEntry *entryGRP=fTender->GetCDBManager()->Get("GRP/GRP/Data",fTender->GetRun());
  if (!entryGRP) {
    AliError("No new GRP entry found");
  } else {
    fGRP = (AliGRPObject*)entryGRP->GetObject();
  }
  
  fGainNew=0x0;
  fGainOld=0x0;
  //
  //find previous entry from the UserInfo
  //
  TTree *tree=((TChain*)fTender->GetInputData(0))->GetTree();
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
  
  Double_t gain = 1;
  if (fGainNew && fGainOld) gain = fGainOld->Eval(time)/fGainNew->Eval(time);
  
  //If there is no new calibration, at least apply correction for pressure
  if (TMath::Abs(gain-1)<1e-20){
    if (fGRP) {
      Double_t pressure=fGRP->GetCavernAtmosPressure()->GetValue(time);
      gain=fGainOld->Eval(time)/(7.03814-0.00459798*pressure)/49.53*48.2;
    }
  }
  
  return gain;
}

