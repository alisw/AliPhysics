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
//  Recalculate VZERO timing and decision using the tender                   //
//  (in case the LHC phase drift is updated in OCDB)                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TList.h>
#include <TObjString.h>
#include <TChain.h>
#include <TF1.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliESDVZERO.h>
#include <AliCDBId.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliTender.h>
#include <AliLHCClockPhase.h>
#include <AliVZEROCalibData.h>
#include <AliVZEROTriggerMask.h>
#include <AliVZEROReconstructor.h>

#include "AliVZEROTenderSupply.h"

ClassImp(AliVZEROTenderSupply)

AliVZEROTenderSupply::AliVZEROTenderSupply() :
  AliTenderSupply(),
  fCalibData(NULL),
  fTimeSlewing(NULL),
  fRecoParam(NULL),
  fLHCClockPhase(0),
  fDebug(kFALSE)
{
  //
  // default ctor
  //
}

//_____________________________________________________
AliVZEROTenderSupply::AliVZEROTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fCalibData(NULL),
  fTimeSlewing(NULL),
  fRecoParam(NULL),
  fLHCClockPhase(0),
  fDebug(kFALSE)
{
  //
  // named ctor
  //
}

//_____________________________________________________
void AliVZEROTenderSupply::Init()
{
  //
  // Initialise VZERO tender
  //
}

//_____________________________________________________
void AliVZEROTenderSupply::ProcessEvent()
{
  //
  // Reapply the LHC-clock phase drift
  //

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
  
  //load gain correction if run has changed
  if (fTender->RunChanged()){
    if (fDebug) printf("AliVZEROTenderSupply::ProcessEvent - Run Changed (%d)\n",fTender->GetRun());
    GetPhaseCorrection();

    AliCDBEntry *entryGeom = fTender->GetCDBManager()->Get("GRP/Geometry/Data",fTender->GetRun());
    if (!entryGeom) {
      AliError("No geometry entry is found");
      return;
    } else {
      if (fDebug) printf("AliVZEROTenderSupply::Used geometry entry: %s\n",entryGeom->GetId().ToString().Data());
    }

    AliCDBEntry *entryCal = fTender->GetCDBManager()->Get("VZERO/Calib/Data",fTender->GetRun());
    if (!entryCal) {
      AliError("No VZERO calibration entry is found");
      fCalibData = NULL;
      return;
    } else {
      fCalibData = (AliVZEROCalibData*)entryCal->GetObject();
      if (fDebug) printf("AliVZEROTenderSupply::Used VZERO calibration entry: %s\n",entryCal->GetId().ToString().Data());
    }

    AliCDBEntry *entrySlew = fTender->GetCDBManager()->Get("VZERO/Calib/TimeSlewing",fTender->GetRun());
    if (!entrySlew) {
      AliError("VZERO time slewing function is not found in OCDB !");
      fTimeSlewing = NULL;
      return;
    } else {
      fTimeSlewing = (TF1*)entrySlew->GetObject();
      if (fDebug) printf("AliVZEROTenderSupply::Used VZERO time slewing entry: %s\n",entrySlew->GetId().ToString().Data());
    }

    AliCDBEntry *entryRecoParam = fTender->GetCDBManager()->Get("VZERO/Calib/RecoParam",fTender->GetRun());
    if (!entryRecoParam) {
      AliError("VZERO reco-param object is not found in OCDB !");
      fRecoParam = NULL;
      return;
    } else {
      fRecoParam = (AliVZERORecoParam*)entryRecoParam->GetObject();
      if (fDebug) printf("AliVZEROTenderSupply::Used VZERO reco-param entry: %s\n",entryRecoParam->GetId().ToString().Data());
    }
  }

  if (!fCalibData || !fTimeSlewing || !fRecoParam) {
    AliWarning("VZERO calibration objects not found!");
    return;
  }
  //
  // correct VZERO time signals and decision
  //
  AliESDVZERO *esdVZERO = event->GetVZEROData();
  if (!esdVZERO) {
    AliError("No VZERO object is found inside ESD!");
    return;
  }
  if (!esdVZERO->TestBit(AliESDVZERO::kDecisionFilled)) {
    AliWarning("VZERO offline trigger decisions were not filled in ESD, the tender supply is disabled");
    return;
  }

  if (fDebug) printf("LHC-clock phase correction: %f\n",fLHCClockPhase);

  if (fDebug) printf("Original VZERO decision %d (%f ns) and %d (%f ns)\n",
		     esdVZERO->GetV0ADecision(),esdVZERO->GetV0ATime(),
		     esdVZERO->GetV0CDecision(),esdVZERO->GetV0CTime());
  Float_t time[64];
  for(Int_t i = 0; i < 64; ++i) {
    time[i] = esdVZERO->GetTime(i);
    if (time[i] > (AliVZEROReconstructor::kInvalidTime + 1e-6))
      time[i] += fLHCClockPhase;
  }
  esdVZERO->SetTime(time);

  {
    AliVZEROTriggerMask triggerMask;
    triggerMask.SetRecoParam(fRecoParam);
    triggerMask.FillMasks(esdVZERO, fCalibData, fTimeSlewing);
  }
  if (fDebug) printf("Modified VZERO decision %d (%f ns) and %d (%f ns)\n",
		     esdVZERO->GetV0ADecision(),esdVZERO->GetV0ATime(),
		     esdVZERO->GetV0CDecision(),esdVZERO->GetV0CTime());

}

//_____________________________________________________
void AliVZEROTenderSupply::GetPhaseCorrection()
{
  //
  // Get Gain splines from OCDB
  //

  AliInfo("Get LHC-clock phase correction");

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

  Float_t oldPhase = 0;

  TIter nextCDB(cdbList);
  TObjString *os=0x0;
  while ( (os=(TObjString*)nextCDB()) ){
    if (!(os->GetString().Contains("GRP/Calib/LHCClockPhase"))) continue;
    AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
    
    AliCDBEntry *entry=fTender->GetCDBManager()->Get(*id);
    if (!entry) {
      AliError("The previous LHC-clock phase entry is not found");
      delete id;
      return;
    }
    
    if (fDebug) printf("AliVZEROTenderSupply::Used old LHC-clock phase entry: %s\n",entry->GetId().ToString().Data());
    
    AliLHCClockPhase *phase = (AliLHCClockPhase*)entry->GetObject();
    if (!phase) {
      AliError("Phase object is not found in the calibration entry");
      delete id;
      return;
    }

    oldPhase = phase->GetMeanPhase();

    delete id;
    break;
  }

  //
  //new LHC-clock phase entry
  //
  Float_t newPhase = 0;
  AliCDBEntry *entryNew=fTender->GetCDBManager()->Get("GRP/Calib/LHCClockPhase",fTender->GetRun());
  if (!entryNew) {
    AliError("No new LHC-clock phase calibration entry is found");
    return;
  }
  if (fDebug) printf("AliVZEROTenderSupply::Used new LHC-clock phase entry: %s\n",entryNew->GetId().ToString().Data());
  
  AliLHCClockPhase *phase2 = (AliLHCClockPhase*)entryNew->GetObject();
  if (!phase2) {
    AliError("Phase object is not found in the calibration entry");
    return;
  }

  newPhase = phase2->GetMeanPhase();

  fLHCClockPhase = newPhase - oldPhase;
}
