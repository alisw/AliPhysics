#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
//
#include "AliAlgSteer.h"
#include "AliAlgDet.h"
#include "AliAlgDetITS.h"
#include "AliAlgDetTPC.h"
#include "AliAlgDetTRD.h"
#include "AliAlgDetTOF.h"
#endif

void alignConf(AliAlgSteer* algSteer);
void ConfigITS(AliAlgSteer* algSteer);
void ConfigTPC(AliAlgSteer* algSteer);
void ConfigTRD(AliAlgSteer* algSteer);
void ConfigTOF(AliAlgSteer* algSteer);
//

void alignConf(AliAlgSteer* algSteer)
{
  //
  algSteer->SetRefOCDBConfigMacro("configRefOCDB.C");
  //algSteer->SetRecoOCDBConfigMacro("configRecoOCDB.C");  
  algSteer->SetRecoOCDBConfigMacro(""); // Use ESD info
  //
  algSteer->AddDetector(AliAlgSteer::kITS);
  algSteer->AddDetector(AliAlgSteer::kTPC);
  algSteer->AddDetector(AliAlgSteer::kTRD);
  algSteer->AddDetector(AliAlgSteer::kTOF);
  algSteer->InitDetectors();
  //
  algSteer->GetDetectorByDetID(AliAlgSteer::kTPC)->SetDisabled();
  //  algSteer->GetDetectorByDetID(AliAlgSteer::kTOF)->SetDisabled();
  //  algSteer->GetDetectorByDetID(AliAlgSteer::kTRD)->SetDisabled();
  //
  ConfigITS(algSteer);
  ConfigTPC(algSteer);
  ConfigTRD(algSteer);
  ConfigTOF(algSteer);
  //
  algSteer->SetVtxMinCont(5);   // accept events with min number of vertexTracks contributors
  algSteer->SetVtxMinContVC(10); // use for vertex constraint only those with min number of contributors
  algSteer->SetMaxDCAforVC(0.1,0.6); // dcaR/Z primary selection to allow vertex constraint
  algSteer->SetMaxChi2forVC(10);     // track-vertex chi2 primary selection to allow vertex constraint

  algSteer->SetCosmicSelStrict(kTRUE); // apply track selection to each leg separately
  //
  algSteer->SetMinDetAccColl(2);       // min number of detectors in track
  algSteer->SetMinDetAccCosm(2);
  //
  algSteer->SetMinPointsColl(6,6);     // min number of points per track Boff/Bon
  algSteer->SetMinPointsCosm(4,4);
  //
  algSteer->SetMPOutType(AliAlgSteer::kMille | AliAlgSteer::kMPRec | AliAlgSteer::kContR);
  //algSteer->SetMPOutType(AliAlgSteer::kMille | AliAlgSteer::kContR);
  //  
  //  algSteer->SetMilleTXT(1);
  //  
}

//======================================================================
void ConfigITS(AliAlgSteer* algSteer)
{
  //
  AliAlgDetITS* det = (AliAlgDetITS*)algSteer->GetDetectorByDetID(AliAlgSteer::kITS);
  if (!det||det->IsDisabled()) return;
  det->SetUseErrorParam(kTRUE);
  //
  det->SetObligatoryColl(kTRUE);
  det->SetObligatoryCosm(kTRUE);
  //
  det->SetTrackFlagSelCosm(AliESDtrack::kITSin);
  det->SetTrackFlagSelColl(AliESDtrack::kITSrefit | AliESDtrack::kTPCrefit);
  //
  det->SetNPointsSelCosm(2);
  det->SetNPointsSelColl(3);
  //
  det->SetITSSelPatternColl(AliAlgDetITS::kSPDAny);
  det->SetITSSelPatternCosm(AliAlgDetITS::kSPDNoSel);
  //
  det->SetAddErrorLr(0,30e-4,200e-4);
  det->SetAddErrorLr(1,30e-4,200e-4);
  det->SetAddErrorLr(2,2000e-4,80e-4);
  det->SetAddErrorLr(3,2000e-4,80e-4);
  det->SetAddErrorLr(4,50e-4,500e-4);
  det->SetAddErrorLr(5,50e-4,500e-4);   
}

//======================================================================
void ConfigTPC(AliAlgSteer* algSteer)
{
  //
  AliAlgDetTPC* det = (AliAlgDetTPC*)algSteer->GetDetectorByDetID(AliAlgSteer::kTPC);
  if (!det||det->IsDisabled()) return;
  //
  det->SetObligatoryColl(kFALSE);
  det->SetObligatoryCosm(kFALSE);
  //
  det->SetTrackFlagSelColl(AliESDtrack::kTPCrefit | AliESDtrack::kITSrefit);
  det->SetTrackFlagSelCosm(AliESDtrack::kTPCin);
  //
  det->SetNPointsSelColl(70);
  det->SetNPointsSelCosm(50);
  //
  det->SetAddError(3,10.); // HUGE errors
}

//======================================================================
void ConfigTRD(AliAlgSteer* algSteer)
{
  //
  AliAlgDetTRD* det = (AliAlgDetTRD*)algSteer->GetDetectorByDetID(AliAlgSteer::kTRD);
  if (!det||det->IsDisabled()) return;
  //
  det->SetObligatoryColl(kFALSE);
  det->SetObligatoryCosm(kFALSE);
  //
  det->SetTrackFlagSelColl(AliESDtrack::kTRDout);
  det->SetTrackFlagSelCosm(AliESDtrack::kTRDout);
  //
  det->SetNPointsSelColl(2);
  det->SetNPointsSelCosm(2);
  //
}

//======================================================================
void ConfigTOF(AliAlgSteer* algSteer)
{
  //
  AliAlgDetTOF* det = (AliAlgDetTOF*)algSteer->GetDetectorByDetID(AliAlgSteer::kTOF);
  if (!det||det->IsDisabled()) return;
  //
  det->SetObligatoryColl(kTRUE);
  det->SetObligatoryCosm(kTRUE);
  //
  det->SetTrackFlagSelColl(AliESDtrack::kTOFout);
  det->SetTrackFlagSelCosm(AliESDtrack::kTOFout);
  //
  det->SetNPointsSelColl(1);
  det->SetNPointsSelCosm(1);
  //
}
