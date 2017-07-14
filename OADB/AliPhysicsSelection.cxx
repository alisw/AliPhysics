/* $Id: AliPhysicsSelection.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//                      Implementation of   Class AliPhysicsSelection
// This class selects collision candidates from data runs, applying selection cuts on triggers 
// and background rejection
//
// Usage:
//
// Create the object:
//   fPhysicsSelection = new AliPhysicsSelection;
//
// For MC data, call
//   fPhysicsSelection->SetAnalyzeMC()
//
// To check if an event is a collision candidate, use:
//   fPhysicsSelection->IsCollisionCandidate(fInputEvent)
//
// After processing save the resulting histograms to a file with (a folder physics_selection 
//   will be created that contains the histograms):
//   fPhysicsSelection->SaveHistograms("physics_selection")
//
// To print statistics after processing use:
//   fPhysicsSelection->Print();
//
// The BX ids corresponding to real bunches crossings p2 are
// automatically selected. You cannot process runs with different
// filling schemes if this option is set. If you want to disable this,
// use: 
//   fPhysicsSelection->SetUseBXNumbers(0);
//
// Usually the class selects the trigger scheme by itself depending on the run number.
// Nevertheless, you can do that manually by calling AddCollisionTriggerClass() and AddBGTriggerClass()
// Example:
// To define the class CINT1B-ABCE-NOPF-ALL as collision trigger (those will be accepted as  
// collision candidates when they pass the selection):
//   AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #769 #3119");
// To select on bunch crossing IDs in addition, use:
//   AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #769 #3119");
// To define the class CINT1A-ABCE-NOPF-ALL as a background trigger (those will only be counted
// for the control histograms):
//   AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL");
// You can also specify more than one trigger class in a string or you can require that some are *not*
// present. The following line would require CSMBA-ABCE-NOPF-ALL, but CSMBB-ABCE-NOPF-ALL is not allowed
// to be present:
//   AddBGTriggerClass("+CSMBA-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL");
//
// The class also supports the triggers used in heavy ion runs
//
//   Origin: Jan Fiete Grosse-Oetringhaus, CERN 
//           Michele Floris, CERN
//-------------------------------------------------------------------------
#include <algorithm>
#include <vector>
#include <iterator>
#include <regex>

#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TIterator.h>
#include <TDirectory.h>
#include <TObjArray.h>
#include <TPRegexp.h>
#include <TParameter.h>
#include <TInterpreter.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,3,0)
#include <v5/TFormula.h>
#else
#include <TFormula.h>
#endif

#include "AliPhysicsSelection.h"

#include "AliTriggerAnalysis.h"
#include "AliLog.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "TPRegexp.h"
#include "TFile.h"
#include "AliOADBContainer.h"
#include "AliOADBPhysicsSelection.h"
#include "AliOADBFillingScheme.h"
#include "AliOADBTriggerAnalysis.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGRPObject.h"
#include "AliCDBEntry.h"
#include "AliVZEROTriggerData.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliITSTriggerConditions.h"
ClassImp(AliPhysicsSelection)

AliPhysicsSelection::AliPhysicsSelection() :
AliAnalysisCuts("AliPhysicsSelection", "AliPhysicsSelection"),
fPassName(""),
fCurrentRun(-1),
fMC(kFALSE),
fPileupCutsEnabled(kFALSE),
fIsPP(kFALSE),
fReadOCDB(kFALSE),
fUseBXNumbers(0),
fUsingCustomClasses(0),
fCollTrigClasses(),
fBGTrigClasses(),
fTriggerAnalysis(),
fHistList(),
fHistStat(0),
fPSOADB(0),
fFillOADB(0),
fTriggerOADB(0),
fTriggerToFormula()
{
  // constructor
  fCollTrigClasses.SetOwner(1);
  fBGTrigClasses.SetOwner(1);
  fTriggerAnalysis.SetOwner(1);
  fHistList.SetOwner(1);
  fTriggerToFormula = new StringToFormula();
  
  AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kWarning);
}

AliPhysicsSelection::AliPhysicsSelection(const char *name) :
 AliAnalysisCuts("AliPhysicsSelection", "AliPhysicsSelection"),
 fPassName(""),
 fCurrentRun(-1),
 fMC(kFALSE),
 fPileupCutsEnabled(kFALSE),
 fIsPP(kFALSE),
 fReadOCDB(kFALSE),
 fUseBXNumbers(0),
 fUsingCustomClasses(0),
 fCollTrigClasses(),
 fBGTrigClasses(),
 fTriggerAnalysis(),
 fHistList(),
 fHistStat(0),
 fPSOADB(0),
 fFillOADB(0),
 fTriggerOADB(0),
 fTriggerToFormula()
 {
   // constructor
   fCollTrigClasses.SetOwner(1);
   fBGTrigClasses.SetOwner(1);
   fTriggerAnalysis.SetOwner(1);
   fHistList.SetOwner(1);
   fTriggerToFormula = new StringToFormula();
   AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kWarning);
 }

AliPhysicsSelection::~AliPhysicsSelection(){
  if (fPSOADB)       delete fPSOADB;
  if (fFillOADB)     delete fFillOADB;
  if (fTriggerOADB)  delete fTriggerOADB;
  delete fTriggerToFormula;
}

UInt_t AliPhysicsSelection::CheckTriggerClass(const AliVEvent* event, const char* trigger, Int_t& triggerLogic) const {
  // checks if the given trigger class(es) are found for the current event
  // format of trigger: +TRIGGER1,TRIGGER1b,TRIGGER1c -TRIGGER2 [#XXX] [&YY] [*ZZ]
  //   requires one out of TRIGGER1,TRIGGER1b,TRIGGER1c and rejects TRIGGER2
  //   in bunch crossing XXX
  //   if successful, YY is returned (for association between entry in fCollTrigClasses and AliVEvent::EOfflineTriggerTypes)
  //   triggerLogic is filled with ZZ, defaults to 0
  
  TString classes = event->GetFiredTriggerClasses();
  
  Bool_t foundBCRequirement = kFALSE;
  Bool_t foundCorrectBC = kFALSE;
  
  UInt_t returnCode = AliVEvent::kUserDefined;
  triggerLogic = 0;
  
  AliDebug(AliLog::kDebug+1, Form("Processing event with triggers %s", event->GetFiredTriggerClasses().Data()));
  
  TString str(trigger);
  TObjArray* tokens = str.Tokenize(" ");
  
  for (Int_t i=0; i < tokens->GetEntries(); i++) {
    TString str2(((TObjString*) tokens->At(i))->String());
    if (str2[0] == '+' || str2[0] == '-') {
      Bool_t flag = (str2[0] == '+');
      str2.Remove(0, 1);
      TObjArray* tokens2 = str2.Tokenize(",");
      Bool_t foundTriggerClass = kFALSE;
      for (Int_t j=0; j < tokens2->GetEntries(); j++) {
        TString str3(((TObjString*) tokens2->At(j))->String());
        if (flag && classes.Contains(str3)) foundTriggerClass = kTRUE;
        if (!flag && classes.Contains(str3)) {
          AliDebug(AliLog::kDebug+1, Form("Rejecting event because trigger class %s is present", str3.Data()));
          delete tokens2;
          delete tokens;
          return kFALSE;
        }
      }
      
      delete tokens2;
      
      if (flag && !foundTriggerClass) {
        AliDebug(AliLog::kDebug+1, Form("Rejecting event because (none of the) trigger class(es) %s is present", str2.Data()));
        delete tokens;
        return kFALSE;
      }
    }
    else if (str2[0] == '#')
    {
      foundBCRequirement = kTRUE;
      
      str2.Remove(0, 1);
      
      Int_t bcNumber = str2.Atoi();
      AliDebug(AliLog::kDebug+1, Form("Checking for bunch crossing number %d", bcNumber));
      
      if (event->GetBunchCrossNumber() == bcNumber)
      {
        foundCorrectBC = kTRUE;
        AliDebug(AliLog::kDebug+1, Form("Found correct bunch crossing %d", bcNumber));
      }
    }
    else if (str2[0] == '&') { str2.Remove(0, 1); returnCode = str2.Atoll();  }
    else if (str2[0] == '*') { str2.Remove(0, 1); triggerLogic = str2.Atoi(); }
    else AliFatal(Form("Invalid trigger syntax: %s", trigger));
  }
  
  delete tokens;
  
  if (foundBCRequirement && !foundCorrectBC) return kFALSE;
  
  return returnCode;
}

/// Evaluate if the given event fulfills a given trigger logic
///
/// \param event Pointer to the current event
/// \param triggerAnalysis Pointer to the TriggerAnlysis class
/// \param triggerLogic Describing trigger logic; e.g. "V0A && V0C && ZDCTime && !TPCHVdip"
/// \param offline Offline analysis(?)
///
/// \return True if the given event matches the trigger logic
Bool_t AliPhysicsSelection::EvaluateTriggerLogic(const AliVEvent* event,
						 AliTriggerAnalysis* triggerAnalysis,
						 const char* triggerLogic, Bool_t offline){
  auto& formula_and_bits = FindForumla(triggerLogic);
  auto& trg_formula = formula_and_bits.first;
  auto& bits = formula_and_bits.second;
  // Get the values for each individual trigger in the trigger logic string;
  // These values are the parameters of the TFormula
  std::vector<Double_t> paras(bits.size());
  auto offline_flag = offline ? AliTriggerAnalysis::kOfflineFlag : 0;
  for (size_t i = 0; i < bits.size(); ++i) {
    typedef AliTriggerAnalysis::Trigger Trigger;
    Trigger bit = static_cast<Trigger>(bits[i] | offline_flag);
    paras[i] = triggerAnalysis->EvaluateTrigger(event, bit);
  }
  Double_t dummy_val[] = {0};
  return trg_formula.EvalPar(dummy_val, paras.data());
}

//______________________________________________________________________________
UInt_t AliPhysicsSelection::IsCollisionCandidate(const AliVEvent* event){
  // checks if the given event is a collision candidate
  // returns a bit word describing the fired offline triggers
  if (fCurrentRun != event->GetRunNumber()) {
    if (!Initialize(event)) AliFatal(Form("Could not initialize for run %d", event->GetRunNumber()));
  }
  
  // check event type; should be PHYSICS = 7 for data and 0 for MC
  Int_t eventType = event->GetHeader()->GetEventType();
  if (fMC) {
    if (eventType != 0) AliFatal(Form("Invalid event type for MC: %d",eventType));
  } else {
    if (eventType != 7) return kFALSE;
  }
  
  UInt_t accept = 0;
  Int_t nColl = fCollTrigClasses.GetEntries();
  Int_t nBG   = fBGTrigClasses.GetEntries();
  for (Int_t i=0; i<nColl+nBG; i++) {
    const char* triggerClass = i<nColl ? fCollTrigClasses.At(i)->GetName() : fBGTrigClasses.At(i-nColl)->GetName();
    AliDebug(AliLog::kDebug+1, Form("Processing trigger class %s", triggerClass));
    
    AliTriggerAnalysis* triggerAnalysis = static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i));
    triggerAnalysis->FillTriggerClasses(event);
    
    Int_t triggerLogic = 0;
    UInt_t singleTriggerResult = CheckTriggerClass(event, triggerClass, triggerLogic);
    if (!singleTriggerResult) continue;
    Bool_t onlineDecision  = EvaluateTriggerLogic(event, triggerAnalysis, fPSOADB->GetHardwareTrigger(triggerLogic), kFALSE);
    Bool_t offlineDecision = EvaluateTriggerLogic(event, triggerAnalysis, fPSOADB->GetOfflineTrigger(triggerLogic), kTRUE);
    triggerAnalysis->FillHistograms(event,onlineDecision,offlineDecision);
    if (!onlineDecision) continue;
    if (!offlineDecision) continue;
    accept |= singleTriggerResult;
  }
  
  if (accept) AliDebug(AliLog::kDebug, Form("Accepted event as collision candidate with bit mask %d", accept));
  return accept;
}


Bool_t AliPhysicsSelection::Initialize(const AliVEvent* event){
  DetectPassName();
  fIsPP = kTRUE;
  if (event->GetDataLayoutType()==AliVEvent::kESD) {
    fIsPP = !TString(((AliESDEvent*) event)->GetESDRun()->GetBeamType()).EqualTo("A-A");
  } else if (fReadOCDB) {
    AliCDBManager* man = AliCDBManager::Instance();
    if (man) {
      if (!man->IsDefaultStorageSet()) man->SetDefaultStorage("raw://");
      man->SetRun(event->GetRunNumber());
      AliGRPObject* grp = (AliGRPObject*) man->Get("GRP/GRP/Data")->GetObject();
      if (grp) fIsPP = !grp->GetBeamType().EqualTo("A-A");
    }
  } else {
    Int_t run = event->GetRunNumber();
    if ((run>=136849 && run<=139517) ||
        (run>=166477 && run<=170593) ||
        (run>=243399 && run<=243984) ||
        (run>=244913 && run<=246994)) fIsPP = kFALSE;
  }
  return Initialize(event->GetRunNumber());
}

Bool_t AliPhysicsSelection::Initialize(Int_t runNumber){
  // initializes the object for the given run  
  AliInfo(Form("Initializing for run %d", runNumber));

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  /// Open OADB file and fetch OADB objects
  TString oadbfilename = AliPhysicsSelection::GetOADBFileName();
  
  TFile * foadb = TFile::Open(oadbfilename);
  if(!foadb->IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));
  
  if(!fPSOADB || !fUsingCustomClasses) { // if it's already set and custom class is required, we use the one provided by the user
    AliInfo("Using Standard OADB");
    AliOADBContainer * psContainer = (AliOADBContainer*) foadb->Get("physSel");
    if (!psContainer) AliFatal("Cannot fetch OADB container for Physics selection");
    fPSOADB = (AliOADBPhysicsSelection*) psContainer->GetObject(runNumber, fIsPP ? "oadbDefaultPP" : "oadbDefaultPbPb",fPassName);
    if (!fPSOADB) AliFatal(Form("Cannot find physics selection object for run %d", runNumber));
  } else {
    AliInfo("Using Custom OADB");
  }
  if(!fFillOADB || !fUsingCustomClasses) { // if it's already set and custom class is required, we use the one provided by the user
    AliOADBContainer * fillContainer = (AliOADBContainer*) foadb->Get("fillScheme");
    if (!fillContainer) AliFatal("Cannot fetch OADB container for filling scheme");
    fFillOADB = (AliOADBFillingScheme*) fillContainer->GetObject(runNumber, "Default",fPassName);
    if (!fFillOADB) AliFatal(Form("Cannot find  filling scheme object for run %d", runNumber));
  }
  if(!fTriggerOADB || !fUsingCustomClasses) { // if it's already set and custom class is required, we use the one provided by the user
    AliOADBContainer * triggerContainer = (AliOADBContainer*) foadb->Get("trigAnalysis");
    if (!triggerContainer) AliFatal("Cannot fetch OADB container for trigger analysis");
    fTriggerOADB = (AliOADBTriggerAnalysis*) triggerContainer->GetObject(runNumber, "Default",fPassName);
    if (!fTriggerOADB) AliFatal(Form("Cannot find  trigger analysis object for run %d", runNumber));
    fTriggerOADB->Print();
  }
  
  if (fMC) {
    // override BX options in case of MC
    fUseBXNumbers = kFALSE;
  }
  
  // initialize first time
  if (fCurrentRun == -1){
    for(UInt_t ibit = 0; ibit < fPSOADB->GetNTriggerBits(); ibit++){
      TObjString * obj = 0;
      TIterator* collIter = fPSOADB->GetCollTrigClass(ibit)->MakeIterator();
      while((obj = (TObjString*) collIter->Next())){
        if (obj->String()!="") fCollTrigClasses.Add(new TObjString(GetTriggerString(obj)));
      }
      if(fMC) continue; // BG classes only make sense for real data
      obj = 0 ;
      TIterator* bgIter = fPSOADB->GetBGTrigClass(ibit)->MakeIterator();
      while((obj = (TObjString*) bgIter->Next())){
        if (obj->String()!= "") fBGTrigClasses.Add(new TObjString(GetTriggerString(obj)));
      }
    }
    
    if (fReadOCDB && 
        (fTriggerOADB->GetV0MOnThreshold()<0 ||
         fTriggerOADB->GetVHMBBAflags()<0  ||
         fTriggerOADB->GetVHMBBCflags()<0  ||
         fTriggerOADB->GetVHMBGAflags()>32 ||
         fTriggerOADB->GetVHMBGCflags()>32)
       ) 
    {
      AliInfo("Setting V0 online thresholds from OCDB");
      AliCDBManager* man = AliCDBManager::Instance();
      if (!man->IsDefaultStorageSet()) man->SetDefaultStorage("raw://");
      man->SetRun(runNumber);
      AliCDBEntry* triggerEntry = man->Get("VZERO/Trigger/Data");
      AliVZEROTriggerData* trigData = triggerEntry ? (AliVZEROTriggerData*) triggerEntry->GetObject() : 0;
      if (trigData) {
        // names in trigData getters are misleading but actual meaning of parameters is ok 
        if (fTriggerOADB->GetV0MOnThreshold()<0) {
          Int_t val = trigData->GetCentralityV0AThrLow();
          AliInfo(Form("  V0MOnThreshold set to %i",val));
          fTriggerOADB->SetV0MOnThreshold(val);
        }
        if (fTriggerOADB->GetVHMBBAflags()<0) {
          Int_t val = trigData->GetMultV0AThrLow();
          AliInfo(Form("  VHMBBAflags set to %i",val));
          fTriggerOADB->SetVHMBBAflags(val);
        }
        if (fTriggerOADB->GetVHMBGAflags()<0) {
          Int_t val = trigData->GetMultV0AThrHigh();
          AliInfo(Form("  VHMBGAflags set to %i",val));
          fTriggerOADB->SetVHMBGAflags(val);
        }
        if (fTriggerOADB->GetVHMBBCflags()>32) {
          Int_t val = trigData->GetMultV0CThrLow();
          AliInfo(Form("  VHMBBCflags set to %i",val));
          fTriggerOADB->SetVHMBBCflags(val);
        }
        if (fTriggerOADB->GetVHMBGCflags()>32) {
          Int_t val = trigData->GetMultV0CThrHigh();
          AliInfo(Form("  VHMBGCflags set to %i",val));
          fTriggerOADB->SetVHMBGCflags(val);
        }
      } else AliError("Failed to set V0 online thresholds from OCDB");
    }
    
    if (fReadOCDB && 
        (fTriggerOADB->GetSH1OuterThreshold()<0 ||
         fTriggerOADB->GetSH2OuterThreshold()<0)
        ){
      AliInfo("Setting FO online thresholds from OCDB");
      AliITSOnlineCalibrationSPDhandler h;
      AliCDBManager* man = AliCDBManager::Instance();
      if (!man->IsDefaultStorageSet()) h.ReadPITConditionsFromDB(runNumber,"raw://");
      else h.ReadPITConditionsFromDB(runNumber,man->GetDefaultStorage()->GetURI().Data());
      AliITSTriggerConditions* tri = h.GetTriggerConditions();
      if (tri) {
        Int_t thresholdInner = tri->GetAlgoParamValueLI("0SH1",0); // algorithm name and param index 0
        if (fTriggerOADB->GetSH1OuterThreshold()<0) {
          Int_t val = tri->GetAlgoParamValueLI("0SH1",1);
          AliInfo(Form("  SH1OuterThreshold set to %i",val));
          fTriggerOADB->SetSH1OuterThreshold(val);
        }
        if (fTriggerOADB->GetSH2OuterThreshold()<0) {
          Int_t val = tri->GetAlgoParamValueLI("0SH2",1);
          AliInfo(Form("  SH2OuterThreshold set to %i",val));
          fTriggerOADB->SetSH2OuterThreshold(val);
        }
      } else  AliError("Failed to set FO online thresholds from OCDB");
    }
    
    for (Int_t i=0; i<fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries(); i++) {
      AliTriggerAnalysis* triggerAnalysis = new AliTriggerAnalysis(Form("TriggerAnalysis%i",i));
      triggerAnalysis->SetParameters(fTriggerOADB);
      triggerAnalysis->SetAnalyzeMC(fMC);
      triggerAnalysis->ApplyPileupCuts(fPileupCutsEnabled);
      triggerAnalysis->EnableHistograms(fIsPP);
      fTriggerAnalysis.Add(triggerAnalysis);
    }
  }
  
  fCurrentRun = runNumber;

  TH1::AddDirectory(oldStatus);
  return kTRUE;
}

void AliPhysicsSelection::FillStatistics(){
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistStat = new TH2F("fHistStat",";;",1,0,1,1,0,1);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,3,0)
  fHistStat->SetCanExtend(TH1::kAllAxes);
#endif
  fHistList.Add(fHistStat);
  TH1::AddDirectory(oldStatus);

  Int_t nColl = fCollTrigClasses.GetEntries();
  for (Int_t i=0; i<nColl; i++) {
    const char* trigger = fCollTrigClasses.At(i)->GetName();
    AliTriggerAnalysis* triggerAnalysis = static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i));
    TH1F* histStat = (TH1F*) triggerAnalysis->GetHistogram("fHistStat");
    if (!histStat) continue;
    Float_t all                 = histStat->GetBinContent(1);
    Float_t accepted            = histStat->GetBinContent(4);
    Float_t v0and = 0;
    Float_t plusNoZDCBG           = 0;
    Float_t plusNoSPDClsVsTrkBG   = 0;
    Float_t plusNoV0C012vsTklBG   = 0;
    Float_t plusNoV0MOnVsOfPileup = 0;
    Float_t plusNoSPDOnVsOfPileup = 0;
    Float_t plusNoSPDVtxPileup    = 0;
    Float_t plusNoV0PFPileup      = 0;
    Float_t plusNoV0Casym         = 0;
    Float_t noSPDClsVsTrkBG   = 0;
    Float_t noV0C012vsTklBG   = 0;
    Float_t noV0MOnVsOfPileup = 0;
    Float_t noSPDOnVsOfPileup = 0;
    Float_t noSPDVtxPileup    = 0;
    Float_t noV0PFPileup      = 0;
    Float_t noV0Casym         = 0;
    Float_t znInTime          = 0;
    Float_t noZNABG           = 0;
    Float_t noZNCBG           = 0;
    Bool_t b[18];
    for (Int_t i=4;i<histStat->GetNbinsX();i++){
      for (Int_t bit=0;bit<18;bit++) b[bit]=i & 1<<bit;
      if (b[ 5]) noSPDClsVsTrkBG  +=histStat->GetBinContent(i+1);
      if (b[ 6]) noV0C012vsTklBG  +=histStat->GetBinContent(i+1);
      if (b[ 7]) noV0MOnVsOfPileup+=histStat->GetBinContent(i+1);
      if (b[ 8]) noSPDOnVsOfPileup+=histStat->GetBinContent(i+1);
      if (b[ 9]) noSPDVtxPileup   +=histStat->GetBinContent(i+1);
      if (b[10]) noV0PFPileup     +=histStat->GetBinContent(i+1);
      if (b[11]) noV0Casym        +=histStat->GetBinContent(i+1);
      if (b[15]) znInTime         +=histStat->GetBinContent(i+1);
      if (b[16]) noZNABG          +=histStat->GetBinContent(i+1);
      if (b[17]) noZNCBG          +=histStat->GetBinContent(i+1);
      if (!b[ 3]) continue;
      if (!b[ 4]) continue;
      v0and+=histStat->GetBinContent(i+1);
      if ( // bad event according to ZDC cuts
          (// Pb-Pb
            !b[15] && ( (fCurrentRun>=136849 && fCurrentRun<=139517) || //10h
                        (fCurrentRun>=166477 && fCurrentRun<=170593) || //11h
                        (fCurrentRun>=243399 && fCurrentRun<=246542) || //15o
                        (fCurrentRun>=246672 && fCurrentRun<=244823) || //15o
                        (fCurrentRun>=244890 && fCurrentRun<=245060) || //15o
                        (fCurrentRun>=245062 && fCurrentRun<=245147) || //15o
                        (fCurrentRun>=245149 && fCurrentRun<=246994)    //15o
                      )
          ) ||
          (// p-Pb
            !b[16] && ( (fCurrentRun>=188144 && fCurrentRun<=188374) || //12g
                        (fCurrentRun>=194713 && fCurrentRun<=196345) || //13bcde
                        (fCurrentRun>=265304 && fCurrentRun<=266318) || //16qr
                        (fCurrentRun>=267132 && fCurrentRun<=267166)    //16t
                      )
          ) ||
          (// Pb-p
            !b[17] && ( (fCurrentRun>=196346 && fCurrentRun<=197411) || //13f
                        (fCurrentRun>=266405 && fCurrentRun<=267131)    //16s
                      )
          )
      ) continue;
      plusNoZDCBG+=histStat->GetBinContent(i+1);
      if (!b[ 5]) continue;
      plusNoSPDClsVsTrkBG+=histStat->GetBinContent(i+1);
      if (!b[ 6]) continue;
      plusNoV0C012vsTklBG+=histStat->GetBinContent(i+1);
      if (!b[ 7]) continue;
      plusNoV0MOnVsOfPileup+=histStat->GetBinContent(i+1);
      if (!b[ 8]) continue;
      plusNoSPDOnVsOfPileup+=histStat->GetBinContent(i+1);
      if (!b[ 9]) continue;
      plusNoSPDVtxPileup+=histStat->GetBinContent(i+1);
      if (!b[10]) continue;
      plusNoV0PFPileup+=histStat->GetBinContent(i+1);
      if (!b[11]) continue;
      plusNoV0Casym+=histStat->GetBinContent(i+1);
    }
    fHistStat->Fill("all",trigger,all);
    fHistStat->Fill("accepted",trigger,accepted);
    fHistStat->Fill("V0A & V0C",trigger,v0and);
    fHistStat->Fill("+ !ZDCBG",trigger,plusNoZDCBG);
    fHistStat->Fill("+ !SPDClsVsTrkBG",trigger,plusNoSPDClsVsTrkBG);
    fHistStat->Fill("+ !V0C012vsTklBG",trigger,plusNoV0C012vsTklBG);
    fHistStat->Fill("+ !V0MOnVsOfPileup",trigger,plusNoV0MOnVsOfPileup);
    fHistStat->Fill("+ !SPDOnVsOfPileup",trigger,plusNoSPDOnVsOfPileup);
    fHistStat->Fill("+ !SPDVtxPileup",trigger,plusNoSPDVtxPileup);
    fHistStat->Fill("+ !V0PFPileup",trigger,plusNoV0PFPileup);
    fHistStat->Fill("+ !V0Casym",trigger,plusNoV0Casym);
    fHistStat->Fill("!SPDClsVsTrkBG",trigger,noSPDClsVsTrkBG);
    fHistStat->Fill("!V0C012vsTklBG",trigger,noV0C012vsTklBG);
    fHistStat->Fill("!V0MOnVsOfPileup",trigger,noV0MOnVsOfPileup);
    fHistStat->Fill("!SPDOnVsOfPileup",trigger,noSPDOnVsOfPileup);
    fHistStat->Fill("!SPDVtxPileup",trigger,noSPDVtxPileup);
    fHistStat->Fill("!V0PFPileup",trigger,noV0PFPileup);
    fHistStat->Fill("!V0Casym",trigger,noV0Casym);
    fHistStat->Fill("ZN time",trigger,znInTime);
  }
  fHistStat->LabelsDeflate("X");
  fHistStat->LabelsDeflate("Y");
}

void AliPhysicsSelection::Print(const Option_t *option) const{
  // print the configuration
  TString msg;
  Printf("Configuration initialized for run %d (MC: %d):", fCurrentRun, fMC);
  msg += Form("Configuration initialized for run %d (MC: %d):\n", fCurrentRun, fMC);
  
  Printf("Collision trigger classes:");
  for (Int_t i=0; i < fCollTrigClasses.GetEntries(); i++)
    Printf("%s", ((TObjString*) fCollTrigClasses.At(i))->String().Data());
  
  Printf("Background trigger classes:");
  for (Int_t i=0; i < fBGTrigClasses.GetEntries(); i++)
    Printf("%s", ((TObjString*) fBGTrigClasses.At(i))->String().Data());
  
  AliTriggerAnalysis* triggerAnalysis = dynamic_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(0));
  
  if (triggerAnalysis) {
    
    Printf("\nTotal available events:");
    
    triggerAnalysis->PrintTriggerClasses();
    // Check if all triggers with counts are known to the physics selection. If not, print a WARNING (only for MC) 
    if(!fMC) {
      TMap * triggers = triggerAnalysis->GetTriggerClasses();
      TIterator* iter = triggers->MakeIterator();
      TObjString* obj = 0;
      static TString alreadyFoundTriggers;
      while ((obj = dynamic_cast<TObjString*> (iter->Next())))
      {
        TString strTrigger = obj->GetString();    
        TParameter<Long64_t>* param = static_cast<TParameter<Long64_t>*> (triggers->GetValue(obj));
        Long_t counts =  (Long_t)param->GetVal();
        TObjArray* tokens = obj->String().Tokenize(" ");
        for (Int_t i=0; i<tokens->GetEntries(); i++)
        {
          TString singleTrigStr = ((TObjString*) tokens->At(i))->String();
          singleTrigStr.Strip(TString::kBoth, ' ' );
          const char * singleTrig = singleTrigStr.Data();
          //	    Printf("%s", singleTrig);
          TString singleTrigStrFull;
          Bool_t found = kFALSE;
          Int_t nCollTriggers = fCollTrigClasses.GetEntries();
          for(Int_t iCollTriggers = 0; iCollTriggers < nCollTriggers; iCollTriggers++){
            singleTrigStrFull = ((TObjString*)fCollTrigClasses.At(iCollTriggers))->String();
            if(singleTrigStrFull.Contains(singleTrigStr)) {
              found = kTRUE;
              break;
            }
            singleTrigStrFull = singleTrigStr;
          }
          Int_t nBGTriggers = fBGTrigClasses.GetEntries();
          for(Int_t iBGTriggers = 0; iBGTriggers < nBGTriggers; iBGTriggers++){
            singleTrigStrFull = ((TObjString*)fBGTrigClasses.At(iBGTriggers))->String();
            if(singleTrigStrFull.Contains(singleTrigStr)) {
              found = kTRUE;
              break;
            }
            singleTrigStrFull = singleTrigStr;
          }
          
          TString blacklist = "CEMC7WU-B-NOPF-ALL, CEMC7WU-AC-NOPF-ALL CEMC7WU-E-NOPF-ALL C0LSR-ABCE-NOPF-TPC CBEAMB-B-NOPF-ALLNOTRD"; // We know we dont support those, so we print no warning	      
          if(counts>0 && !found && !blacklist.Contains(singleTrig) && !singleTrigStr.Contains("WU") && !alreadyFoundTriggers.Contains(singleTrig)) {
            Printf("WARNING: Found unknown trigger [%s] with %ld counts!", singleTrig, counts);
            alreadyFoundTriggers += singleTrig; // Avoid printing warning twice for the same trigger
          }
        }
        delete tokens;	
      }
      delete iter;
    }
  }
  
  if (fUsingCustomClasses) AliWarning("Using custom trigger classes!");
  TString opt(option);
  opt.ToUpper();
  if (opt == "STAT") {
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (mgr) mgr->AddStatisticsMsg(msg);
  }
}

Long64_t AliPhysicsSelection::Merge(TCollection* list) {
  // Merge a list of AliPhysicsSelection objects.
  // Returns the number of merged objects (including this).
  if (!list) return 0;
  if (list->IsEmpty()) return 0;
  
  TIterator* iter = list->MakeIterator();
  TObject* obj;
  TList collTriggerAnalysis;
  TList collHistList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    AliPhysicsSelection* entry = dynamic_cast<AliPhysicsSelection*> (obj);
    if (entry == 0) continue;
    collTriggerAnalysis.Add(&(entry->fTriggerAnalysis));
    collHistList.Add(&(entry->fHistList));
    count++;
  }
  fTriggerAnalysis.Merge(&collTriggerAnalysis);
  fHistList.Merge(&collHistList);
  
  delete iter;
  return count+1;
}

void AliPhysicsSelection::SaveHistograms(const char* folder) {
  // write histograms to current directory
  if (folder) {
    gDirectory->mkdir(folder);
    gDirectory->cd(folder);
  }
  for (Int_t i=0; i < fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries(); i++) {
    TString triggerClass = "trigger_histograms_";
    if (i < fCollTrigClasses.GetEntries())
      triggerClass += ((TObjString*) fCollTrigClasses.At(i))->GetName();
    else
      triggerClass += ((TObjString*) fBGTrigClasses.At(i - fCollTrigClasses.GetEntries()))->String();
    
    gDirectory->mkdir(triggerClass);
    gDirectory->cd(triggerClass);
    
    static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i))->SaveHistograms();
    
    gDirectory->cd("..");
  }
  fHistList.Write();
  if (folder) gDirectory->cd("..");
}


const char * AliPhysicsSelection::GetTriggerString(TObjString * obj) { 
  // Returns a formed objstring
  static TString retString;
  
  retString.Form("%s%s", 
      obj->String().Data(), 
      fUseBXNumbers ? fFillOADB->GetBXIDs(fPSOADB->GetBeamSide(obj->String().Data())) : ""
  );
  
  if (fMC) {
    TPRegexp stripClasses("\\+\\S* ");
    stripClasses.Substitute(retString,"","g");
    stripClasses=TPRegexp("\\-\\S* ");
    stripClasses.Substitute(retString,"","g");
  }
  
  return retString.Data();
}


void AliPhysicsSelection::DetectPassName(){
  if (fMC) return;
  AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!handler) return;
  TObject* prodInfoData = handler->GetUserInfo()->FindObject("alirootVersion");
  TString filePath = handler->GetTree()->GetCurrentFile()->GetName();
  if (prodInfoData) {
    // try to take filePath from UserInfo - available only from ~LHC12d period
    TString str(prodInfoData->GetTitle());
    TObjArray* tokens = str.Tokenize(";");
    for (Int_t i=0;i<=tokens->GetLast();i++) {
      TObjString* stObj = (TObjString*) tokens->At(i);
      TString s = stObj->GetString();
      if (s.Contains("OutputDir")) {
        filePath = s;
        break;
      }
    }
    delete tokens;
  } 
  
  TString passName="";
  
  TObjArray* tokens = filePath.Tokenize("/");
  for (Int_t i=0;i<=tokens->GetLast();i++) {
    TObjString* stObj = (TObjString*) tokens->At(i);
    TString s = stObj->GetString();
    if (s.Contains("pass")) {
      passName = s;
      break;
    }
  }
  delete tokens;
  //
  // temporary patch for LEGO train runners
  //
  if (passName.Contains("_pass")){ // try with "_" as a fallback (as it is the case in the test data of the LEGO train) and do further tokenize
    TObjArray* tokens2 = filePath.Tokenize("_");
    for (Int_t i=0;i<=tokens2->GetLast();i++) {
      TObjString* stObj = (TObjString*) tokens2->At(i);
      TString s = stObj->GetString();
      if (s.Contains("pass")) {
        passName = s;
        break;
      }
    }
    delete tokens2;
  }
  
  
  if (!passName.Contains("pass")){
    AliError(" Failed to find reconstruction pass name:");
    AliError(" --> If these are MC data: please set kTRUE first argument of AddTaskPhysicsSelection");
    AliError(" --> If these are real data: please insert pass name inside the path of your local file, e.g. /your_path/pass2/AliESDs.root");
    AliError(" --> Setting default pass name...");
    passName = "default";
  }
  
  AliInfo(Form("pass name: %s\n",passName.Data()));
  fPassName = passName;
}

FormulaAndBits& AliPhysicsSelection::FindForumla(const char* triggerLogic) {
  // Do we have this logic cached? If not, set it up
  auto it = fTriggerToFormula->find(triggerLogic);
  if (it == fTriggerToFormula->end()) {
    std::string trg_logic_formated;
    std::vector<AliTriggerAnalysis::Trigger> bits;

    TString trigger(triggerLogic);
    TArrayI pos;
    Int_t b = 0;
    Int_t e = 0;
    TPRegexp trigger_regexp("[[:alpha:]][[:alnum:]]*");

    // Go through the matches and construct a new string where each
    // match is replaced by a TFormula parameter
    while (trigger_regexp.Match(trigger, "", e, 1, &pos) != 0) {
      b = e;
      e = pos[0];
      std::string unmatched(trigger.Data() + b, trigger.Data() + e);

      trg_logic_formated.append(unmatched);
      // Each param in the TFormula will be set as the value behind the trigger bit
      trg_logic_formated.append(Form("int([%i])", (Int_t)bits.size()));

      b = e;
      e = pos[1];
      std::string matched(trigger.Data() + b, trigger.Data() + e);

      TInterpreter::EErrorCode error;
      Int_t bit = gInterpreter->ProcessLine(Form("AliTriggerAnalysis::k%s;", matched.c_str()), &error);

      if (error > 0)
	AliFatal(Form("Trigger token %s unknown", matched.c_str()));

      bits.push_back(static_cast<AliTriggerAnalysis::Trigger>(bit));
    }
    trg_logic_formated.append({trigger.Data() + e, trigger.Data() + trigger.Length()});

    R5TFormula formula(Form("dummy_name_%zu", fTriggerToFormula->size()), trg_logic_formated.c_str());

    if (formula.Compile() > 0) {
      AliFatal(Form("Could not evaluate trigger logic %s (evaluated to %s)",
		    triggerLogic, trg_logic_formated.c_str()));
    }
    // Have the iterator point at the newly inserted element so that
    // we don't have to look it up in the return statement
    it = fTriggerToFormula->emplace(std::string(triggerLogic),
				    std::make_pair(formula, std::move(bits))).first;
  }
  return it->second;
}
