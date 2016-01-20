#include <TRandom.h>
#include <TTree.h>
#include <TList.h>
#include <TH1.h>
#include <TFile.h>
#include <TString.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliTriggerIR.h"
#include "AliESDVertex.h"
#include "AliTriggerIR.h"

#include "AliAnalysisTaskDG.h"
#include "AliRawEventHeaderBase.h"
#include "AliESDVZERO.h"
#include "AliESDAD.h"
#include "AliESDtrackCuts.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"

ClassImp(AliAnalysisTaskDG);
ClassImp(AliAnalysisTaskDG::TreeData);
ClassImp(AliAnalysisTaskDG::TrackData);

void AliAnalysisTaskDG::EventInfo::Fill(const AliESDHeader* esdHeader) {
   fClassMask       = esdHeader->GetTriggerMask();
   fClassMaskNext50 = esdHeader->GetTriggerMaskNext50();

   fL0Inputs        = esdHeader->GetL0TriggerInputs();
   fL1Inputs        = esdHeader->GetL1TriggerInputs();
   fL2Inputs        = esdHeader->GetL2TriggerInputs();

   fBCID            = esdHeader->GetBunchCrossNumber();
   fOrbitID         = esdHeader->GetOrbitNumber();
   fPeriod          = esdHeader->GetPeriodNumber();
   fTimeStamp       = esdHeader->GetTimeStamp();
}

void AliAnalysisTaskDG::ADV0::FillInvalid() {
  fTime[0] = fTime[1] = -10240.0f;
  fBB[0] = fBG[0] = fBB[1] = fBG[1] = -1;
}

void AliAnalysisTaskDG::ADV0::FillV0(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
  fDecisionOnline[0]  = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kCSide, kFALSE);
  fDecisionOnline[1]  = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kASide, kFALSE);

  fDecisionOffline[0] = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kCSide, kTRUE);
  fDecisionOffline[1] = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kASide, kTRUE);

  const AliESDAD *esdAD = esdEvent->GetADData();
  if (NULL == esdAD) {
    FillInvalid();
    return;
  }
  fTime[0] = esdAD->GetADCTime();
  fTime[1] = esdAD->GetADATime();

  fBB[0] = fBB[1] = fBG[0] = fBG[1] = 0;
  for (Int_t ch=0; ch<16; ++ch) {
    fBB[ch/8] += esdAD->GetBBFlag(ch);
    fBG[ch/8] += esdAD->GetBGFlag(ch);
  }
}
void AliAnalysisTaskDG::ADV0::FillAD(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
  fDecisionOnline[0]  = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kCSide, kFALSE);
  fDecisionOnline[1]  = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kASide, kFALSE);

  fDecisionOffline[0] = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kCSide, kTRUE);
  fDecisionOffline[1] = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kASide, kTRUE);

  const AliESDVZERO *esdV0 = esdEvent->GetVZEROData();
  if (NULL == esdV0) {
    FillInvalid();
    return;
  }

  fTime[0] = esdV0->GetV0CTime();
  fTime[1] = esdV0->GetV0ATime();

  fBB[0] = fBB[1] = fBG[0] = fBG[1] = 0;
  for (Int_t ch=0; ch<64; ++ch) {
    fBB[ch/32] += esdV0->GetBBFlag(ch);
    fBG[ch/32] += esdV0->GetBGFlag(ch);
  }
}


void AliAnalysisTaskDG::TrackData::Fill(AliESDtrack *tr) {
  if (NULL == tr)
    return;
  sign = Int_t(tr->GetSign());
  px   = tr->Px();
  py   = tr->Py();  
  pz   = tr->Pz();  
  itsSignal = tr->GetITSsignal();
  tpcSignal = tr->GetTPCsignal();
}

AliAnalysisTaskDG::AliAnalysisTaskDG(const char *name)
  : AliAnalysisTaskSE(name)
  , fTrackCutType("TPCOnly")
  , fTriggerSelection("")
  , fTriggerAnalysis()
  , fList(NULL)
  , fTE(NULL)
  , fIR1InteractionMap()
  , fIR2InteractionMap()
  , fFastOrMap()
  , fFiredChipMap()
  , fVertexSPD()
  , fVertexTPC()
  , fVertexTracks()
  , fTOFHeader()
  , fTriggerIRs("AliTriggerIR", 10)
  , fTrackData("AliAnalysisTaskDG::TrackData", 4)
  , fTrackCuts(NULL)
  , fUseTriggerMask(kFALSE)
  , fClassMask(0ULL)
  , fClassMaskNext50(0ULL)
{  
  for (Int_t i=0; i<kNHist;++i) {
    fHist[i] = NULL;
  }
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

 AliAnalysisTaskDG::~AliAnalysisTaskDG()
{
  const AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (NULL != man && man->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  if (NULL != fList)
    delete fList;
  fList = NULL;

  if (NULL != fTE)
    delete fTE;
  fTE = NULL;

  if (NULL != fTrackCuts)
    delete fTrackCuts;
  fTrackCuts = NULL;

}

void AliAnalysisTaskDG::SetBranches(TTree* t) {
  t->Branch("AliAnalysisTaskDG::TreeData",            &fTreeData);

  if (fTreeBranchNames.Contains("IRMap")) {
    t->Branch("IR1InteractionMap", &fIR1InteractionMap, 32000, 0);
    t->Branch("IR2InteractionMap", &fIR2InteractionMap, 32000, 0);
  }
  if (fTreeBranchNames.Contains("SPDMaps")) {
    t->Branch("FastOrMap",         &fFastOrMap,    32000, 0);
    t->Branch("FiredChipMap",      &fFiredChipMap, 32000, 0);
  }
  if (fTreeBranchNames.Contains("VertexSPD")) {
    t->Branch("VertexSPD",         &fVertexSPD, 32000, 0);
  }
  if (fTreeBranchNames.Contains("VertexTPC")) {
    t->Branch("VertexTPC",         &fVertexTPC, 32000, 0);
  }
  if (fTreeBranchNames.Contains("VertexTracks")) {
    t->Branch("VertexTracks",      &fVertexTracks, 32000, 0);
  }
//   if (fTreeBranchNames.Contains("TOFHeader")) {
//     t->Branch("TOFHeader",         &fTOFHeader, 32000, 0);
//   }
  if (fTreeBranchNames.Contains("TriggerIR")) {
    t->Branch("TriggerIRs",        &fTriggerIRs, 32000, 0);
  }

  if (fTreeBranchNames.Contains("Tracks")) {
    t->Branch("AliAnalysisTaskDG::TrackData", &fTrackData);
  }
}

void AliAnalysisTaskDG::UserCreateOutputObjects()
{
  if (fTrackCutType == "ITSPureSA")
    fTrackCuts  = AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(kTRUE, kFALSE);

  if (fTrackCutType == "TPCOnly")
    fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  
  if (fTrackCutType == "ITSTPC2011")
    fTrackCuts  = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1);

  if (NULL == fTrackCuts) {
    AliFatal(Form("NULL == fTrackCuts (%s)", fTrackCutType.Data()));
  }

  fList = new TList;
  fList->SetOwner(kTRUE);
  fList->SetName(GetListName());
  fHist[kHistTrig] = new TH1D("HTrig", ";trigger class index", 102, -1.5, 100.5);
  fHist[kHistTrig]->SetStats(0);
  fList->Add(fHist[kHistTrig]);
  PostData(1, fList);

  TDirectory *owd = gDirectory;
  TFile *fSave = OpenFile(1);
  fTE = new TTree;
  fTE->SetName(GetTreeName());
  SetBranches(fTE);
  PostData(2, fTE);
  owd->cd();
}

void AliAnalysisTaskDG::NotifyRun()
{
  AliInfo(Form("run %d %s", fCurrentRunNumber, fCDBStorage.Data()));
  fUseTriggerMask = (fTriggerSelection != "");
  if (fUseTriggerMask) {
    fClassMask = fClassMaskNext50 = 0ULL;

    AliCDBManager *man = AliCDBManager::Instance();
    if (NULL == man) {
      AliFatal("NULL == man");
    }
    man->SetDefaultStorage(fCDBStorage);
    man->SetRun(fCurrentRunNumber);
    AliCDBEntry *entry = man->Get("GRP/CTP/Config");
    const AliTriggerConfiguration *triggerConfig  = dynamic_cast<const AliTriggerConfiguration*>(entry->GetObject());
    if (NULL == triggerConfig) {
      AliFatal("NULL == triggerConfig");
    }
 
    TObjArray *split = fTriggerSelection.Tokenize("|");
    const TObjArray &classes = triggerConfig->GetClasses();
    for (Int_t i=0, n=classes.GetEntries(); i<n; ++i) {
      const AliTriggerClass *c = dynamic_cast<const AliTriggerClass*>(classes.At(i));
      for (Int_t j=0, m=split->GetEntries(); j<m; ++j) {
	if (TString(c->GetName()).Contains(split->At(j)->GetName())) {
	  if (c->GetIndex() < 51)
	    fClassMask       |= (1ULL << ULong64_t(c->GetIndex()- 1));
	  else
	    fClassMaskNext50 |= (1ULL << ULong64_t(c->GetIndex()-51));
	  AliInfo(Form("Selected trigger class %s %lld %lld", c->GetName(), fClassMask, fClassMaskNext50));
	}
      }
    }
    delete split;
  }
}

void AliAnalysisTaskDG::UserExec(Option_t *)
{
  AliVEvent *event = InputEvent();
  if (NULL == event) {
     Error("UserExec", "Could not retrieve event");
     return;
  }

  // ESD Event
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (NULL == esdEvent) {
    AliWarning("ESDEvent is not available");
    return;
  }
  
  // input handler
  const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
  if (NULL == man) {
    AliWarning("AliAnalysisManager is not available");
    return;
  }

  AliESDInputHandler* inputHandler(dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler()));  
  if (NULL == inputHandler) {
    AliWarning("AliESDInputHandler is not available");
    return;
  }

  const AliESDHeader *esdHeader = esdEvent->GetHeader();
  if (NULL == esdHeader)
    return;

  if (esdHeader->GetEventType() != AliRawEventHeaderBase::kPhysicsEvent)
    return;

  const AliMultiplicity *mult = esdEvent->GetMultiplicity();
  if (NULL == mult)
    return;


  fHist[kHistTrig]->Fill(-1); // # analyzed events in underflow bin

  for (Int_t i=0; i<50; ++i) {
    const ULong64_t mask(1ULL<<i);
    if ((esdHeader->GetTriggerMask() & mask) == mask)
      fHist[kHistTrig]->Fill(i);
    if ((esdHeader->GetTriggerMaskNext50() & mask) == mask)
      fHist[kHistTrig]->Fill(50+i);
  }

  PostData(1, fList);

  AliInfo(Form("fUseTriggerMask=%d %lld %lld", fUseTriggerMask, fClassMask, fClassMaskNext50));
  if (fUseTriggerMask) {
    if ((esdHeader->GetTriggerMask()       & fClassMask)       == 0LL &&
	(esdHeader->GetTriggerMaskNext50() & fClassMaskNext50) == 0LL) {
      AliInfo(Form("not selected: %s", esdEvent->GetFiredTriggerClasses().Data()));
      return;
    } else {
      AliInfo(Form("selected: %s", esdEvent->GetFiredTriggerClasses().Data()));
    }
  }

  fTreeData.fEventInfo.Fill(esdHeader);

  fTreeData.fIsIncompleteDAQ = esdEvent->IsIncompleteDAQ();

  fTreeData.fV0Info.FillV0(esdEvent, fTriggerAnalysis);
  fTreeData.fADInfo.FillAD(esdEvent, fTriggerAnalysis);

  fVertexSPD    = *(esdEvent->GetPrimaryVertexSPD());
  fVertexTPC    = *(esdEvent->GetPrimaryVertexTPC());
  fVertexTracks = *(esdEvent->GetPrimaryVertexTracks());

//   fTOFHeader    = *(esdEvent->GetTOFHeader());

  fTriggerIRs.Delete();
  for (Int_t i=0; i<esdHeader->GetTriggerIREntries(); ++i) {
    new(fTriggerIRs[i]) AliTriggerIR(*(esdHeader->GetTriggerIR(i)));
  }

  fFastOrMap    = mult->GetFastOrFiredChips();
  fFiredChipMap = mult->GetFiredChipMap();

  fIR1InteractionMap = esdHeader->GetIRInt1InteractionMap();
  fIR2InteractionMap = esdHeader->GetIRInt2InteractionMap();

  fTreeData.fEventInfo.fnTrklet = mult->GetNumberOfTracklets();

  TObjArray* oa = fTrackCuts->GetAcceptedTracks(esdEvent);
  fTreeData.fEventInfo.fnTrk = oa->GetEntries();
  fTreeData.fEventInfo.fCharge = 0;
  for (Int_t i=0, n=oa->GetEntries(); i<n; ++i)
    fTreeData.fEventInfo.fCharge += Int_t(dynamic_cast<AliESDtrack*>(oa->At(i))->GetSign());

  fTrackData.Delete();
  if (oa->GetEntries() <= 4)  {
    Printf("TRK %d", oa->GetEntries());
    for (Int_t i=0, n=TMath::Min(oa->GetEntries(), 10); i<n; ++i)
      new(fTrackData[i]) TrackData(dynamic_cast<AliESDtrack*>(oa->At(i)));
  }
  delete oa;

  fTE->Fill();
  PostData(2, fTE);
}

void AliAnalysisTaskDG::Terminate(Option_t*)
{
  fList  = dynamic_cast<TList*>(GetOutputData(1));
  if (NULL == fList)
    Error("Terminate","fList is not available");

  fTE  = dynamic_cast<TTree*>(GetOutputData(2));
  if (NULL == fTE)
    Error("Terminate","fTE is not available");

}
