#include <memory>

#include <TRandom.h>
#include <TTree.h>
#include <TList.h>
#include <TH1.h>
#include <TFile.h>
#include <TString.h>
#include <TLorentzVector.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
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

void AliAnalysisTaskDG::EventInfo::Fill(const AliESDEvent* esdEvent) {
  const AliESDHeader *esdHeader = esdEvent->GetHeader();
  if (NULL == esdHeader) // this is already dealt with in UserExec
    return;

  fClassMask       = esdHeader->GetTriggerMask();
  fClassMaskNext50 = esdHeader->GetTriggerMaskNext50();
  
  fRunNumber       = esdEvent->GetRunNumber();
  
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
  for (Int_t bc=0; bc<21; ++bc)
    fPFBBA[bc] = fPFBBC[bc] = fPFBGA[bc] = fPFBGC[bc] = 0;
}

void AliAnalysisTaskDG::ADV0::FillAD(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
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
  for (Int_t ch=0; ch<4; ++ch) {
    fBB[0] += (esdAD->GetBBFlag(ch  ) && esdAD->GetBBFlag(ch+ 4));
    fBB[1] += (esdAD->GetBBFlag(ch+8) && esdAD->GetBBFlag(ch+12));
    fBG[0] += (esdAD->GetBGFlag(ch  ) && esdAD->GetBGFlag(ch+ 4));
    fBG[1] += (esdAD->GetBGFlag(ch+8) && esdAD->GetBGFlag(ch+12));
  }

  for (Int_t bc=0; bc<21; ++bc) {
    fPFBBA[bc] = fPFBBC[bc] = fPFBGA[bc] = fPFBGC[bc] = 0;
    for (Int_t ch=0; ch<4; ++ch) {
      fPFBBC[bc] += (esdAD->GetPFBBFlag(ch, bc) && esdAD->GetPFBBFlag(ch+4, bc));
      fPFBGC[bc] += (esdAD->GetPFBGFlag(ch, bc) && esdAD->GetPFBGFlag(ch+4, bc));

      fPFBBA[bc] += (esdAD->GetPFBBFlag(ch+8, bc) && esdAD->GetPFBBFlag(ch+12, bc));
      fPFBGA[bc] += (esdAD->GetPFBGFlag(ch+8, bc) && esdAD->GetPFBGFlag(ch+12, bc));
    }
  }
}
void AliAnalysisTaskDG::ADV0::FillV0(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
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

  for (Int_t bc=0; bc<21; ++bc) {
    fPFBBA[bc] = fPFBBC[bc] = fPFBGA[bc] = fPFBGC[bc] = 0;
    for (Int_t ch=0; ch<32; ++ch) {
      fPFBBC[bc] += esdV0->GetPFBBFlag(ch,    bc);
      fPFBGC[bc] += esdV0->GetPFBGFlag(ch,    bc);
      fPFBBA[bc] += esdV0->GetPFBBFlag(ch+32, bc);
      fPFBGA[bc] += esdV0->GetPFBGFlag(ch+32, bc);
    }
  }
}

void AliAnalysisTaskDG::FMD::Fill(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
  fA = trigAna.FMDTrigger(esdEvent, AliTriggerAnalysis::kASide);
  fC = trigAna.FMDTrigger(esdEvent, AliTriggerAnalysis::kCSide);
}

void AliAnalysisTaskDG::TrackData::Fill(AliESDtrack *tr, AliPIDResponse *pidResponse=NULL) {
  if (NULL == tr || NULL == pidResponse) {
    AliError(Form("tr=%p pidResponse=%p", tr, pidResponse));
    return;
  }
  fSign = tr->GetSign();
  fPx   = tr->Px();
  fPy   = tr->Py();  
  fPz   = tr->Pz();  
  fITSsignal = tr->GetITSsignal();
  fTPCsignal = tr->GetTPCsignal();
  fTOFsignal = tr->GetTOFsignal();

  fPIDStatus[0] = pidResponse->CheckPIDStatus(AliPIDResponse::kITS, tr);
  fPIDStatus[1] = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, tr);
  fPIDStatus[2] = pidResponse->CheckPIDStatus(AliPIDResponse::kTOF, tr);

  for (Int_t i=0; i<AliPID::kSPECIES; ++i) {
    const AliPID::EParticleType particleType(static_cast<const AliPID::EParticleType>(i));
    fNumSigmaITS[i] = pidResponse->NumberOfSigmasITS(tr, particleType);
    fNumSigmaTPC[i] = pidResponse->NumberOfSigmasTPC(tr, particleType);
    fNumSigmaTOF[i] = pidResponse->NumberOfSigmasTOF(tr, particleType);
    AliInfo(Form("%d %f %f %f", i, fNumSigmaITS[i], fNumSigmaTPC[i], fNumSigmaTOF[i]));
  }
}

AliAnalysisTaskDG::AliAnalysisTaskDG(const char *name)
  : AliAnalysisTaskSE(name)
  , fIsMC(kFALSE)
  , fTreeBranchNames("")
  , fTrackCutType("TPCOnly")
  , fTriggerSelection("")
  , fCDBStorage("raw://")
  , fMaxTrackSave(4)
  , fTriggerAnalysis()
  , fAnalysisUtils()
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
  , fTriggerIRs("AliTriggerIR", 3)
  , fTrackData("AliAnalysisTaskDG::TrackData", fMaxTrackSave)
  , fMCTracks("TLorentzVector", 2)
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
  if (fTreeBranchNames.Contains("TOFHeader")) {
    t->Branch("TOFHeader",         &fTOFHeader, 32000, 0);
  }
  if (fTreeBranchNames.Contains("TriggerIR")) {
    t->Branch("TriggerIRs",        &fTriggerIRs, 32000, 0);
  }

  if (fTreeBranchNames.Contains("Tracks")) {
    t->Branch("AliAnalysisTaskDG::TrackData", &fTrackData);
  }

  if (fIsMC) {
    t->Branch("TLorentzVector", &fMCTracks, 32000, 0);
  }

}

void AliAnalysisTaskDG::UserCreateOutputObjects()
{
  if (fTrackCutType == "ITSPureSA") {
    fTrackCuts  = AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(kTRUE, kFALSE);
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kBoth);
  }

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
  if (fUseTriggerMask && fCDBStorage != "NONE") {
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
 
    std::unique_ptr<const TObjArray> split(fTriggerSelection.Tokenize("|"));
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
  }
}

void AliAnalysisTaskDG::UserExec(Option_t *)
{
  AliVEvent *event = InputEvent();
  if (NULL == event) {
    AliFatal("NULL == event");
    return;
  }

  // ESD Event
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (NULL == esdEvent) {
    AliFatal("NULL == esdEvent");
    return;
  }
  
  // input handler
  const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
  if (NULL == man) {
    AliFatal("NULL == man");
    return;
  }

  AliESDInputHandler* inputHandler(dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler()));  
  if (NULL == inputHandler) {
    AliFatal("NULL == inputHandler");
    return;
  }

  AliPIDResponse* pidResponse = inputHandler->GetPIDResponse();
  if (NULL == pidResponse) {
    AliFatal("NULL == pidResponse");
    return;
  }

  const AliESDHeader *esdHeader = esdEvent->GetHeader();
  if (NULL == esdHeader) {
    AliFatal("NULL == esdHeader");
    return;
  }

  if (kFALSE == fIsMC && esdHeader->GetEventType() != AliRawEventHeaderBase::kPhysicsEvent)
    return;

  const AliMultiplicity *mult = esdEvent->GetMultiplicity();
  if (NULL == mult) {
    AliFatal("NULL == mult");
    return;
  }

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
    if (fCDBStorage != "NONE") { // OCDB used
      if ((esdHeader->GetTriggerMask()       & fClassMask)       == 0LL &&
	  (esdHeader->GetTriggerMaskNext50() & fClassMaskNext50) == 0LL) {
	AliInfo(Form("not selected: %s", esdEvent->GetFiredTriggerClasses().Data()));
	return;
      } else {
	AliInfo(Form("selected: %s", esdEvent->GetFiredTriggerClasses().Data()));
      }
    } else { // no OCDB used
      std::unique_ptr<const TObjArray> split(fTriggerSelection.Tokenize("|"));
      Bool_t selected = kFALSE;
      for (Int_t i=0, n=split->GetEntries(); i<n && !selected; ++i)
	selected = esdEvent->GetFiredTriggerClasses().Contains(split->At(i)->GetName());
      
      AliInfo(Form("selected: %d %s", selected, esdEvent->GetFiredTriggerClasses().Data()));
      if (!selected)
	return;
    }
  }

  fTreeData.fEventInfo.Fill(esdEvent);

  fTreeData.fIsIncompleteDAQ          = esdEvent->IsIncompleteDAQ();
  fTreeData.fIsSPDClusterVsTrackletBG = fAnalysisUtils.IsSPDClusterVsTrackletBG(esdEvent);
  fTreeData.fIskMB                    = (inputHandler->IsEventSelected() & AliVEvent::kMB);
  AliInfo(Form("inputHandler->IsEventSelected() = %d", inputHandler->IsEventSelected()));
  fTreeData.fV0Info.FillV0(esdEvent, fTriggerAnalysis);
  fTreeData.fADInfo.FillAD(esdEvent, fTriggerAnalysis);

  fVertexSPD    = *(esdEvent->GetPrimaryVertexSPD());
  fVertexTPC    = *(esdEvent->GetPrimaryVertexTPC());
  fVertexTracks = *(esdEvent->GetPrimaryVertexTracks());

  fTOFHeader    = *(esdEvent->GetTOFHeader());

  // store trigger IR for up to +-1 orbits around the event
  fTriggerIRs.Delete();
  for (Int_t i=0,j=0,n=esdHeader->GetTriggerIREntries(); i<n; ++i) {
    const AliTriggerIR *ir = esdHeader->GetTriggerIR(i);
    if (!ir || TMath::Abs(Int_t(ir->GetOrbit()&0xFFFF) - Int_t(fTreeData.fEventInfo.fOrbitID)) > 1)
      continue;
    new(fTriggerIRs[j++]) AliTriggerIR(*ir);
  }

  fFastOrMap    = mult->GetFastOrFiredChips();
  fFiredChipMap = mult->GetFiredChipMap();

  fIR1InteractionMap = esdHeader->GetIRInt1InteractionMap();
  fIR2InteractionMap = esdHeader->GetIRInt2InteractionMap();

  fTreeData.fEventInfo.fnTrklet = mult->GetNumberOfTracklets();

  std::unique_ptr<const TObjArray> oa(fTrackCuts->GetAcceptedTracks(esdEvent));
  fTreeData.fEventInfo.fnTrk = oa->GetEntries();
  fTreeData.fEventInfo.fCharge = 0;
  for (Int_t i=0, n=oa->GetEntries(); i<n; ++i)
    fTreeData.fEventInfo.fCharge += Int_t(dynamic_cast<AliESDtrack*>(oa->At(i))->GetSign());

  fTrackData.Delete();
  if (oa->GetEntries() <= fMaxTrackSave)  {
    for (Int_t i=0, n=TMath::Min(oa->GetEntries(), fMaxTrackSave); i<n; ++i)
      new(fTrackData[i]) TrackData(dynamic_cast<AliESDtrack*>(oa->At(i)), pidResponse);
  }

  if (fIsMC) {
    AliMCEvent *mcEvent = MCEvent();
    if (NULL == mcEvent)
      AliFatal("NULL ==mcEvent");

    fMCTracks.Delete();
    Int_t counter = 0;
    for(Int_t i=0, n=mcEvent->GetNumberOfTracks(); i<n && counter<2; ++i) {
      AliMCParticle *p = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(i));
      if (NULL == p) continue;
      p->Print();
      Printf("A");
      p->Particle()->Print();
      Printf("B");
      new(fMCTracks[counter]) TLorentzVector;
      TLorentzVector *v = dynamic_cast<TLorentzVector*>(fMCTracks.At(counter));
      Printf("B %p", v);
      p->Particle()->Momentum(*v);
      ++counter;
    }
  } // fIsMC

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
