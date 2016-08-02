#include <memory>

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TLorentzVector.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliRawEventHeaderBase.h"
#include "AliESDVZERO.h"
#include "AliESDAD.h"

#include "AliAnalysisTaskDiffCrossSections.h"

ClassImp(AliAnalysisTaskDiffCrossSections);
ClassImp(AliAnalysisTaskDiffCrossSections::TreeData);
ClassImp(AliAnalysisTaskDiffCrossSections::MCInfo);

void AliAnalysisTaskDiffCrossSections::MCInfo::Fill(const AliMCEvent* mcEvent, TString mcType) {
  fEventType = kInvalid;
  if (!mcEvent)
    AliFatal("NULL == mcEvent");

  const AliGenEventHeader* h(mcEvent->GenEventHeader());
  if (!h)
    AliFatal("NULL == h");

  Int_t pdgDiff_p = -1;
  if (mcType.Contains("Pythia")) {
    pdgDiff_p = 9902210;
    const AliGenPythiaEventHeader* ph = dynamic_cast<const AliGenPythiaEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case 92:
      fEventType = kSDR;
      break;
    case 93:
      fEventType = kSDL;
      break;
    case 94:
      fEventType = kDD;
      break;
    case 91:
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
   }
  }
  if (mcType.Contains("Pythia8")) {
    pdgDiff_p = 9902210;
    const AliGenPythiaEventHeader* ph = dynamic_cast<const AliGenPythiaEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case 103:
      fEventType = kSDR;
      break;
    case 104:
      fEventType = kSDL;
      break;
    case 105:
      fEventType = kDD;
      break;
    case 106:
      fEventType = kCD;
      break;
    case 102:
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
    }
  }
  if (mcType.Contains("PHOJET")) {
    const AliGenDPMjetEventHeader* ph = dynamic_cast<const AliGenDPMjetEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case  5: 
      fEventType = kSDR;
      break;
    case  6: 
      fEventType = kSDL;
      break;
    case  7: 
      fEventType = kDD;
      break;
    case  4: 
      fEventType = kCD;
      break;
    case  2: 
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
    }    
  }

  const AliStack* stack  = dynamic_cast<const AliStack *>(const_cast<AliMCEvent*>(mcEvent)->Stack());
  if (!stack)
    AliFatal("NULL == stack");
  const Int_t npart = stack->GetNprimary();

  TLorentzVector v;
  fDiffMass[0] = fDiffMass[1] = -1.0f;
  for (Int_t i=0, counter=0, n=stack->GetNprimary(); i<n; ++i) {
    const TParticle *p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p)
      continue;
    if (pdgDiff_p > 0 && p->GetPdgCode() == pdgDiff_p) {
      p->Momentum(v);
      AliInfoF("found diff(%d): m=%f eta=%f", counter, v.M(), v.Eta());
      fDiffMass[v.Eta()>0] = v.M(); // eta>0,<0 -> SDL,SDR to be verified
    }
  }
}

void AliAnalysisTaskDiffCrossSections::EventInfo::Fill(const AliESDEvent* esdEvent) {
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

void AliAnalysisTaskDiffCrossSections::VtxInfo::Fill(const AliESDVertex *vtx) {
  if (vtx) {
    fZ      = vtx->GetZ();
    fNcontr = vtx->GetNContributors();
    AliErrorClass("NULL == vtx");
  } else {
    fZ      =  0.0f;
    fNcontr = -4;
  }
}

void AliAnalysisTaskDiffCrossSections::ADV0::FillInvalid() {
  fTime[0] = fTime[1] = -10240.0f;
  fCharge[0] = fCharge[1] = fBB[0] = fBG[0] = fBB[1] = fBG[1] = -1.0f;
}

void AliAnalysisTaskDiffCrossSections::ADV0::FillAD(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
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
  fCharge[0] = fCharge[1] = 0.0;
  for (Int_t ch=0; ch<16; ++ch)
    fCharge[ch/8] += esdAD->GetMultiplicity(ch);
}
void AliAnalysisTaskDiffCrossSections::ADV0::FillV0(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
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
  fCharge[0] = fCharge[1] = 0.0;
  for (Int_t ch=0; ch<64; ++ch)
    fCharge[ch/32] += esdV0->GetMultiplicity(ch);
}


AliAnalysisTaskDiffCrossSections::AliAnalysisTaskDiffCrossSections(const char *name)
  : AliAnalysisTaskSE(name)
  , fIsMC(kFALSE)
  , fMCType("")
  , fTriggerSelection("")
  , fTriggerAnalysis()
  , fAnalysisUtils()
  , fTE(NULL)
  , fFastOrMap()
  , fFiredChipMap()
  , fTreeData()
  , fMCInfo()
{  
  fTriggerAnalysis.SetAnalyzeMC(fIsMC);

  DefineOutput(1, TTree::Class());
}

AliAnalysisTaskDiffCrossSections::~AliAnalysisTaskDiffCrossSections()
{
  const AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (NULL != man && man->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  if (NULL != fTE)
    delete fTE;
  fTE = NULL;
}

void AliAnalysisTaskDiffCrossSections::SetBranches(TTree* t) {
  t->Branch("AliAnalysisTaskDG::TreeData", &fTreeData);
  t->Branch("FastOrMap",    &fFastOrMap,    32000, 0);
  //  t->Branch("FiredChipMap", &fFiredChipMap, 32000, 0);
  if (fIsMC)
    t->Branch("AliAnalysisTaskDG::MCInfo", &fMCInfo);
}

void AliAnalysisTaskDiffCrossSections::UserCreateOutputObjects()
{
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTE = new TTree;
  fTE->SetName(GetTreeName());
  SetBranches(fTE);
  PostData(1, fTE);
  owd->cd();
}

void AliAnalysisTaskDiffCrossSections::NotifyRun() {
}

void AliAnalysisTaskDiffCrossSections::UserExec(Option_t *)
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

  if (!fIsMC) {
    std::unique_ptr<const TObjArray> split(fTriggerSelection.Tokenize("|"));
    Bool_t selected = kFALSE;
    for (Int_t i=0, n=split->GetEntries(); i<n && !selected; ++i)
      selected = esdEvent->GetFiredTriggerClasses().Contains(split->At(i)->GetName());
    
    AliInfo(Form("selected: %d %s", selected, esdEvent->GetFiredTriggerClasses().Data()));
    if (!selected)
      return;
  }

  fTreeData.fEventInfo.Fill(esdEvent);

  fTreeData.fPhysSelBits              = inputHandler->IsEventSelected();
  fTreeData.fIsIncompleteDAQ          = esdEvent->IsIncompleteDAQ();
  fTreeData.fIsSPDClusterVsTrackletBG = fAnalysisUtils.IsSPDClusterVsTrackletBG(esdEvent);
  fTreeData.fV0Info.FillV0(esdEvent, fTriggerAnalysis);
  fTreeData.fADInfo.FillAD(esdEvent, fTriggerAnalysis);

  fTreeData.fVtxInfo.Fill(esdEvent->GetPrimaryVertexSPD());

  fFastOrMap    = mult->GetFastOrFiredChips();
  fFiredChipMap = mult->GetFiredChipMap();

  fTreeData.fEventInfo.fnTrklet         = mult->GetNumberOfTracklets();
  fTreeData.fEventInfo.fnSPDClusters[0] = mult->GetNumberOfITSClusters(0);
  fTreeData.fEventInfo.fnSPDClusters[1] = mult->GetNumberOfITSClusters(1);

  if (fIsMC)
    fMCInfo.Fill(MCEvent(), fMCType);

  fTE->Fill();
  PostData(1, fTE);
}

void AliAnalysisTaskDiffCrossSections::Terminate(Option_t*)
{
  fTE  = dynamic_cast<TTree*>(GetOutputData(1));
  if (NULL == fTE)
    Error("Terminate","fTE is not available");
}
