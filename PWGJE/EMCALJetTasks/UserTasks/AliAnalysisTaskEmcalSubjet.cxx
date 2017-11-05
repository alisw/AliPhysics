#include <AliVParticle.h>

#include "THistManager.h"

#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"

#include "AliJetContainer.h"
#include "AliEmcalJetFinder.h"

#include "AliAnalysisTaskEmcalSubjet.h"
//=============================================================================

ClassImp(AliAnalysisTaskEmcalSubjet)

//_____________________________________________________________________________
AliAnalysisTaskEmcalSubjet::AliAnalysisTaskEmcalSubjet() :
AliAnalysisTaskEmcalJet(),
fSubjetRadius(0.),
fSubjetAlgorithm(0),
fSubjetFinder(nullptr),
fHistMgr(nullptr)
{
//
//  AliAnalysisTaskEmcalSubjet::AliAnalysisTaskEmcalSubjet() :
//
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalSubjet::AliAnalysisTaskEmcalSubjet(const char *name, const Bool_t bHistos) :
AliAnalysisTaskEmcalJet(name,bHistos),
fSubjetRadius(0.1),
fSubjetAlgorithm(1),
fSubjetFinder(0),
fHistMgr(nullptr)
{
//
//  AliAnalysisTaskEmcalSubjet::AliAnalysisTaskEmcalSubjet(const char *name, const Bool_t bHistos) :
//

  AliAnalysisTaskEmcal::fGeneralHistograms = bHistos;
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalSubjet::~AliAnalysisTaskEmcalSubjet()
{
//
//  AliAnalysisTaskEmcalSubjet::~AliAnalysisTaskEmcalSubjet
//

  if (fHistMgr)      { delete fHistMgr;      fHistMgr      = nullptr; }
  if (fSubjetFinder) { delete fSubjetFinder; fSubjetFinder = nullptr; }
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalSubjet::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskEmcalSubjet::UserCreateOutputObjects
//

  fCreateHisto = kTRUE;
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
//=============================================================================

  if (fSubjetFinder) { delete fSubjetFinder; fSubjetFinder = nullptr; }

  fSubjetFinder = new AliEmcalJetFinder();
  fSubjetFinder->SetRadius(fSubjetRadius);
  fSubjetFinder->SetJetAlgorithm(fSubjetAlgorithm);
//=============================================================================

  if (fHistMgr) { delete fHistMgr; fHistMgr = nullptr; }
  fHistMgr = new THistManager(GetName());

  CreateHistoJets();
  CreateHistoSubjets();
  CreateHistoJetConstis();
//=============================================================================

  TObject *p(nullptr);
  TIter next(fHistMgr->GetListOfHistograms());
  while ((p = next())) fOutput->Add(p);
//=============================================================================

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalSubjet::Run()
{
//
//  Bool_t AliAnalysisTaskEmcalSubjet::Run()
//

  if (!AliAnalysisTaskEmcalJet::Run()) return kFALSE;
//=============================================================================

  TIter next(&fJetCollArray);
  AliJetContainer *pc(nullptr);
  while ((pc = static_cast<AliJetContainer*>(next()))) LoopJets(pc);
//=============================================================================

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalSubjet::LoopJets(AliJetContainer const *pc)
{
//
//  void AliAnalysisTaskEmcalSubjet::LoopJets(AliJetContainer const *pc)
//
  if (!pc) return;
//=============================================================================

  const TString s(pc->GetName());
  const auto d(pc->GetRhoName().IsNull() ? -1. : pc->GetRhoVal());
//=============================================================================

  UInt_t w(0);
  for (auto pj : pc->accepted()) if (!(pj->IsGhost())) {

    w += 1;
    const auto da(pj->Area());
    const auto dj(d>0. ? pj->Pt() - d*da : pj->Pt());
    fHistMgr->FillTH1(Form("%s/hJetPt_%d",s.Data(),fCentBin), dj);
    fHistMgr->FillTH2(Form("%s/hArea_%d", s.Data(),fCentBin), dj, da);
//=============================================================================

    LoopJetConstis(pj, pc);
    LoopSubjets(pj, pc);
  }
//=============================================================================

  fHistMgr->FillTH1(Form("%s/hNjets_%d",s.Data(),fCentBin), w);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalSubjet::LoopSubjets(AliEmcalJet const *pj, AliJetContainer const *pc)
{
//
//  void AliAnalysisTaskEmcalSubjet::LoopSubjets(AliEmcalJet const *pj, AliJetContainer const *pc)
//

  if (!(pj && pc)) return;
  if (!(fSubjetFinder->Filter(const_cast<AliEmcalJet*>(pj),const_cast<AliJetContainer*>(pc),fVertex))) return;
//=============================================================================

  const TString s(Form("%s_sj",pc->GetName()));
  const auto d(pc->GetRhoName().IsNull() ? -1. : pc->GetRhoVal());
  const auto dj((d>0.) ? (pj->Pt() - d*(pj->Area())) : pj->Pt());
//=============================================================================

  UInt_t w(0);
  auto d1(-9999.), d2(-9999.);
  for (auto i=0; i<fSubjetFinder->GetNumberOfJets(); ++i) {
    auto p(fSubjetFinder->GetJet(i)); if (!p) continue;
//  if (p->IsGhost()) continue;

    w += 1;
    const auto da(p->Area());
    const auto dd(d>0. ? (p->Pt() - d*da) : p->Pt());
    if (dd>d1) { d2 = d1; d1 = dd; } else if (dd>d2) { d2 = dd; }
//=============================================================================

    fHistMgr->FillTH2(Form("%s/hSubjetPt_%d",  s.Data(),fCentBin), dd, dj);
    fHistMgr->FillTH2(Form("%s/hSubjetArea_%d",s.Data(),fCentBin), dd, da);
  }
//=============================================================================

  fHistMgr->FillTH1(Form("%s/hNsubjets_%d",s.Data(),fCentBin), w);
  fHistMgr->FillTH2(Form("%s/hL1stSjPt_%d",s.Data(),fCentBin), d1, dj);
  fHistMgr->FillTH2(Form("%s/hL2ndSjPt_%d",s.Data(),fCentBin), d2, dj);
  fHistMgr->FillTH2(Form("%s/hDelta_Sj_%d",s.Data(),fCentBin), d1 - d2, dj);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalSubjet::LoopJetConstis(AliEmcalJet const *pj, AliJetContainer const *pc)
{
//
//  void AliAnalysisTaskEmcalSubjet::LoopJetConstis(AliEmcalJet const *pj, AliJetContainer const *pc)
//

  if (!(pj && pc)) return;
//=============================================================================

  const TString s(Form("%s_consti",pc->GetName()));
  const auto d(pc->GetRhoName().IsNull() ? -1. : pc->GetRhoVal());
  const auto dj((d>0.) ? (pj->Pt() - d*(pj->Area())) : pj->Pt());
//=============================================================================

  UInt_t w(0);
  auto d1(-9999.), d2(-9999.);
  for (auto i=0; i<pj->GetNumberOfTracks(); ++i) {
    const auto p(pj->Track(i)); if (!p) continue;

    w += 1;
    const auto dd(p->Pt());
    if (dd>d1) { d2 = d1; d1 = dd; } else if (dd>d2) { d2 = dd; }
  }
//=============================================================================

/*TLorentzVector v;
  for (auto i=0; i<pj->GetNumberOfClusters(); ++i) {
    auto p(pj->Cluster(i)); if (!p) continue;
    p->GetMomentum(v,fVertex);

    w += 1;
    auto dd(v.Pt());
    if (dd>d1) { d2 = d1; d1 = dd; } else if (dd>d2) { d2 = dd; }
  }*/
//=============================================================================

  fHistMgr->FillTH1(Form("%s/hNtracks_%d",s.Data(),fCentBin), w);
  fHistMgr->FillTH2(Form("%s/hL1stTrk_%d",s.Data(),fCentBin), d1, dj);
  fHistMgr->FillTH2(Form("%s/hL2ndTrk_%d",s.Data(),fCentBin), d2, dj);
  fHistMgr->FillTH2(Form("%s/hDeltaTk_%d",s.Data(),fCentBin), d1-d2, dj);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalSubjet::CreateHistoJets()
{
//
//  void AliAnalysisTaskEmcalSubjet::CreateHistoJets()
//

  TIter next(&fJetCollArray);
  AliJetContainer *p(nullptr);
  while ((p = static_cast<AliJetContainer*>(next()))) {
    const TString s(p->GetName());

    if (fHistMgr->FindObject(s)) {
       AliWarning(Form("Found group name %s in hist manager", s.Data()));
       continue;
    }
//=============================================================================

    fHistMgr->CreateHistoGroup(s);
    const auto bRho(!(p->GetRhoName().IsNull()));
    const auto dMin(bRho ? -fMaxBinPt/2. : fMinBinPt);
    const auto dMax(bRho ?  fMaxBinPt/2. : fMaxBinPt);
//=============================================================================

    const auto nbs(2*fNbins);
    for (auto i=0; i<fNcentBins; ++i) {
      fHistMgr->CreateTH1(Form("%s/hNjets_%d",s.Data(),i), "", 500,  -0.5,  499.5, "s");
      fHistMgr->CreateTH1(Form("%s/hJetPt_%d",s.Data(),i), "", fNbins, dMin, dMax, "s");
      fHistMgr->CreateTH2(Form("%s/hArea_%d", s.Data(),i), "", fNbins, dMin, dMax, nbs, 0., 3., "s");
    }
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalSubjet::CreateHistoSubjets()
{
//
//  void AliAnalysisTaskEmcalSubjet::CreateHistoSubjets()
//

  TIter next(&fJetCollArray);
  AliJetContainer *p(nullptr);
  while ((p = static_cast<AliJetContainer*>(next()))) {
    const TString s(Form("%s_sj",p->GetName()));

    if (fHistMgr->FindObject(s)) {
       AliWarning(Form("Found group name %s in hist manager", s.Data()));
      continue;
    }
//=============================================================================

    fHistMgr->CreateHistoGroup(s);
    const auto bRho(!(p->GetRhoName().IsNull()));
    const auto dMin(bRho ? -fMaxBinPt/2. : fMinBinPt);
    const auto dMax(bRho ?  fMaxBinPt/2. : fMaxBinPt);
//=============================================================================

    const auto nbs(2*fNbins);
    for (auto i=0; i<fNcentBins; ++i) {
      fHistMgr->CreateTH1(Form("%s/hNsubjets_%d",  s.Data(),i), "", 500, -0.5, 499.5, "s");
      fHistMgr->CreateTH2(Form("%s/hSubjetPt_%d",  s.Data(),i), "", nbs, dMin, dMax, fNbins, dMin, dMax, "s");
      fHistMgr->CreateTH2(Form("%s/hL1stSjPt_%d",  s.Data(),i), "", nbs, dMin, dMax, fNbins, dMin, dMax, "s");
      fHistMgr->CreateTH2(Form("%s/hL2ndSjPt_%d",  s.Data(),i), "", nbs, dMin, dMax, fNbins, dMin, dMax, "s");
      fHistMgr->CreateTH2(Form("%s/hDelta_Sj_%d",  s.Data(),i), "", nbs, dMin, dMax, fNbins, dMin, dMax, "s");
      fHistMgr->CreateTH2(Form("%s/hSubjetArea_%d",s.Data(),i), "", nbs, dMin, dMax, nbs,  0., 3., "s");
    }
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalSubjet::CreateHistoJetConstis()
{
//
//  void AliAnalysisTaskEmcalSubjet::CreateHistoJetConstis()
//

  TIter next(&fJetCollArray);
  AliJetContainer *p(nullptr);
  while ((p = static_cast<AliJetContainer*>(next()))) {
    const TString s(Form("%s_consti",p->GetName()));

    if (fHistMgr->FindObject(s)) {
       AliWarning(Form("Found group name %s in hist manager", s.Data()));
      continue;
    }
//=============================================================================

    fHistMgr->CreateHistoGroup(s);
    const auto bRho(!(p->GetRhoName().IsNull()));
    const auto dMin(bRho ? -fMaxBinPt/2. : fMinBinPt);
    const auto dMax(bRho ?  fMaxBinPt/2. : fMaxBinPt);
//=============================================================================

    const auto nbs(2*fNbins);
    for (auto i=0; i<fNcentBins; ++i) {
      fHistMgr->CreateTH1(Form("%s/hNtracks_%d",s.Data(),i), "", 500, -0.5, 499.5, "s");
      fHistMgr->CreateTH2(Form("%s/hL1stTrk_%d",s.Data(),i), "", nbs, dMin, dMax, fNbins, dMin, dMax, "s");
      fHistMgr->CreateTH2(Form("%s/hL2ndTrk_%d",s.Data(),i), "", nbs, dMin, dMax, fNbins, dMin, dMax, "s");
      fHistMgr->CreateTH2(Form("%s/hDeltaTk_%d",s.Data(),i), "", nbs, dMin, dMax, fNbins, dMin, dMax, "s");
    }
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalSubjet *AliAnalysisTaskEmcalSubjet::AddTask(const TString sTrks,
                                                                const TString sClus,
                                                                const TString sCells)
{
//
//  AliAnalysisTaskEmcalSubjet *AliAnalysisTaskEmcalSubjet::AddTask()
//

  auto mgr(AliAnalysisManager::GetAnalysisManager());

  if (!mgr) {
    ::Error(Form("AliAnalysisTaskEmcalSubjet::%s",__func__), "No analysis manager to connect to");
    return nullptr;
  }
//=============================================================================

  AliVEventHandler *pH(mgr->GetInputEventHandler());

  if (!pH) {
    ::Error(Form("AliAnalysisTaskEmcalSubjet::%s",__func__), "This task requires an input event handler");
    return nullptr;
  }
//=============================================================================

  enum EDataType_t { kUnknown, kESD, kAOD };

  EDataType_t wType(kUnknown);
  if (pH->InheritsFrom("AliESDInputHandler")) wType = kESD;
  if (pH->InheritsFrom("AliAODInputHandler")) wType = kAOD;

  if (wType==kUnknown) {
    ::Error(Form("AliAnalysisTaskEmcalSubjet::%s",__func__), "Unkown data input");
    return nullptr;
  }
//=============================================================================

  auto task(new AliAnalysisTaskEmcalSubjet("AliAnalysisTaskEmcalSubjet",kTRUE));
//=============================================================================

  TString sTrkName(sTrks);
  if (sTrks=="usedefault") {
    if (wType==kESD) sTrkName = "Tracks";
    if (wType==kAOD) sTrkName = "tracks";
  }

  if (sTrkName=="mcparticles") {
    task->AddMCParticleContainer(sTrkName);
  } else if ((sTrkName=="tracks") || (sTrkName=="Tracks")) {
    task->AddTrackContainer(sTrkName);
  } else if (!sTrkName.IsNull()) {
    task->AddParticleContainer(sTrkName);
  }
//=============================================================================

  TString sClsName(sClus);
  if (sClus=="usedefault") {
    if (wType==kESD) sClsName = "CaloClusters";
    if (wType==kAOD) sClsName = "caloClusters";
  }

  task->AddClusterContainer(sClsName);
//=============================================================================

  TString sCellName(sCells);
  if (sCells=="usedefault") {
    if (wType==kESD) sCellName = "EMCALCells";
    if (wType==kAOD) sCellName = "emcalCells";
  }

  task->SetCaloCellsName(sCellName);
//=============================================================================

  mgr->AddTask(task);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("listSubjet",
                                                   AliEmcalList::Class(),
                                                   AliAnalysisManager::kOutputContainer,
                                                   AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
