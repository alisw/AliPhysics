#include <TObjString.h>
#include <TClonesArray.h>

#include "AliAnalysisManager.h"

#include "AliAODHandler.h"

#include "AliEmcalJet.h"
#include "AliJetContainer.h"

#include "AliPicoV0RD.h"
#include "AliPicoV0MC.h"

#include "AliPicoJet.h"
#include "AliPicoHeaderJet.h"

#include "AliAnalysisTaskEmcalJetV0Filter.h"

ClassImp(AliAnalysisTaskEmcalJetV0Filter)

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0Filter::AliAnalysisTaskEmcalJetV0Filter() :
AliAnalysisTaskEmcalJet(),
fMult(""),
fIsMC(kFALSE),
fV0sName(""),
fV0s(nullptr),
fPicoHeader(nullptr),
fMapJets(),
fPicoV0sClArr(nullptr),
fListUserOutputs(nullptr)
{
//
//  AliAnalysisTaskEmcalJetV0Filter::AliAnalysisTaskEmcalJetV0Filter
//
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0Filter::AliAnalysisTaskEmcalJetV0Filter(const char *name, Bool_t bHistos) :
AliAnalysisTaskEmcalJet(name,bHistos),
fMult(""),
fIsMC(kFALSE),
fV0sName("PicoV0s"),
fV0s(nullptr),
fPicoHeader(nullptr),
fMapJets(),
fPicoV0sClArr(nullptr),
fListUserOutputs(nullptr)
{
//
//  AliAnalysisTaskEmcalJetV0Filter::AliAnalysisTaskEmcalJetV0Filter
//

  AliAnalysisTaskEmcal::fGeneralHistograms = bHistos;

//DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0Filter::~AliAnalysisTaskEmcalJetV0Filter()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::~AliAnalysisTaskEmcalJetV0Filter
//

  if (fV0s) { delete fV0s; fV0s = nullptr; }

  if (fPicoHeader)   { delete fPicoHeader;   fPicoHeader   = nullptr; }
  if (fPicoV0sClArr) { delete fPicoV0sClArr; fPicoV0sClArr = nullptr; }

  if (fListUserOutputs) { delete fListUserOutputs; fListUserOutputs = nullptr; }
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0Filter::Init()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::Init
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0Filter::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::UserCreateOutputObjects
//

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
//=============================================================================

  if (fPicoHeader) {
    delete fPicoHeader;
    fPicoHeader = nullptr;
  }

  fPicoHeader = new AliPicoHeaderJet();
  fPicoHeader->SetName("PicoHeaderJet");

  if (!fMult.IsNull()) {
    const auto aMult(fMult.Tokenize(":"));

    TIter next(aMult);
    TObjString *ps(nullptr);
    while ((ps = static_cast<TObjString*>(next()))) {
      fPicoHeader->AddMultEstimator(ps->String());
    }
  }

  AddAODBranch("AliPicoHeaderJet", &fPicoHeader);
//=============================================================================

  fMapJets.DeleteAll();
  if (fJetCollArray.GetEntriesFast()>0) {
    TIter next(&fJetCollArray);
    AliEmcalContainer *pc(nullptr);
    while ((pc = static_cast<AliEmcalContainer*>(next()))) if (pc) {
      const TString s(pc->GetName());
      auto pa(new TClonesArray("AliPicoJet"));
      pa->SetName(Form("Pico%s",s.Data()));
      fMapJets.Add(new TObjString(s.Data()), pa);
      AddAODBranch("TClonesArray", &pa);
    }
  }
//=============================================================================

  if (!fV0sName.IsNull()) {
    if (fPicoV0sClArr) {
      delete fPicoV0sClArr;
      fPicoV0sClArr = nullptr;
    }

    if (fIsMC) {
      fPicoV0sClArr = new TClonesArray("AliPicoV0MC");
      fPicoV0sClArr->SetName("PicoV0sMC");
    } else {
      fPicoV0sClArr = new TClonesArray("AliPicoV0RD");
      fPicoV0sClArr->SetName("PicoV0sRD");
    }

    AddAODBranch("TClonesArray", &fPicoV0sClArr);
  }
//=============================================================================

/*if (fListUserOutputs) {
    delete fListUserOutputs;
    fListUserOutputs = nullptr;
  }

  fListUserOutputs = new TList();
  fListUserOutputs->SetOwner();

  const auto b(TH1::AddDirectoryStatus());
  TH1::AddDirectory(kFALSE);
//TODO
  TH1::AddDirectory(b);
  PostData(1, fListUserOutputs);*/
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0Filter::Terminate(Option_t *opt)
{
//
//  AliAnalysisTaskEmcalJetV0Filter::Terminate
//

  AliAnalysisTaskEmcalJet::Terminate(opt);

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0Filter::Run()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::Run
//

  fPicoHeader->Reset();

  TObjString *ps(nullptr);
  const auto na(fMapJets.MakeIterator());
  while ((ps = static_cast<TObjString*>((*na)()))) {
    const auto s(ps->String());
    const auto p(static_cast<TClonesArray*>(fMapJets(s.Data())));
    if (p) p->Delete();
  }

  if (fPicoV0sClArr) fPicoV0sClArr->Delete();
//=============================================================================

  if (!AliAnalysisTaskEmcalJet::Run()) return kFALSE;
//=============================================================================

  fPicoHeader->SetEventInfo(fInputHandler);
  if (fRho) fPicoHeader->BackgroundRho(fRhoVal);
//=============================================================================

  TIter nc(&fJetCollArray);
  AliEmcalContainer *pc(nullptr);
  while ((pc = static_cast<AliEmcalContainer*>(nc()))) {
    auto pa(static_cast<TClonesArray*>(fMapJets(pc->GetName()))); if (!pa) continue;

    auto l(pa->GetEntriesFast());
    for (auto p : pc->accepted()) if (p) {
      const auto pJet(static_cast<AliEmcalJet*>(p));
      if (pJet) new ((*pa)[l++]) AliPicoJet(pJet,GetLeadingHadronPt(pJet));
    }
  }
//=============================================================================

  if (fV0s) {
    auto l(fPicoV0sClArr->GetEntriesFast());
    for (Int_t i=0; i<fV0s->GetEntriesFast(); ++i) {
      if (fIsMC) {
        const auto pV0(dynamic_cast<AliPicoV0MC*>(fV0s->At(i)));

        if (!pV0) continue;
        new ((*fPicoV0sClArr)[l++]) AliPicoV0MC(*pV0);
      } else {
        const auto pV0(dynamic_cast<AliPicoV0RD*>(fV0s->At(i)));

        if (!pV0) continue;
        new ((*fPicoV0sClArr)[l++]) AliPicoV0RD(*pV0);
      }
    }
  }
//=============================================================================

  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
//=============================================================================

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0Filter::RetrieveEventObjects()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::RetrieveEventObjects
//

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0Filter::IsEventSelected()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::IsEventSelected
//

  if (!AliAnalysisTaskEmcalJet::IsEventSelected()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0Filter::FillHistograms()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::FillHistograms
//

  if (!AliAnalysisTaskEmcalJet::FillHistograms()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0Filter::FillGeneralHistograms()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::FillGeneralHistograms
//

  if (!AliAnalysisTaskEmcalJet::FillGeneralHistograms()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0Filter::ExecOnce()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::ExecOnce
//

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (!fInitialized) return;
//=============================================================================

  if ((!fV0sName.IsNull()) && (!fV0s)) {
    fV0s = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fV0sName.Data()));

    if (!fV0s) {
      AliError(Form("%s: Could not retrieve V0 %s!", GetName(), fV0sName.Data()));
      fInitialized = kFALSE;
      return;
    }
  }
//=============================================================================

  return;
}
