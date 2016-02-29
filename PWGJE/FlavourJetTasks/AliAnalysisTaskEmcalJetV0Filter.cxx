#include <TH1.h>
#include <TString.h>
#include <TClonesArray.h>

#include "AliAnalysisManager.h"

#include "AliAODHandler.h"
#include "AliAODEvent.h"

#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliPicoHeaderCJ.h"
#include "AliPicoV0RD.h"
#include "AliPicoV0MC.h"
#include "AliPicoJet.h"
#include "AliAnalysisTaskEmcalJetV0Filter.h"

ClassImp(AliAnalysisTaskEmcalJetV0Filter)

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0Filter::AliAnalysisTaskEmcalJetV0Filter() :
AliAnalysisTaskEmcalJet(),
fNameJetRD02(""),
fNameJetRD03(""),
fNameJetRD04(""),
fNameJetMC02(""),
fNameJetMC03(""),
fNameJetMC04(""),
fIsAnaPicoV0(kFALSE),
fAnaPicoV0MC(kFALSE),
fJetsContRD02(0),
fTracksContRD02(0),
fCaloClustersContRD02(0),
fJetsContRD03(0),
fTracksContRD03(0),
fCaloClustersContRD03(0),
fJetsContRD04(0),
fTracksContRD04(0),
fCaloClustersContRD04(0),
fJetsContMC02(0),
fTracksContMC02(0),
fJetsContMC03(0),
fTracksContMC03(0),
fJetsContMC04(0),
fTracksContMC04(0),
fV0s(0),
fPicoHeaderCJ(0),
fPicoJetsClArrRD02(0),
fPicoJetsClArrRD03(0),
fPicoJetsClArrRD04(0),
fPicoJetsClArrMC02(0),
fPicoJetsClArrMC03(0),
fPicoJetsClArrMC04(0),
fPicoV0sClArr(0),
fListUserOutputs(0)
{
//
//  AliAnalysisTaskEmcalJetV0Filter::AliAnalysisTaskEmcalJetV0Filter
//
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0Filter::AliAnalysisTaskEmcalJetV0Filter(const char *name, Bool_t bHistos) :
AliAnalysisTaskEmcalJet(name,bHistos),
fNameJetRD02(""),
fNameJetRD03(""),
fNameJetRD04(""),
fNameJetMC02(""),
fNameJetMC03(""),
fNameJetMC04(""),
fIsAnaPicoV0(kFALSE),
fAnaPicoV0MC(kFALSE),
fJetsContRD02(0),
fTracksContRD02(0),
fCaloClustersContRD02(0),
fJetsContRD03(0),
fTracksContRD03(0),
fCaloClustersContRD03(0),
fJetsContRD04(0),
fTracksContRD04(0),
fCaloClustersContRD04(0),
fJetsContMC02(0),
fTracksContMC02(0),
fJetsContMC03(0),
fTracksContMC03(0),
fJetsContMC04(0),
fTracksContMC04(0),
fV0s(0),
fPicoHeaderCJ(0),
fPicoJetsClArrRD02(0),
fPicoJetsClArrRD03(0),
fPicoJetsClArrRD04(0),
fPicoJetsClArrMC02(0),
fPicoJetsClArrMC03(0),
fPicoJetsClArrMC04(0),
fPicoV0sClArr(0),
fListUserOutputs(0)
{
//
//  AliAnalysisTaskEmcalJetV0Filter::AliAnalysisTaskEmcalJetV0Filter
//

  AliAnalysisTaskEmcal::fGeneralHistograms = bHistos;

//DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0Filter::~AliAnalysisTaskEmcalJetV0Filter()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::~AliAnalysisTaskEmcalJetV0Filter
//

  if (fJetsContRD02)         { delete fJetsContRD02;         fJetsContRD02         = 0; }
  if (fTracksContRD02)       { delete fTracksContRD02;       fTracksContRD02       = 0; }
  if (fCaloClustersContRD02) { delete fCaloClustersContRD02; fCaloClustersContRD02 = 0; }

  if (fJetsContRD03)         { delete fJetsContRD03;         fJetsContRD03         = 0; }
  if (fTracksContRD03)       { delete fTracksContRD03;       fTracksContRD03       = 0; }
  if (fCaloClustersContRD03) { delete fCaloClustersContRD03; fCaloClustersContRD03 = 0; }

  if (fJetsContRD04)         { delete fJetsContRD04;         fJetsContRD04         = 0; }
  if (fTracksContRD04)       { delete fTracksContRD04;       fTracksContRD04       = 0; }
  if (fCaloClustersContRD04) { delete fCaloClustersContRD04; fCaloClustersContRD04 = 0; }

  if (fJetsContMC02)   { delete fJetsContMC02;   fJetsContMC02   = 0; }
  if (fTracksContMC02) { delete fTracksContMC02; fTracksContMC02 = 0; }

  if (fJetsContMC03)   { delete fJetsContMC03;   fJetsContMC03   = 0; }
  if (fTracksContMC03) { delete fTracksContMC03; fTracksContMC03 = 0; }

  if (fJetsContMC04)   { delete fJetsContMC04;   fJetsContMC04   = 0; }
  if (fTracksContMC04) { delete fTracksContMC04; fTracksContMC04 = 0; }

  if (fV0s) { delete fV0s; fV0s = 0; }

  if (fPicoHeaderCJ) { delete fPicoHeaderCJ; fPicoHeaderCJ = 0; }

  if (fPicoJetsClArrRD02) { delete fPicoJetsClArrRD02; fPicoJetsClArrRD02 = 0; }
  if (fPicoJetsClArrRD03) { delete fPicoJetsClArrRD03; fPicoJetsClArrRD03 = 0; }
  if (fPicoJetsClArrRD04) { delete fPicoJetsClArrRD04; fPicoJetsClArrRD04 = 0; }

  if (fPicoJetsClArrMC02) { delete fPicoJetsClArrMC02; fPicoJetsClArrMC02 = 0; }
  if (fPicoJetsClArrMC03) { delete fPicoJetsClArrMC03; fPicoJetsClArrMC03 = 0; }
  if (fPicoJetsClArrMC04) { delete fPicoJetsClArrMC04; fPicoJetsClArrMC04 = 0; }

  if (fPicoV0sClArr)    { delete fPicoV0sClArr;    fPicoV0sClArr    = 0; }

  if (fListUserOutputs) { delete fListUserOutputs; fListUserOutputs = 0; }
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

  if (!fNameJetRD02.IsNull()) {
    fJetsContRD02 = GetJetContainer(fNameJetRD02.Data());

    if (fJetsContRD02) {
      fTracksContRD02       = fJetsContRD02->GetParticleContainer();
      fCaloClustersContRD02 = fJetsContRD02->GetClusterContainer();
    }
  }

  if (!fNameJetRD03.IsNull()) {
    fJetsContRD03 = GetJetContainer(fNameJetRD03.Data());

    if (fJetsContRD03) {
      fTracksContRD03       = fJetsContRD03->GetParticleContainer();
      fCaloClustersContRD03 = fJetsContRD03->GetClusterContainer();
    }
  }

  if (!fNameJetRD04.IsNull()) {
    fJetsContRD04 = GetJetContainer(fNameJetRD04.Data());

    if (fJetsContRD04) {
      fTracksContRD04       = fJetsContRD04->GetParticleContainer();
      fCaloClustersContRD04 = fJetsContRD04->GetClusterContainer();
    }
  }
//=============================================================================

  if (!fNameJetMC02.IsNull()) {
    fJetsContMC02 = GetJetContainer(fNameJetMC02.Data());
    if (fJetsContMC02) fTracksContMC02 = fJetsContMC02->GetParticleContainer();
  }

  if (!fNameJetMC03.IsNull()) {
    fJetsContMC03 = GetJetContainer(fNameJetMC03.Data());
    if (fJetsContMC03) fTracksContMC03 = fJetsContMC03->GetParticleContainer();
  }

  if (!fNameJetMC04.IsNull()) {
    fJetsContMC04 = GetJetContainer(fNameJetMC04.Data());
    if (fJetsContMC04) fTracksContMC04 = fJetsContMC04->GetParticleContainer();
  }
//=============================================================================

  fPicoHeaderCJ = new  AliPicoHeaderCJ();
  fPicoHeaderCJ->SetName("PicoHeaderCJ");
  AddAODBranch("AliPicoHeaderCJ", &fPicoHeaderCJ);

  if (fJetsContRD02) {
    fPicoJetsClArrRD02 = new TClonesArray("AliPicoJet");
    fPicoJetsClArrRD02->SetName("PicoJetsRD02");
    AddAODBranch("TClonesArray", &fPicoJetsClArrRD02);
  }

  if (fJetsContRD03) {
    fPicoJetsClArrRD03 = new TClonesArray("AliPicoJet");
    fPicoJetsClArrRD03->SetName("PicoJetsRD03");
    AddAODBranch("TClonesArray", &fPicoJetsClArrRD03);
  }

  if (fJetsContRD04) {
    fPicoJetsClArrRD04 = new TClonesArray("AliPicoJet");
    fPicoJetsClArrRD04->SetName("PicoJetsRD04");
    AddAODBranch("TClonesArray", &fPicoJetsClArrRD04);
  }
//=============================================================================

  if (fJetsContMC02) {
    fPicoJetsClArrMC02 = new TClonesArray("AliPicoJet");
    fPicoJetsClArrMC02->SetName("PicoJetsMC02");
    AddAODBranch("TClonesArray", &fPicoJetsClArrMC02);
  }

  if (fJetsContMC03) {
    fPicoJetsClArrMC03 = new TClonesArray("AliPicoJet");
    fPicoJetsClArrMC03->SetName("PicoJetsMC03");
    AddAODBranch("TClonesArray", &fPicoJetsClArrMC03);
  }

  if (fJetsContMC04) {
    fPicoJetsClArrMC04 = new TClonesArray("AliPicoJet");
    fPicoJetsClArrMC04->SetName("PicoJetsMC04");
    AddAODBranch("TClonesArray", &fPicoJetsClArrMC04);
  }
//=============================================================================

  if (fIsAnaPicoV0) {
    if (fAnaPicoV0MC) {
      fPicoV0sClArr = new TClonesArray("AliPicoV0MC");
      fPicoV0sClArr->SetName("PicoV0sMC");
    } else {
      fPicoV0sClArr = new TClonesArray("AliPicoV0RD");
      fPicoV0sClArr->SetName("PicoV0sRD");
    }

    AddAODBranch("TClonesArray", &fPicoV0sClArr);
  }

/*fListUserOutputs = new TList();
  fListUserOutputs->SetOwner();
  CreateUserOutputHistograms();
  PostData(2, fListUserOutputs);*/
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

  Int_t ncs = 0;
  fPicoHeaderCJ->Reset();

  if (fPicoJetsClArrRD02) fPicoJetsClArrRD02->Delete();
  if (fPicoJetsClArrRD03) fPicoJetsClArrRD03->Delete();
  if (fPicoJetsClArrRD04) fPicoJetsClArrRD04->Delete();

  if (fPicoJetsClArrMC02) fPicoJetsClArrMC02->Delete();
  if (fPicoJetsClArrMC03) fPicoJetsClArrMC03->Delete();
  if (fPicoJetsClArrMC04) fPicoJetsClArrMC04->Delete();

  if (fPicoV0sClArr) fPicoV0sClArr->Delete();
  else return kFALSE; // Should not happen, make Coverity happy

  if (!AliAnalysisTaskEmcalJet::Run()) return kFALSE;
//=============================================================================

  fPicoHeaderCJ->SetEventInfo(fInputHandler);
//=============================================================================

  if (fJetsContRD02) {
    ncs = fPicoJetsClArrRD02->GetEntriesFast();
    fPicoHeaderCJ->BackgroundRhoRD02(fJetsContRD02->GetRhoVal());
    fJetsContRD02->ResetCurrentID();
    AliEmcalJet *pJet = fJetsContRD02->GetNextAcceptJet(); while (pJet) {
      new ((*fPicoJetsClArrRD02)[ncs++]) AliPicoJet(pJet, GetLeadingHadronPt(pJet));
      pJet = fJetsContRD02->GetNextAcceptJet();
    }
  }

  if (fJetsContRD03) {
    ncs = fPicoJetsClArrRD03->GetEntriesFast();
    fPicoHeaderCJ->BackgroundRhoRD03(fJetsContRD03->GetRhoVal());
    fJetsContRD03->ResetCurrentID();
    AliEmcalJet *pJet = fJetsContRD03->GetNextAcceptJet(); while (pJet) {
      new ((*fPicoJetsClArrRD03)[ncs++]) AliPicoJet(pJet, GetLeadingHadronPt(pJet));
      pJet = fJetsContRD03->GetNextAcceptJet();
    }
  }

  if (fJetsContRD04) {
    ncs = fPicoJetsClArrRD04->GetEntriesFast();
    fPicoHeaderCJ->BackgroundRhoRD04(fJetsContRD04->GetRhoVal());
    fJetsContRD04->ResetCurrentID();
    AliEmcalJet *pJet = fJetsContRD04->GetNextAcceptJet(); while (pJet) {
      new ((*fPicoJetsClArrRD04)[ncs++]) AliPicoJet(pJet, GetLeadingHadronPt(pJet));
      pJet = fJetsContRD04->GetNextAcceptJet();
    }
  }
//=============================================================================

  if (fJetsContMC02) {
    ncs = fPicoJetsClArrMC02->GetEntriesFast();
    fPicoHeaderCJ->BackgroundRhoMC02(fJetsContMC02->GetRhoVal());
    fJetsContMC02->ResetCurrentID();
    AliEmcalJet *pJet = fJetsContMC02->GetNextAcceptJet(); while (pJet) {
      new ((*fPicoJetsClArrMC02)[ncs++]) AliPicoJet(pJet, GetLeadingHadronPt(pJet));
      pJet = fJetsContMC02->GetNextAcceptJet();
    }
  }

  if (fJetsContMC03) {
    ncs = fPicoJetsClArrMC03->GetEntriesFast();
    fPicoHeaderCJ->BackgroundRhoMC03(fJetsContMC03->GetRhoVal());
    fJetsContMC03->ResetCurrentID();
    AliEmcalJet *pJet = fJetsContMC03->GetNextAcceptJet(); while (pJet) {
      new ((*fPicoJetsClArrMC03)[ncs++]) AliPicoJet(pJet, GetLeadingHadronPt(pJet));
      pJet = fJetsContMC03->GetNextAcceptJet();
    }
  }

  if (fJetsContMC04) {
    ncs = fPicoJetsClArrMC04->GetEntriesFast();
    fPicoHeaderCJ->BackgroundRhoMC04(fJetsContMC04->GetRhoVal());
    fJetsContMC04->ResetCurrentID();
    AliEmcalJet *pJet = fJetsContMC04->GetNextAcceptJet(); while (pJet) {
      new ((*fPicoJetsClArrMC04)[ncs++]) AliPicoJet(pJet, GetLeadingHadronPt(pJet));
      pJet = fJetsContMC04->GetNextAcceptJet();
    }
  }
//=============================================================================

  if (fV0s) {
    AliPicoV0RD *pV0RD = 0;
    AliPicoV0MC *pV0MC = 0;
    ncs = fPicoV0sClArr->GetEntriesFast();
    for (Int_t i=0; i<fV0s->GetEntriesFast(); i++) {

      if (fAnaPicoV0MC) {
        pV0MC = static_cast<AliPicoV0MC*>(fV0s->At(i)); if (!pV0MC) continue;
        new ((*fPicoV0sClArr)[ncs++]) AliPicoV0MC(*pV0MC); pV0MC = 0;
      } else {
        pV0RD = static_cast<AliPicoV0RD*>(fV0s->At(i)); if (!pV0RD) continue;
        new ((*fPicoV0sClArr)[ncs++]) AliPicoV0RD(*pV0RD); pV0RD = 0;
      }
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0Filter::RetrieveEventObjects()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::RetrieveEventObjects
//

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects()) return kFALSE;
  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);

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

  if (fIsAnaPicoV0 && (!fV0s)) {
    fV0s = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("PicoV0s"));

    if (!fV0s) {
      AliError(Form("%s: Could not retrieve V0 %s!", GetName(), "PicoV0s"));
      fInitialized = kFALSE;
      return;
    }
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0Filter::CreateUserOutputHistograms()
{
//
//  AliAnalysisTaskEmcalJetV0Filter::CreateUserOutputHistograms
//

  if (!fListUserOutputs) return;

  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1::AddDirectory(bStatusTmpH);
  return;
}
