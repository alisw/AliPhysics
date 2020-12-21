#include <THnSparse.h>
#include <TClonesArray.h>
#include <TParticle.h>

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"

#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliPicoV0MC.h"
#include "AliAnalysisTaskEmcalJetV0CF.h"

ClassImp(AliAnalysisTaskEmcalJetV0CF)

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0CF::AliAnalysisTaskEmcalJetV0CF() :
AliAnalysisTaskEmcalJet(),
fKaCutNS(0.),
fLaCutNS(0.),
fV0CutMinEta(0.),
fV0CutMaxEta(0.),
fEventAOD(nullptr),
fEventESD(nullptr),
fCentInfo(nullptr),
fJetsContRD(nullptr),
fTracksContRD(nullptr),
fCaloClustersContRD(nullptr),
fJetsContMC(nullptr),
fTracksContMC(nullptr),
fV0s(nullptr),
fHistoKshortInvM(nullptr),
fHistoLambdaInvM(nullptr),
fHistoAntiLaInvM(nullptr),
fListUserOutputs(nullptr)
{
//
//  AliAnalysisTaskEmcalJetV0CF::AliAnalysisTaskEmcalJetV0CF
//
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0CF::AliAnalysisTaskEmcalJetV0CF(const char *name, Bool_t bHistos) :
AliAnalysisTaskEmcalJet(name,bHistos),
fKaCutNS(6.),
fLaCutNS(6.),
fV0CutMinEta(-10.),
fV0CutMaxEta(10.),
fEventAOD(nullptr),
fEventESD(nullptr),
fCentInfo(nullptr),
fJetsContRD(nullptr),
fTracksContRD(nullptr),
fCaloClustersContRD(nullptr),
fJetsContMC(nullptr),
fTracksContMC(nullptr),
fV0s(nullptr),
fHistoKshortInvM(nullptr),
fHistoLambdaInvM(nullptr),
fHistoAntiLaInvM(nullptr),
fListUserOutputs(nullptr)
{
//
//  AliAnalysisTaskEmcalJetV0CF::AliAnalysisTaskEmcalJetV0CF
//

  AliAnalysisTaskEmcal::fGeneralHistograms = bHistos;

  DefineOutput(bHistos ? 2 : 1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0CF::~AliAnalysisTaskEmcalJetV0CF()
{
//
//  AliAnalysisTaskEmcalJetV0CF::~AliAnalysisTaskEmcalJetV0CF
//

  if (fEventAOD) { delete fEventAOD; fEventAOD = nullptr; }
  if (fEventESD) { delete fEventESD; fEventESD = nullptr; }
  if (fCentInfo) { delete fCentInfo; fCentInfo = nullptr; }

  if (fJetsContRD)         { delete fJetsContRD;         fJetsContRD         = nullptr; }
  if (fTracksContRD)       { delete fTracksContRD;       fTracksContRD       = nullptr; }
  if (fCaloClustersContRD) { delete fCaloClustersContRD; fCaloClustersContRD = nullptr; }

  if (fJetsContMC)   { delete fJetsContMC;   fJetsContMC   = nullptr; }
  if (fTracksContMC) { delete fTracksContMC; fTracksContMC = nullptr; }

  if (fV0s) { delete fV0s; fV0s = nullptr; }

  if (fHistoKshortInvM) { delete fHistoKshortInvM; fHistoKshortInvM = nullptr; }
  if (fHistoLambdaInvM) { delete fHistoLambdaInvM; fHistoLambdaInvM = nullptr; }
  if (fHistoAntiLaInvM) { delete fHistoAntiLaInvM; fHistoAntiLaInvM = nullptr; }

  if (fListUserOutputs) { delete fListUserOutputs; fListUserOutputs = nullptr; }
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0CF::Init()
{
//
//  AliAnalysisTaskEmcalJetV0CF::Init
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0CF::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskEmcalJetV0CF::UserCreateOutputObjects
//

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
//=============================================================================

  fJetsContRD = GetJetContainer(0);
  if (fJetsContRD) {
    fTracksContRD       = fJetsContRD->GetParticleContainer();
    fCaloClustersContRD = fJetsContRD->GetClusterContainer();
  }

  fJetsContMC = GetJetContainer(1);
  if (fJetsContMC) fTracksContMC = fJetsContMC->GetParticleContainer();
//=============================================================================

  fListUserOutputs = new TList();
  fListUserOutputs->SetOwner();
  CreateUserOutputHistograms();
  PostData(fCreateHisto ? 2 : 1, fListUserOutputs);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0CF::Terminate(Option_t *opt)
{
//
//  AliAnalysisTaskEmcalJetV0CF::Terminate
//

  AliAnalysisTaskEmcalJet::Terminate(opt);

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::Run()
{
//
//  AliAnalysisTaskEmcalJetV0CF::Run
//

  if (!AliAnalysisTaskEmcalJet::Run()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::RetrieveEventObjects()
{
//
//  AliAnalysisTaskEmcalJetV0CF::RetrieveEventObjects
//

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects()) return kFALSE;

  fCentInfo = InputEvent()->GetCentrality(); if (!fCentInfo) return kFALSE;

  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fEventESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if ((!fEventAOD) && (!fEventESD)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::IsEventSelected()
{
//
//  AliAnalysisTaskEmcalJetV0CF::IsEventSelected
//

  if (!AliAnalysisTaskEmcalJet::IsEventSelected()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::FillHistograms()
{
//
//  AliAnalysisTaskEmcalJetV0CF::FillHistograms
//

  if (!AliAnalysisTaskEmcalJet::FillHistograms()) return kFALSE;

  if (FillRecoInfo()) return kFALSE;
  if (FillKineInfo()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::FillGeneralHistograms()
{
//
//  AliAnalysisTaskEmcalJetV0CF::FillGeneralHistograms
//

  if (!AliAnalysisTaskEmcalJet::FillGeneralHistograms()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetV0CF::ExecOnce()
{
//
//  AliAnalysisTaskEmcalJetV0CF::ExecOnce
//

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (!fInitialized) return;

  if (!fV0s) {
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
void AliAnalysisTaskEmcalJetV0CF::CreateUserOutputHistograms()
{
//
//  AliAnalysisTaskEmcalJetV0CF::CreateUserOutputHistograms
//

  if (!fListUserOutputs) return;

  const auto b(TH1::AddDirectoryStatus());
  TH1::AddDirectory(kFALSE);

  const Int_t nV0(8); // 0: particle type
                      //    ==0, Kshort
                      //    ==1, Lambda
                      //    ==2, AntiLa
                      // 1: Jet pT bin
                      //    == 0.5, jet pT>10.
                      //    == 1.5, jet pT>15.
                      //    == 2.5, jet pT>20.
                      //    == 3.5, jet pT>25.
                      // 2: V0M
                      // 3: V0A
                      // 4: CL1
                      // 5: ZNA
                      // 6: eta
                      // 7: Pt
  const Int_t    nV0Bin[nV0] = {  3,   4,  210,  210,  210,  210, 100, 1000  };
  const Double_t dV0Min[nV0] = { -0.5, 0., -10., -10., -10., -10., -5.,   0. };
  const Double_t dV0Max[nV0] = {  2.5, 4., 200., 200., 200., 200.,  5., 100. };
  fListUserOutputs->Add(new THnSparseD("hsReco", "", nV0, nV0Bin, dV0Min, dV0Max));
  fListUserOutputs->Add(new THnSparseD("hsKine", "", nV0, nV0Bin, dV0Min, dV0Max));

  TH1::AddDirectory(b);
  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::FillRecoInfo()
{
//
//  AliAnalysisTaskEmcalJetV0CF::FillRecoInfo
//

  if (!fV0s) return kTRUE;
//=============================================================================

  const auto dV0M(fCentInfo->GetCentralityPercentile("V0M"));
  const auto dV0A(fCentInfo->GetCentralityPercentile("V0A"));
  const auto dCL1(fCentInfo->GetCentralityPercentile("CL1"));
  const auto dZNA(fCentInfo->GetCentralityPercentile("ZNA"));

  auto hs(static_cast<THnSparseD*>(fListUserOutputs->FindObject("hsReco")));
//=============================================================================

  for (auto i=0; i<fV0s->GetEntriesFast(); ++i) {
    const auto pV0(static_cast<AliPicoV0MC*>(fV0s->At(i))); if (!pV0) continue;

    if (!pV0->IsV0PhysicalPrimary()) continue;
    if (!pV0->IsV0InEtaAcc(fV0CutMinEta,fV0CutMaxEta)) continue;

    const auto bKshort(pV0->IsKshortMC());
    const auto bLambda(pV0->IsLambdaMC());
    const auto bAntiLa(pV0->IsAntiLaMC());
    if (!(bKshort || bLambda || bAntiLa)) continue;
//=============================================================================

    const auto vV0(pV0->KineMC().Vect());
    const auto dPt(vV0.Pt());

    auto histo(fHistoKshortInvM);
    if (bLambda) histo = fHistoLambdaInvM;
    if (bAntiLa) histo = fHistoAntiLaInvM;

    const auto k(histo->FindBin(dPt)); if (k<=0) continue;
    const auto dMean(histo->GetBinContent(k));
    const auto dSigma(histo->GetBinError(k));

    const auto dCutNS(bKshort ? fKaCutNS : fLaCutNS);
    const auto dUpperL(dMean - (dCutNS * dSigma));
    const auto dLowerR(dMean + (dCutNS * dSigma));

    auto dInvM(pV0->KineKshort().M());
    if (bLambda) dInvM = pV0->KineLambda().M();
    if (bAntiLa) dInvM = pV0->KineAntiLa().M();
    if ((dInvM<dUpperL) || (dInvM>=dLowerR)) continue;
//=============================================================================

    Double_t dVar[8];
    if (bKshort) dVar[0] = 0.;
    if (bLambda) dVar[0] = 1.;
    if (bAntiLa) dVar[0] = 2.;

    dVar[1] = -1.;
    if (IsV0InJet(vV0,10.)) dVar[1] = 0.5;
    if (IsV0InJet(vV0,15.)) dVar[1] = 1.5;
    if (IsV0InJet(vV0,20.)) dVar[1] = 2.5;
    if (IsV0InJet(vV0,25.)) dVar[1] = 3.5;
    if (dVar[1]<0.) continue;

    dVar[2] = dV0M;
    dVar[3] = dV0A;
    dVar[4] = dCL1;
    dVar[5] = dZNA;

    dVar[6] = vV0.Eta();
    dVar[7] = dPt;

    hs->Fill(dVar);
  }
//=============================================================================

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::FillKineInfo()
{
//
//  AliAnalysisTaskEmcalJetV0CF::FillKineInfo
//

  const auto dV0M(fCentInfo->GetCentralityPercentile("V0M"));
  const auto dV0A(fCentInfo->GetCentralityPercentile("V0A"));
  const auto dCL1(fCentInfo->GetCentralityPercentile("CL1"));
  const auto dZNA(fCentInfo->GetCentralityPercentile("ZNA"));

  auto hs(dynamic_cast<THnSparseD*>(fListUserOutputs->FindObject("hsKine")));
//=============================================================================

  AliStack *pStack(nullptr);

  if (fEventESD) {
    pStack = MCEvent()->Stack();
    if (!pStack) return kTRUE;
  }
//=============================================================================

  for (Int_t i=0; i<MCEvent()->GetNumberOfTracks(); ++i) {
    TParticle *pESD(nullptr);
    AliMCParticle *pTmp(nullptr);
    AliAODMCParticle *pAOD(nullptr);

    if (fEventAOD) {
      pAOD = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(i));
      if (!pAOD) continue;
    }

    if (fEventESD) {
      pTmp = static_cast<AliMCParticle*>(MCEvent()->GetTrack(i));
      if (!pTmp) continue;

      pESD = pTmp->Particle();
      if (!pESD) continue;
    }
//=============================================================================

    const auto bPhy(pAOD ? pAOD->IsPhysicalPrimary() : pStack->IsPhysicalPrimary(i));
    if (bPhy) continue;

    const auto dEta(pAOD ? pAOD->Eta() : pESD->Eta());
    if ((dEta<fV0CutMinEta) || (dEta>=fV0CutMaxEta)) continue;
//=============================================================================

    Double_t dVar[8];
    const auto id(pAOD ? pAOD->GetPdgCode() : pESD->GetPdgCode());

    dVar[0] = -1.;
    if (id== 310 ) dVar[0] = 0.;
    if (id== 3122) dVar[0] = 1.;
    if (id==-3122) dVar[0] = 2.;
    if (dVar[0]<-0.5) continue;
//=============================================================================

    TVector3 vV0;
    if (pAOD) vV0.SetXYZ(pAOD->Px(), pAOD->Py(), pAOD->Pz());
    if (pESD) vV0.SetXYZ(pESD->Px(), pESD->Py(), pESD->Pz());

    dVar[1] = -1.;
    if (IsV0InJet(vV0,10.)) dVar[1] = 0.5;
    if (IsV0InJet(vV0,15.)) dVar[1] = 1.5;
    if (IsV0InJet(vV0,20.)) dVar[1] = 2.5;
    if (IsV0InJet(vV0,25.)) dVar[1] = 3.5;
    if (dVar[1]<0.) continue;
//=============================================================================

    dVar[2] = dV0M;
    dVar[3] = dV0A;
    dVar[4] = dCL1;
    dVar[5] = dZNA;
    dVar[6] = dEta;
    if (pAOD) dVar[7] = pAOD->Pt();
    if (pESD) dVar[7] = pESD->Pt();
    hs->Fill(dVar);
  }
//=============================================================================

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::IsV0InJet(TVector3 vV0, Double_t dJetPtMin)
{
//
//  AliAnalysisTaskEmcalJetV0CF::IsV0InJet
//

  if (!fJetsContRD) return kFALSE;
//=============================================================================

  TVector3 vJet;
  const auto dJetRadius(fJetsContRD->GetJetRadius());

  fJetsContRD->ResetCurrentID();
  auto pJet(fJetsContRD->GetNextAcceptJet()); while (pJet) {
    const auto dPt(fJetsContRD->GetJetPtCorr(fJetsContRD->GetCurrentID()));
    if (dPt<dJetPtMin) { pJet = fJetsContRD->GetNextAcceptJet(); continue; }

    vJet.SetPtEtaPhi(dPt, pJet->Eta(), pJet->Phi());
    if (vJet.DeltaR(vV0)<dJetRadius) return kTRUE;
    pJet = fJetsContRD->GetNextAcceptJet();
  }
//=============================================================================

  return kFALSE;
}
