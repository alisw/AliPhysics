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
fEventAOD(0),
fEventESD(0),
fCentInfo(0),
fJetsContRD(0),
fTracksContRD(0),
fCaloClustersContRD(0),
fJetsContMC(0),
fTracksContMC(0),
fV0s(0),
fHistoKshortInvM(0),
fHistoLambdaInvM(0),
fHistoAntiLaInvM(0),
fListUserOutputs(0)
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
fEventAOD(0),
fEventESD(0),
fCentInfo(0),
fJetsContRD(0),
fTracksContRD(0),
fCaloClustersContRD(0),
fJetsContMC(0),
fTracksContMC(0),
fV0s(0),
fHistoKshortInvM(0),
fHistoLambdaInvM(0),
fHistoAntiLaInvM(0),
fListUserOutputs(0)
{
//
//  AliAnalysisTaskEmcalJetV0CF::AliAnalysisTaskEmcalJetV0CF
//

  AliAnalysisTaskEmcal::fGeneralHistograms = bHistos;

  DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetV0CF::~AliAnalysisTaskEmcalJetV0CF()
{
//
//  AliAnalysisTaskEmcalJetV0CF::~AliAnalysisTaskEmcalJetV0CF
//

  if (fEventAOD) { delete fEventAOD; fEventAOD = 0; }
  if (fEventESD) { delete fEventESD; fEventESD = 0; }
  if (fCentInfo) { delete fCentInfo; fCentInfo = 0; }

  if (fJetsContRD)         { delete fJetsContRD;         fJetsContRD         = 0; }
  if (fTracksContRD)       { delete fTracksContRD;       fTracksContRD       = 0; }
  if (fCaloClustersContRD) { delete fCaloClustersContRD; fCaloClustersContRD = 0; }

  if (fJetsContMC)   { delete fJetsContMC;   fJetsContMC   = 0; }
  if (fTracksContMC) { delete fTracksContMC; fTracksContMC = 0; }

  if (fV0s) { delete fV0s; fV0s = 0; }

  if (fHistoKshortInvM) { delete fHistoKshortInvM; fHistoKshortInvM = 0; }
  if (fHistoLambdaInvM) { delete fHistoLambdaInvM; fHistoLambdaInvM = 0; }
  if (fHistoAntiLaInvM) { delete fHistoAntiLaInvM; fHistoAntiLaInvM = 0; }

  if (fListUserOutputs) { delete fListUserOutputs; fListUserOutputs = 0; }
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
  PostData(2, fListUserOutputs);
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

  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nV0 = 8; // 0: particle type
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

  THnSparseD *hsReco = new THnSparseD("hsReco", "", nV0, nV0Bin, dV0Min, dV0Max); fListUserOutputs->Add(hsReco);
  THnSparseD *hsKine = new THnSparseD("hsKine", "", nV0, nV0Bin, dV0Min, dV0Max); fListUserOutputs->Add(hsKine);

  TH1::AddDirectory(bStatusTmpH);
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

  Double_t dV0M = fCentInfo->GetCentralityPercentile("V0M");
  Double_t dV0A = fCentInfo->GetCentralityPercentile("V0A");
  Double_t dCL1 = fCentInfo->GetCentralityPercentile("CL1");
  Double_t dZNA = fCentInfo->GetCentralityPercentile("ZNA");

  THnSparseD *hs = dynamic_cast<THnSparseD*>(fListUserOutputs->FindObject("hsReco"));
  if (!hs)
    return kTRUE;  // should not happen; make Coverity happen
//=============================================================================

  AliPicoV0MC *pV0 = 0;
  for (Int_t i=0; i<fV0s->GetEntriesFast(); i++) {
    pV0 = static_cast<AliPicoV0MC*>(fV0s->At(i)); if (!pV0) continue;

    if (!pV0->IsV0PhysicalPrimary()) { pV0 = 0; continue; }
    if (!pV0->IsV0InEtaAcc(fV0CutMinEta,fV0CutMaxEta)) { pV0 = 0; continue; }

    TVector3 vV0 = pV0->KineMC().Vect();
    Double_t dPt = vV0.Pt();
    Double_t dVar[8];

    dVar[0] = -1.;
    if (pV0->IsKshort()) {
      Int_t k = fHistoKshortInvM->FindBin(dPt); if (k<=0) { pV0 = 0; continue; }

      Double_t dMean  = fHistoKshortInvM->GetBinContent(k);
      Double_t dSigma = fHistoKshortInvM->GetBinError(k);

      Double_t dUpperL = dMean - (fKaCutNS * dSigma);
      Double_t dLowerR = dMean + (fKaCutNS * dSigma);

      Double_t dInvM = pV0->KineKshort().M();
      if ((dInvM<dUpperL) || (dInvM>=dLowerR)) { pV0 = 0; continue; }
      
      dVar[0] = 0.;
    }

    if (pV0->IsLambda()) {
      Int_t k = fHistoLambdaInvM->FindBin(dPt); if (k<=0) { pV0 = 0; continue; }

      Double_t dMean  = fHistoLambdaInvM->GetBinContent(k);
      Double_t dSigma = fHistoLambdaInvM->GetBinError(k);

      Double_t dUpperL = dMean - (fLaCutNS * dSigma);
      Double_t dLowerR = dMean + (fLaCutNS * dSigma);

      Double_t dInvM = pV0->KineLambda().M();
      if ((dInvM<dUpperL) || (dInvM>=dLowerR)) { pV0 = 0; continue; }
      
      dVar[0] = 1.;
    }

    if (pV0->IsAntiLa()) {
      Int_t k = fHistoAntiLaInvM->FindBin(dPt); if (k<=0) { pV0 = 0; continue; }

      Double_t dMean  = fHistoAntiLaInvM->GetBinContent(k);
      Double_t dSigma = fHistoAntiLaInvM->GetBinError(k);

      Double_t dUpperL = dMean - (fLaCutNS * dSigma);
      Double_t dLowerR = dMean + (fLaCutNS * dSigma);

      Double_t dInvM = pV0->KineAntiLa().M();
      if ((dInvM<dUpperL) || (dInvM>=dLowerR)) { pV0 = 0; continue; }

      dVar[0] = 2.;
    }

    if (dVar[0]<-0.5) { pV0 = 0; continue; }

    dVar[1] = -1.;
    if (IsV0InJet(vV0,10.)) dVar[1] = 0.5;
    if (IsV0InJet(vV0,15.)) dVar[1] = 1.5;
    if (IsV0InJet(vV0,20.)) dVar[1] = 2.5;
    if (IsV0InJet(vV0,25.)) dVar[1] = 3.5;
    if (dVar[1]<0.) { pV0 = 0; continue; }

    dVar[2] = dV0M;
    dVar[3] = dV0A;
    dVar[4] = dCL1;
    dVar[5] = dZNA;

    dVar[6] = vV0.Eta();
    dVar[7] = dPt;

    hs->Fill(dVar);

    pV0 = 0;
  }

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetV0CF::FillKineInfo()
{
//
//  AliAnalysisTaskEmcalJetV0CF::FillKineInfo
//

  Double_t dV0M = fCentInfo->GetCentralityPercentile("V0M");
  Double_t dV0A = fCentInfo->GetCentralityPercentile("V0A");
  Double_t dCL1 = fCentInfo->GetCentralityPercentile("CL1");
  Double_t dZNA = fCentInfo->GetCentralityPercentile("ZNA");

  THnSparseD *hs = dynamic_cast<THnSparseD*>(fListUserOutputs->FindObject("hsKine"));
  if (!hs) {
    return kTRUE; // Should not happen; make Coverity happy
  }
//=============================================================================

  AliStack *pStack = 0;
  if (fEventESD) { pStack = MCEvent()->Stack(); if (!pStack) return kTRUE; }
//=============================================================================

  TParticle        *pESD = 0;
  AliAODMCParticle *pAOD = 0;
  for (Int_t i=0; i<MCEvent()->GetNumberOfTracks(); i++) {
    if (fEventAOD) { pAOD = (AliAODMCParticle*)MCEvent()->GetTrack(i);              if (!pAOD) continue; }
    if (fEventESD) { pESD =   ((AliMCParticle*)MCEvent()->GetTrack(i))->Particle(); if (!pESD) continue; }

    Bool_t bPhy = kFALSE;
    if (pAOD) bPhy =   pAOD->IsPhysicalPrimary();
    if (pESD) bPhy = pStack->IsPhysicalPrimary(i);
    if (!bPhy) { pAOD=0; pESD=0; continue; }

    Double_t  dEta = 0.;
    if (pAOD) dEta = pAOD->Eta();
    if (pESD) dEta = pESD->Eta();
    if ((dEta<fV0CutMinEta) || (dEta>=fV0CutMaxEta)) { pAOD=0; pESD=0; continue; }

    Int_t id = 0;
    if (pAOD) id = pAOD->GetPdgCode();
    if (pESD) id = pESD->GetPdgCode();

    Double_t dVar[8]; dVar[0] = -1.;
    if (id== 310 ) dVar[0] = 0.;
    if (id== 3122) dVar[0] = 1.;
    if (id==-3122) dVar[0] = 2.;
    if (dVar[0]<-0.5) { pAOD=0; pESD=0; continue; }

    TVector3 vV0;
    if (pAOD) vV0.SetXYZ(pAOD->Px(), pAOD->Py(), pAOD->Pz());
    if (pESD) vV0.SetXYZ(pESD->Px(), pESD->Py(), pESD->Pz());

    dVar[1] = -1.;
    if (IsV0InJet(vV0,10.)) dVar[1] = 0.5;
    if (IsV0InJet(vV0,15.)) dVar[1] = 1.5;
    if (IsV0InJet(vV0,20.)) dVar[1] = 2.5;
    if (IsV0InJet(vV0,25.)) dVar[1] = 3.5;
    if (dVar[1]<0.) { pAOD=0; pESD=0; continue; }

    dVar[2] = dV0M;
    dVar[3] = dV0A;
    dVar[4] = dCL1;
    dVar[5] = dZNA;
    dVar[6] = dEta;
    if (pAOD) dVar[7] = pAOD->Pt();
    if (pESD) dVar[7] = pESD->Pt();
    hs->Fill(dVar);

    pAOD = 0;
    pESD = 0;
  }

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
  Double_t dJetRadius = fJetsContRD->GetJetRadius();
  fJetsContRD->ResetCurrentID();
  AliEmcalJet *pJet = fJetsContRD->GetNextAcceptJet(); while (pJet) {
    Double_t dPt = fJetsContRD->GetJetPtCorr(fJetsContRD->GetCurrentID());
    if (dPt<dJetPtMin) { pJet = fJetsContRD->GetNextAcceptJet(); continue; }

    vJet.SetPtEtaPhi(dPt, pJet->Eta(), pJet->Phi());
    if (vJet.DeltaR(vV0)<dJetRadius) return kTRUE;
    pJet = fJetsContRD->GetNextAcceptJet();
  }

  return kFALSE;
}
