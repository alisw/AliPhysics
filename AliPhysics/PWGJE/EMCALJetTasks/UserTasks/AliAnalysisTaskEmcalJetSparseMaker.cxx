#include <TH2D.h>
#include <TVector2.h>
#include <THnSparse.h>
#include <TLorentzVector.h>

#include "AliLog.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVParticle.h"
#include "AliPicoTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskEmcalJetSparseMaker.h"

ClassImp(AliAnalysisTaskEmcalJetSparseMaker)

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetSparseMaker::AliAnalysisTaskEmcalJetSparseMaker() :
AliAnalysisTaskEmcalJet(),
fMtCh(0.),
fMtEm(0.),
fRhoV(0.),
fAPhi(0.),
fNameJet(""),
fContJets(0),
fContTrks(0),
fContClus(0),
fHnsEveH(0),
fHnsJets(0),
fListOutputEvH(0),
fListOutputJet(0)
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::AliAnalysisTaskEmcalJetSparseMaker
//
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetSparseMaker::AliAnalysisTaskEmcalJetSparseMaker(const char *name, const Bool_t bHistos) :
AliAnalysisTaskEmcalJet(name, bHistos),
fMtCh(0.),
fMtEm(0.),
fRhoV(0.),
fAPhi(0.),
fNameJet("cJets"),
fContJets(0),
fContTrks(0),
fContClus(0),
fHnsEveH(0),
fHnsJets(0),
fListOutputEvH(0),
fListOutputJet(0)
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::AliAnalysisTaskEmcalJetSparseMaker
//

  SetMakeGeneralHistograms(bHistos);

  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetSparseMaker::~AliAnalysisTaskEmcalJetSparseMaker()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::~AliAnalysisTaskEmcalJetSparseMaker
//

  if (fContJets) { delete fContJets; fContJets = 0; }
  if (fContTrks) { delete fContTrks; fContTrks = 0; }
  if (fContClus) { delete fContClus; fContClus = 0; }

  if (fHnsEveH) { delete fHnsEveH; fHnsEveH = 0; }
  if (fHnsJets) { delete fHnsJets; fHnsJets = 0; }

  if (fListOutputEvH) { delete fListOutputEvH; fListOutputEvH = 0; }
  if (fListOutputJet) { delete fListOutputJet; fListOutputJet = 0; }
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetSparseMaker::Init()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::Init
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetSparseMaker::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::UserCreateOutputObjects
//

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
//=============================================================================

  fContJets = GetJetContainer(fNameJet.Data());

  fContTrks = GetParticleContainer(0); if (fContTrks) fContTrks->SetClassName("AliVTrack");
  fContClus = GetClusterContainer(0);  if (fContClus) fContClus->SetClassName("AliVCluster");
//=============================================================================

  if (fListOutputEvH) { delete fListOutputEvH; fListOutputEvH = 0; } fListOutputEvH = new TList();
  if (fListOutputJet) { delete fListOutputJet; fListOutputJet = 0; } fListOutputJet = new TList();

  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  MakeSparse();

  TH1::AddDirectory(bStatusTmpH);
//=============================================================================

  PostData(2, fListOutputEvH);
  PostData(3, fListOutputJet);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetSparseMaker::Terminate(Option_t *opt)
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::Terminate
//

  AliAnalysisTaskEmcalJet::Terminate(opt);

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSparseMaker::Run()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::Run
//

  FillSparse();
//=============================================================================

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSparseMaker::FillHistograms()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::FillHistograms
//

  if (!AliAnalysisTaskEmcalJet::FillHistograms()) return kFALSE;
//=============================================================================

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSparseMaker::FillGeneralHistograms()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::FillGeneralHistograms
//

  if (!AliAnalysisTaskEmcalJet::FillGeneralHistograms()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetSparseMaker::ExecOnce()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::ExecOnce
//

  AliAnalysisTaskEmcalJet::ExecOnce();
  if (!fLocalInitialized) return;
//=============================================================================

  if (fContJets) if (fContJets->GetArray()==0) fContJets = 0x0;

  if (fContTrks) if (fContTrks->GetArray()==0) fContTrks = 0x0;
  if (fContClus) if (fContClus->GetArray()==0) fContClus = 0x0;

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSparseMaker::RetrieveEventObjects()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::RetrieveEventObjects
//

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects()) return kFALSE;
//=============================================================================

  if (fContTrks) fAPhi = CalcAysPlane();
  if (fContJets) fRhoV = fContJets->GetRhoVal();
  if (fContTrks) fMtCh = fContTrks->GetNAcceptedParticles();
  if (fContClus) fMtEm = fContClus->GetNAcceptedClusters();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSparseMaker::IsEventSelected()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::IsEventSelected
//

  if (!AliAnalysisTaskEmcalJet::IsEventSelected()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetSparseMaker::FillSparse()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::FillSparseJets
//

  if ((!fContJets) || (!fHnsEveH) || (!fHnsJets)) return;
//=============================================================================

  Bool_t bFillEH = kFALSE;
  const Double_t dR = fContJets->GetJetRadius();
  const Double_t dAreaRef = dR * dR * TMath::Pi();

  AliEmcalJet *pJet = 0;
  fContJets->ResetCurrentID();
  while ((pJet=fContJets->GetNextAcceptJet())) {
    Double_t dJet[] = { fCent, fVertex[2], fEPV0, fAPhi,
                        CalcRelPhiEP(pJet->Phi()),
                        fContJets->GetLeadingHadronPt(pJet),
                        pJet->Pt(), fContJets->GetJetPtCorr(fContJets->GetCurrentID()),
                        pJet->Eta(), pJet->Phi(), pJet->Area()/dAreaRef };

    fHnsJets->Fill(dJet);
    if (!bFillEH) bFillEH = kTRUE;
  }
//=============================================================================

  if (bFillEH) {
    Double_t dEve[] = { fCent, fMtCh, fMtEm, fVertex[2], fEPV0, fAPhi, fRhoV };
    fHnsEveH->Fill(dEve);
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetSparseMaker::MakeSparse()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::MakeSparseJets
//

  if (!fListOutputEvH) return;
  if (!fListOutputJet) return;
//=============================================================================

  const Int_t nEveV = 7;
  const Double_t dPhiNeg = -1.*TMath::Pi();
  const Double_t dPhiPos =  1.*TMath::Pi()  + TMath::Pi()/50.;
  const Double_t dPhiMax = TMath::TwoPi()   + TMath::Pi()/50.;
  const Double_t dPhiMin = TMath::PiOver2() + TMath::Pi()/50.;
  const TString  sEveBin[] = { "aCent", "aMtCh", "aMtEm", "aVz",   "aEPhi", "aAPhi", "aRho" };
  const Int_t    nEveBin[] = {     100,     600,     600,    20,     101,       101,    500 };
  const Double_t dEveMin[] = {      0.,      0.,      0.,  -10., dPhiNeg,        0.,     0. };
  const Double_t dEveMax[] = {    100.,   3000.,   3000.,   10., dPhiPos,   dPhiMax,   500. };

  if (fHnsEveH) { delete fHnsEveH; fHnsEveH = 0; }
  fHnsEveH = new THnSparseD("hsEveH", "", nEveV, nEveBin, dEveMin, dEveMax);
  for (Int_t i=0; i<nEveV; i++) fHnsEveH->GetAxis(i)->SetName(sEveBin[i]);
  fListOutputEvH->Add(fHnsEveH);
//=============================================================================

  const Int_t nJetV = 11;
  const TString  sJetBin[] = { "aCent", "aVz",   "aEPhi", "aAPhi", "aRPhi", "aLeading", "aPtb", "aPta", "aEta",  "aPhi", "aArea" };
  const Int_t    nJetBin[] = {     100,    20,     101,       101,      26,        500,    500,   1000,    200,     101,     100 };
  const Double_t dJetMin[] = {      0.,  -10., dPhiNeg,        0.,      0.,         0.,     0.,  -500.,    -1.,      0.,      0. };
  const Double_t dJetMax[] = {    100.,   10., dPhiPos,   dPhiMax, dPhiMin,       250.,   500.,   500.,     1., dPhiMax,      5. };

  if (fHnsJets) { delete fHnsJets; fHnsJets = 0; }
  fHnsJets = new THnSparseD("hsJets", "", nJetV, nJetBin, dJetMin, dJetMax);
  for (Int_t i=0; i<nJetV; i++) fHnsJets->GetAxis(i)->SetName(sJetBin[i]);
  fListOutputJet->Add(fHnsJets);
//=============================================================================

  return;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetSparseMaker::CalcRelPhiEP(Double_t dPhi)
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::CalcRelPhiEP
//

  Double_t dPhiSA = TVector2::Phi_mpi_pi(dPhi);
  Double_t dPhiEP = TVector2::Phi_mpi_pi(dPhiSA - fEPV0);

  dPhiEP = TMath::Abs(dPhiEP);
  if (dPhiEP>TMath::PiOver2()) dPhiEP = TMath::Pi() - dPhiEP;

  return dPhiEP;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetSparseMaker::CalcAysPlane()
{
//
//  AliAnalysisTaskEmcalJetSparseMaker::CalcAysPlane
//

  if (!fContTrks) return -9999.;
//=============================================================================

  AliVParticle *pTrk = 0;
  Double_t dQx = 0, dQy = 0.;
  fContTrks->ResetCurrentID();
  while ((pTrk=fContTrks->GetNextAcceptParticle())) {
    dQx += TMath::Cos(pTrk->Phi());
    dQy += TMath::Sin(pTrk->Phi());
  } TVector2 vQ(dQx,dQy);

  return (vQ.Phi());
}
