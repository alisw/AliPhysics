#include <TH2D.h>
#include <TVector2.h>
#include <THnSparse.h>
#include <TLorentzVector.h>

#include "AliLog.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliPicoTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskEmcalTmpSparseMaker.h"

ClassImp(AliAnalysisTaskEmcalTmpSparseMaker)

//_____________________________________________________________________________
AliAnalysisTaskEmcalTmpSparseMaker::AliAnalysisTaskEmcalTmpSparseMaker() :
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
fHnsTrks(0),
fHnsClus(0),
fHnsJets(0),
fListOutputEvH(0),
fListOutputTrk(0),
fListOutputClu(0),
fListOutputJet(0)
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::AliAnalysisTaskEmcalTmpSparseMaker
//
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalTmpSparseMaker::AliAnalysisTaskEmcalTmpSparseMaker(const char *name, const Bool_t bHistos) :
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
fHnsTrks(0),
fHnsClus(0),
fHnsJets(0),
fListOutputEvH(0),
fListOutputTrk(0),
fListOutputClu(0),
fListOutputJet(0)
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::AliAnalysisTaskEmcalTmpSparseMaker
//

  SetMakeGeneralHistograms(bHistos);

  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskEmcalTmpSparseMaker::~AliAnalysisTaskEmcalTmpSparseMaker()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::~AliAnalysisTaskEmcalTmpSparseMaker
//

  if (fContJets) { delete fContJets; fContJets = 0; }
  if (fContTrks) { delete fContTrks; fContTrks = 0; }
  if (fContClus) { delete fContClus; fContClus = 0; }

  if (fHnsEveH) { delete fHnsEveH; fHnsEveH = 0; }
  if (fHnsTrks) { delete fHnsTrks; fHnsTrks = 0; }
  if (fHnsClus) { delete fHnsClus; fHnsClus = 0; }
  if (fHnsJets) { delete fHnsJets; fHnsJets = 0; }

  if (fListOutputEvH) { delete fListOutputEvH; fListOutputEvH = 0; }
  if (fListOutputTrk) { delete fListOutputTrk; fListOutputTrk = 0; }
  if (fListOutputClu) { delete fListOutputClu; fListOutputClu = 0; }
  if (fListOutputJet) { delete fListOutputJet; fListOutputJet = 0; }
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::Init()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::Init
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::UserCreateOutputObjects
//

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
//=============================================================================

  fContJets = GetJetContainer(fNameJet.Data());

  fContTrks = GetParticleContainer(0); if (fContTrks) fContTrks->SetClassName("AliVTrack");
  fContClus = GetClusterContainer(0);  if (fContClus) fContClus->SetClassName("AliVCluster");
//=============================================================================

  if (fListOutputEvH) { delete fListOutputEvH; fListOutputEvH = 0; }
  if (fListOutputTrk) { delete fListOutputTrk; fListOutputTrk = 0; }
  if (fListOutputClu) { delete fListOutputClu; fListOutputClu = 0; }
  if (fListOutputJet) { delete fListOutputJet; fListOutputJet = 0; }

  fListOutputEvH = new TList();
  if (fContTrks) fListOutputTrk = new TList();
  if (fContClus) fListOutputClu = new TList();
  if (fContJets) fListOutputJet = new TList();
//=============================================================================

  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  MakeSparseEveH();
  if (fContTrks) MakeSparseTrks();
  if (fContClus) MakeSparseClus();
  if (fContJets) MakeSparseJets();

  TH1::AddDirectory(bStatusTmpH);
//=============================================================================

  PostData(2, fListOutputEvH);
  if (fListOutputTrk) PostData(3, fListOutputTrk);
  if (fListOutputClu) PostData(4, fListOutputClu);
  if (fListOutputJet) PostData(5, fListOutputJet);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::Terminate(Option_t *opt)
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::Terminate
//

  AliAnalysisTaskEmcalJet::Terminate(opt);

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalTmpSparseMaker::Run()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::Run
//

  if (fContTrks) FillSparseTrks();
  if (fContClus) FillSparseClus();
  if (fContJets) FillSparseJets();
//=============================================================================

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalTmpSparseMaker::FillHistograms()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::FillHistograms
//

  if (!AliAnalysisTaskEmcalJet::FillHistograms()) return kFALSE;
//=============================================================================

  FillSparseEveH();
//=============================================================================

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEmcalTmpSparseMaker::FillGeneralHistograms()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::FillGeneralHistograms
//

  if (!AliAnalysisTaskEmcalJet::FillGeneralHistograms()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::ExecOnce()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::ExecOnce
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
Bool_t AliAnalysisTaskEmcalTmpSparseMaker::RetrieveEventObjects()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::RetrieveEventObjects
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
Bool_t AliAnalysisTaskEmcalTmpSparseMaker::IsEventSelected()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::IsEventSelected
//

  if (!AliAnalysisTaskEmcalJet::IsEventSelected()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::FillSparseEveH()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::FillSparseEveH
//

  if (!fHnsEveH) return;
//=============================================================================

  Double_t dV[] = { fCent, fMtCh, fMtEm, fVertex[2], fEPV0, fAPhi, fRhoV };
  fHnsEveH->Fill(dV);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::FillSparseTrks()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::FillSparseTrks
//

  if ((!fContTrks) || (!fHnsTrks)) return;
//=============================================================================

  AliVTrack *pTrk  = 0;
  fContTrks->ResetCurrentID();
  while ((pTrk=static_cast<AliVTrack*>(fContTrks->GetNextAcceptParticle()))) {
    Double_t dV[] = { fCent, fMtCh, fVertex[2], fEPV0, fAPhi,
                      CalcRelPhiEP(pTrk->Phi()),
                      pTrk->Pt(), pTrk->Eta(), pTrk->Phi() };

    fHnsTrks->Fill(dV);
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::FillSparseClus()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::FillSparseClus
//

  if ((!fContClus) || (!fHnsClus)) return;
//=============================================================================

  AliVCluster *pClu = 0;
  fContClus->ResetCurrentID();
  while ((pClu=static_cast<AliVCluster*>(fContClus->GetNextAcceptCluster()))) {

    TLorentzVector vClu;
    pClu->GetMomentum(vClu, fVertex);
    Double_t dV[] = { fCent, fMtEm, fVertex[2], fEPV0, fAPhi,
                      CalcRelPhiEP(vClu.Phi()),
                      vClu.Pt(), vClu.Eta(), vClu.Phi() };

    fHnsClus->Fill(dV);
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::FillSparseJets()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::FillSparseJets
//

  if ((!fContJets) || (!fHnsJets)) return;
//=============================================================================

  Int_t    kLjetID[2] = { -1,  -1  };
  Double_t dLjetPt[2] = { -1., -1. };

  AliEmcalJet *pJet = 0;
  fContJets->ResetCurrentID();
  while ((pJet=fContJets->GetNextAcceptJet())) {
    Double_t dPtJ = fContJets->GetJetPtCorr(fContJets->GetCurrentID());

    if (dPtJ>dLjetPt[0]) {
      dLjetPt[1] = dLjetPt[0];
      kLjetID[1] = kLjetID[0];

      dLjetPt[0] = dPtJ;
      kLjetID[0] = fContJets->GetCurrentID();
    } else if (dPtJ>dLjetPt[1]) {
      dLjetPt[1] = dPtJ;
      kLjetID[1] = fContJets->GetCurrentID();
    }
  }
//=============================================================================

  for (Int_t j=0; j<2; j++) {
    if ((kLjetID[j]<0) || (dLjetPt[j]<0.)) continue;
    pJet = fContJets->GetAcceptJet(kLjetID[j]); if (!pJet) continue;
//=============================================================================

    AliParticleContainer *pContTrk = fContJets->GetParticleContainer();
    if (pContTrk) for (Int_t i=0; i<pJet->GetNumberOfTracks(); i++) {
      AliVParticle *pTrk = pJet->TrackAt(i, pContTrk->GetArray());

      if (pTrk) {
        Double_t dV[] = { fCent, fMtCh, fVertex[2], fEPV0, fAPhi,
                          CalcRelPhiEP(pTrk->Phi()), dLjetPt[j],
                          pTrk->Pt(), pTrk->Eta(), pTrk->Phi(), 2.*j };

        fHnsJets->Fill(dV);
      }
    }
//=============================================================================

    AliClusterContainer *pContClus = fContJets->GetClusterContainer();
    if (pContClus) for (Int_t i=0; i<pJet->GetNumberOfClusters(); i++) {
      AliVCluster *pClu = pJet->ClusterAt(i, pContClus->GetArray());

      if (pClu) {
        TLorentzVector vClu;
        pClu->GetMomentum(vClu, fVertex);
        Double_t dV[] = { fCent, fMtEm, fVertex[2], fEPV0, fAPhi,
                          CalcRelPhiEP(vClu.Phi()), dLjetPt[j],
                          vClu.Pt(), vClu.Eta(), vClu.Phi(), 2.*j+1. };

        fHnsJets->Fill(dV);
      }
    }
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::MakeSparseEveH()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::MakeSparseEveH
//

  if (!fListOutputEvH) return;
//=============================================================================

  const Int_t nV = 7;
  const Double_t dPhiNeg = -1.*TMath::Pi();
  const Double_t dPhiPos =  1.*TMath::Pi()  + TMath::Pi()/50.;
  const Double_t dPhiMax = TMath::TwoPi()   + TMath::Pi()/50.;
  const TString  sBin[] = { "aCent", "aMtCh", "aMtEm", "aVz", "aEPhi", "aAPhi", "aRho" };
  const Int_t    nBin[] = {     100,     600,     600,    20,     101,     101,    500 };
  const Double_t dMin[] = {      0.,      0.,      0.,  -10., dPhiNeg,      0.,     0. };
  const Double_t dMax[] = {    100.,   3000.,   3000.,   10., dPhiPos, dPhiMax,   500. };

  if (fHnsEveH) { delete fHnsEveH; fHnsEveH = 0; }
  fHnsEveH = new THnSparseD("hsEveH", "", nV, nBin, dMin, dMax);
  for (Int_t i=0; i<nV; i++) fHnsEveH->GetAxis(i)->SetName(sBin[i]);
  fListOutputEvH->Add(fHnsEveH);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::MakeSparseTrks()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::MakeSparseTrks
//

  if (!fListOutputTrk) return;
//=============================================================================

  const Int_t nV = 9;
  const Double_t dPhiNeg = -1.*TMath::Pi();
  const Double_t dPhiPos =  1.*TMath::Pi()  + TMath::Pi()/50.;
  const Double_t dPhiMax = TMath::TwoPi()   + TMath::Pi()/50.;
  const Double_t dPhiMin = TMath::PiOver2() + TMath::Pi()/50.;
  const TString  sBin[] = { "aCent", "aMult", "aVz", "aEPhi", "aAPhi", "aRPhi", "aPt", "aEta",  "aPhi" };
  const Int_t    nBin[] = {     100,     600,    20,     101,     101,      26,  1000,    200,     101 };
  const Double_t dMin[] = {      0.,      0.,  -10., dPhiNeg,      0.,      0.,    0.,    -1.,      0. };
  const Double_t dMax[] = {    100.,   3000.,   10., dPhiPos, dPhiMax, dPhiMin,  100.,     1., dPhiMax };

  if (fHnsTrks) { delete fHnsTrks; fHnsTrks = 0; }
  fHnsTrks = new THnSparseD("hsTrks", "", nV, nBin, dMin, dMax);
  for (Int_t i=0; i<nV; i++) fHnsTrks->GetAxis(i)->SetName(sBin[i]);
  fListOutputTrk->Add(fHnsTrks);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::MakeSparseClus()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::MakeSparseClus
//

  if (!fListOutputClu) return;
//=============================================================================

  const Int_t nV = 9;
  const Double_t dPhiNeg = -1.*TMath::Pi();
  const Double_t dPhiPos =  1.*TMath::Pi()  + TMath::Pi()/50.;
  const Double_t dPhiMax = TMath::TwoPi()   + TMath::Pi()/50.;
  const Double_t dPhiMin = TMath::PiOver2() + TMath::Pi()/50.;
  const TString  sBin[] = { "aCent", "aMult", "aVz",   "aEPhi", "aAPhi", "aRPhi", "aPt", "aEta",  "aPhi" };
  const Int_t    nBin[] = {     100,     600,    20,       101,     101,      26,  1000,    200,     101 };
  const Double_t dMin[] = {      0.,      0.,  -10.,   dPhiNeg,      0.,      0.,    0.,    -1.,      0. };
  const Double_t dMax[] = {    100.,   3000.,   10.,   dPhiPos, dPhiMax, dPhiMin,  100.,     1., dPhiMax };

  if (fHnsClus) { delete fHnsClus; fHnsClus = 0; }
  fHnsClus = new THnSparseD("hsClus", "", nV, nBin, dMin, dMax);
  for (Int_t i=0; i<nV; i++) fHnsClus->GetAxis(i)->SetName(sBin[i]);
  fListOutputClu->Add(fHnsClus);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskEmcalTmpSparseMaker::MakeSparseJets()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::MakeSparseJets
//

  if (!fListOutputJet) return;
//=============================================================================

  const Int_t nV = 11;
  const Double_t dPhiNeg = -1.*TMath::Pi();
  const Double_t dPhiPos =  1.*TMath::Pi()  + TMath::Pi()/50.;
  const Double_t dPhiMax = TMath::TwoPi()   + TMath::Pi()/50.;
  const Double_t dPhiMin = TMath::PiOver2() + TMath::Pi()/50.;
  const TString  sBin[] = { "aCent", "aMult", "aVz", "aEPhi", "aAPhi", "aRPhi", "aJet", "aPt", "aEta",  "aPhi", "aType" };
  const Int_t    nBin[] = {     100,     600,    20,     101,     101,      26,    500,  1000,    200,     101,       4 };
  const Double_t dMin[] = {      0.,      0.,  -10., dPhiNeg,      0.,      0.,     0.,    0.,    -1.,      0.,    -0.5 };
  const Double_t dMax[] = {    100.,   3000.,   10., dPhiPos, dPhiMax, dPhiMin,    250,  100.,     1., dPhiMax,     3.5 };

  if (fHnsJets) { delete fHnsJets; fHnsJets = 0; }
  fHnsJets = new THnSparseD("hsJets", "", nV, nBin, dMin, dMax);
  for (Int_t i=0; i<nV; i++) fHnsJets->GetAxis(i)->SetName(sBin[i]);
  fListOutputJet->Add(fHnsJets);
//=============================================================================

  return;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskEmcalTmpSparseMaker::CalcRelPhiEP(Double_t dPhi)
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::CalcRelPhiEP
//

  Double_t dPhiSA = TVector2::Phi_mpi_pi(dPhi);
  Double_t dPhiEP = TVector2::Phi_mpi_pi(dPhiSA - fEPV0);

  dPhiEP = TMath::Abs(dPhiEP);
  if (dPhiEP>TMath::PiOver2()) dPhiEP = TMath::Pi() - dPhiEP;

  return dPhiEP;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskEmcalTmpSparseMaker::CalcAysPlane()
{
//
//  AliAnalysisTaskEmcalTmpSparseMaker::CalcAysPlane
//

  if (!fContTrks) return -9999.;

  AliVParticle *pTrk = 0;
  Double_t dQx = 0, dQy = 0.;
  fContTrks->ResetCurrentID();
  while ((pTrk=fContTrks->GetNextAcceptParticle())) {
    dQx += TMath::Cos(pTrk->Phi());
    dQy += TMath::Sin(pTrk->Phi());
  } TVector2 vQ(dQx,dQy);

  return (vQ.Phi());
}
