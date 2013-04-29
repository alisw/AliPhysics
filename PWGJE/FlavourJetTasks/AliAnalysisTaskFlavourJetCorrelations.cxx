// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for Dmesons - jet correlations analysis
//
// Author: Xiaoming Zhang, xmzhang@lbl.gov
///////////////////////////////////////////////////////////////

#include <iostream>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>

#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliEmcalJet.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisTaskFlavourJetCorrelations.h"

ClassImp(AliAnalysisTaskFlavourJetCorrelations)

//_____________________________________________________________________________
AliAnalysisTaskFlavourJetCorrelations::AliAnalysisTaskFlavourJetCorrelations() :
AliAnalysisTaskEmcalJet(),
fUsedDzeros(0),
fUsedD0bars(0),
fUsedDstars(0),
fListControlHistos(0),
fListAnDzeroHistos(0),
fListAnDstarHistos(0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliAnalysisTaskFlavourJetCorrelations::AliAnalysisTaskFlavourJetCorrelations(const char *name, Bool_t bIsHisto) :
AliAnalysisTaskEmcalJet(name, bIsHisto),
fUsedDzeros(0),
fUsedD0bars(0),
fUsedDstars(0),
fListControlHistos(0),
fListAnDzeroHistos(0),
fListAnDstarHistos(0)
{
//
// Constructor
//

  SetMakeGeneralHistograms(bIsHisto);
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskFlavourJetCorrelations::~AliAnalysisTaskFlavourJetCorrelations()
{
//
// Default destructor
//

  if (fUsedDzeros) { delete fUsedDzeros; fUsedDzeros=0; }
  if (fUsedD0bars) { delete fUsedD0bars; fUsedD0bars=0; }
  if (fUsedDstars) { delete fUsedDstars; fUsedDstars=0; }

  if (fListControlHistos) { delete fListControlHistos; fListControlHistos=0; }
  if (fListAnDzeroHistos) { delete fListAnDzeroHistos; fListAnDzeroHistos=0; }
  if (fListAnDstarHistos) { delete fListAnDstarHistos; fListAnDstarHistos=0; }
}

//_____________________________________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::UserCreateOutputObjects()
{
//
// AliAnalysisTaskFlavourJetCorrelations::UserCreateOutputObjects
//


  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if (!fListControlHistos) fListControlHistos = new TList(); fListControlHistos->SetOwner();
  if (!fListAnDzeroHistos) fListAnDzeroHistos = new TList(); fListAnDzeroHistos->SetOwner();
  if (!fListAnDstarHistos) fListAnDstarHistos = new TList(); fListAnDstarHistos->SetOwner();

  MakeControlHistograms();
  CreateDzeroHistograms();
  CreateDstarHistograms();

  PostData(2, fListControlHistos);
  PostData(3, fListAnDzeroHistos);
  PostData(4, fListAnDstarHistos);
  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlavourJetCorrelations::Run()
{
//
// AliAnalysisTaskFlavourJetCorrelations::Run
//

  const Int_t    nJetPtBin   = 7;
  const Double_t dJetPtBin[] = { 0., 5., 10., 20., 30., 50., 80., 120. };
//=============================================================================

  AliEmcalJet *pJet = 0; 
  for (Int_t iJet=0; iJet<fJets->GetEntriesFast(); iJet++) {
    pJet = static_cast<AliEmcalJet*>(fJets->At(iJet)); if (!pJet) continue;
    if (!AcceptJet(pJet)) { pJet=0; continue; }

    Int_t iJetPtBin = TMath::BinarySearch(nJetPtBin,dJetPtBin,pJet->Pt());
    if (iJetPtBin<0 || iJetPtBin>=nJetPtBin) continue;

    if (fUsedDzeros) RunDzeroJet(pJet, iJetPtBin, kTRUE);
    if (fUsedD0bars) RunDzeroJet(pJet, iJetPtBin, kFALSE);
    if (fUsedDstars) RunDstarJet(pJet, iJetPtBin);
    pJet = 0;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::RunDzeroJet(AliEmcalJet const *pJet, const Int_t iJetPtBin, const Bool_t bIsD0)
{
//
// AliAnalysisTaskFlavourJetCorrelations::RunDzeroJet
//

  Int_t nCands = 0;
  if (bIsD0) nCands = fUsedDzeros->GetEntriesFast();
  else nCands = fUsedD0bars->GetEntriesFast(); if (nCands==0) return;
//=============================================================================

  const Int_t nms = kDzeroMatchType;
  const TString sMatched[nms] = { "MatchConeCandi", "MatchAreaCandi", "MatchConeProng", "MatchAreaProng" };
//=============================================================================

  Double_t dJetPt = pJet->Pt();
  TVector3 vJet, vRec, vec0, vec1;
  vJet.SetPtEtaPhi(dJetPt, pJet->Eta(), pJet->Phi());

  Double_t dMag1    = vJet.Mag();
  Double_t dMag2    = vJet.Mag2();
  Double_t dJetArea = pJet->Area();
//=============================================================================

  Bool_t bIsMatch[nms], bIsUsed[nms];
  for (Int_t i=0; i<nms; i++) { bIsMatch[i]=kFALSE; bIsUsed[i]=kFALSE; }

  AliAODRecoDecayHF2Prong *candi = 0;
  for (Int_t iCand=0; iCand<nCands; iCand++) {
    if (bIsD0)
      candi = dynamic_cast<AliAODRecoDecayHF2Prong*>(fUsedDzeros->At(iCand));
    else
      candi = dynamic_cast<AliAODRecoDecayHF2Prong*>(fUsedD0bars->At(iCand));
    if (!candi) continue;
    if (candi->Pt()>dJetPt) continue;
//=============================================================================

    Double_t dMassDzero = candi->InvMassD0();
    Double_t dMassD0bar = candi->InvMassD0bar();
    vRec.SetXYZ(candi->Px(), candi->Pt(), candi->Pz());
    vec0.SetXYZ(candi->PxProng(0), candi->PyProng(0), candi->PzProng(0));
    vec1.SetXYZ(candi->PxProng(1), candi->PyProng(1), candi->PzProng(1));

    Double_t delRecoR = vRec.DeltaR(vJet);
    Double_t delTrk0R = vec0.DeltaR(vJet);
    Double_t delTrk1R = vec1.DeltaR(vJet);

    bIsMatch[kMatchConeCandi] =  (delRecoR<fJetRadius);
    bIsMatch[kMatchAreaCandi] =  (delRecoR<dJetArea);
    bIsMatch[kMatchConeProng] = ((delTrk0R<fJetRadius) && (delTrk1R<fJetRadius));
    bIsMatch[kMatchAreaProng] = ((delTrk0R<dJetArea)   && (delTrk1R<dJetArea));

    Bool_t bNoMatch = kTRUE;
    for (Int_t i=0; i<nms; i++) if (bIsMatch[i]) { bIsUsed[i]=kTRUE; bNoMatch=kFALSE; } if (bNoMatch) continue;
//=============================================================================

    Double_t dRecoFFs = vRec.Dot(vJet) / dMag2;
    Double_t dRecoRel = vRec.Cross(vJet).Mag() / dMag1;
    Double_t dRecoRho = 0.5 / TMath::Pi() / ((delRecoR<1e-12) ? 1e-12 : delRecoR);

    Double_t dTrk0FFs = vec0.Dot(vJet) / dMag2;
    Double_t dTrk0Rel = vec0.Cross(vJet).Mag() / dMag1;
    Double_t dTrk0Rho = 0.5 / TMath::Pi() / ((delTrk0R<1e-12) ? 1e-12 : delTrk0R);

    Double_t dTrk1FFs = vec1.Dot(vJet) / dMag2;
    Double_t dTrk1Rel = vec1.Cross(vJet).Mag() / dMag1;
    Double_t dTrk1Rho = 0.5 / TMath::Pi() / ((delTrk1R<1e-12) ? 1e-12 : delTrk1R);

    for (Int_t i=0; i<nms; i++) if (bIsMatch[i]) {
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_InvM_JetPt_%s",sMatched[i].Data())))->Fill(dMassDzero,dJetPt);
  
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Reco_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dRecoFFs,dMassDzero);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Reco_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dRecoRel,dMassDzero);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Reco_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delRecoR,dMassDzero,dRecoRho);

      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Trks_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk0FFs,dMassDzero);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Trks_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk0Rel,dMassDzero);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Trks_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delTrk0R,dMassDzero,dTrk0Rho);
  
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Trks_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk1FFs,dMassDzero);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Trks_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk1Rel,dMassDzero);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_Trks_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delTrk1R,dMassDzero,dTrk1Rho);

      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hDzero_2ProngsRho_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delTrk0R,delTrk1R);
//=============================================================================

      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_InvM_JetPt_%s",sMatched[i].Data())))->Fill(dMassD0bar,dJetPt);

      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Reco_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dRecoFFs,dMassD0bar);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Reco_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dRecoRel,dMassD0bar);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Reco_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delRecoR,dMassD0bar,dRecoRho);

      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Trks_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk0FFs,dMassD0bar);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Trks_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk0Rel,dMassD0bar);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Trks_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delTrk0R,dMassD0bar,dTrk0Rho);

      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Trks_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk1FFs,dMassD0bar);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Trks_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk1Rel,dMassD0bar);
      ((TH2D*)fListAnDzeroHistos->FindObject(Form("hD0bar_Trks_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delTrk1R,dMassD0bar,dTrk1Rho);
    }

    candi = 0;
  }

  for (Int_t i=0; i<nms; i++) if (bIsUsed[i]) ((TH1D*)fListAnDzeroHistos->FindObject(Form("hDzero_UsedJetPt_%s",sMatched[i].Data())))->Fill(dJetPt);

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::RunDstarJet(AliEmcalJet const *pJet, const Int_t iJetPtBin)
{
//
// AliAnalysisTaskFlavourJetCorrelations::RunDstarJet
//

  Int_t nCands = fUsedDstars->GetEntriesFast(); if (nCands==0) return;
//=============================================================================

  const Int_t nms = kDzeroMatchType;
  const TString sMatched[nms] = { "MatchConeCandi", "MatchAreaCandi", "MatchConeProng", "MatchAreaProng" };
//=============================================================================

  Double_t dJetPt = pJet->Pt();
  TVector3 vJet, vRec, vec0, vec1;
  vJet.SetPtEtaPhi(dJetPt, pJet->Eta(), pJet->Phi());

  Double_t dMag1    = vJet.Mag();
  Double_t dMag2    = vJet.Mag2();
  Double_t dJetArea = pJet->Area();
//=============================================================================

  Bool_t bIsMatch[nms], bIsUsed[nms];
  for (Int_t i=0; i<nms; i++) { bIsMatch[i]=kFALSE; bIsUsed[i]=kFALSE; }

//AliAODTrack             *bache = 0
  AliAODRecoCascadeHF     *candi = 0;
  AliAODRecoDecayHF2Prong *dZero = 0;
  for (Int_t iCand=0; iCand<nCands; iCand++) {
    candi = dynamic_cast<AliAODRecoCascadeHF*>(fUsedDstars->At(iCand)); if (!candi) continue;
    dZero = candi->Get2Prong();   if (!dZero) continue;
//  bache = candi->GetBachelor(); if (!bache) continue;
    if (candi->Pt()>dJetPt) continue;
//=============================================================================

    Double_t dInvM = candi->DeltaInvMass();
    vRec.SetXYZ(candi->Px(), candi->Pt(), candi->Pz());
    vec0.SetXYZ(dZero->PxProng(0), dZero->PyProng(0), dZero->PzProng(0));
    vec1.SetXYZ(dZero->PxProng(1), dZero->PyProng(1), dZero->PzProng(1));
 
    Double_t delRecoR = vRec.DeltaR(vJet);
    Double_t delTrk0R = vec0.DeltaR(vJet);
    Double_t delTrk1R = vec1.DeltaR(vJet);

    bIsMatch[kMatchConeCandi] =  (delRecoR<fJetRadius);
    bIsMatch[kMatchAreaCandi] =  (delRecoR<dJetArea);
    bIsMatch[kMatchConeProng] = ((delTrk0R<fJetRadius) && (delTrk1R<fJetRadius));
    bIsMatch[kMatchAreaProng] = ((delTrk0R<dJetArea)   && (delTrk1R<dJetArea));

    Bool_t bNoMatch = kTRUE;
    for (Int_t i=0; i<nms; i++) if (bIsMatch[i]) { bIsUsed[i]=kTRUE; bNoMatch=kFALSE; } if (bNoMatch) continue;
//=============================================================================

    Double_t dRecoFFs = vRec.Dot(vJet) / dMag2;
    Double_t dRecoRel = vRec.Cross(vJet).Mag() / dMag1;
    Double_t dRecoRho = 0.5 / TMath::Pi() / ((delRecoR<1e-12) ? 1e-12 : delRecoR);
    
    Double_t dTrk0FFs = vec0.Dot(vJet) / dMag2;
    Double_t dTrk0Rel = vec0.Cross(vJet).Mag() / dMag1;
    Double_t dTrk0Rho = 0.5 / TMath::Pi() / ((delTrk0R<1e-12) ? 1e-12 : delTrk0R);
    
    Double_t dTrk1FFs = vec1.Dot(vJet) / dMag2;
    Double_t dTrk1Rel = vec1.Cross(vJet).Mag() / dMag1;
    Double_t dTrk1Rho = 0.5 / TMath::Pi() / ((delTrk1R<1e-12) ? 1e-12 : delTrk1R);
    
    for (Int_t i=0; i<nms; i++) if (bIsMatch[i]) {
      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_InvM_JetPt_%s",sMatched[i].Data())))->Fill(dInvM,dJetPt);

      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Reco_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dRecoFFs,dInvM);
      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Reco_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dRecoRel,dInvM);
      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Reco_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delRecoR,dInvM,dRecoRho);

      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Trks_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk0FFs,dInvM);
      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Trks_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk0Rel,dInvM);
      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Trks_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delTrk0R,dInvM,dTrk0Rho);

      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Trks_FFs_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk1FFs,dInvM);
      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Trks_Rel_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(dTrk1Rel,dInvM);
      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_Trks_Rho_InvM_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delTrk1R,dInvM,dTrk1Rho);

      ((TH2D*)fListAnDstarHistos->FindObject(Form("hDstar_2ProngsRho_%s_JetPtBin%d",sMatched[i].Data(),iJetPtBin)))->Fill(delTrk0R,delTrk1R);
    }

    dZero = 0;
    candi = 0;
  }

  for (Int_t i=0; i<nms; i++) if (bIsUsed[i]) ((TH1D*)fListAnDstarHistos->FindObject(Form("hDstar_UsedJetPt_%s",sMatched[i].Data())))->Fill(dJetPt);

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlavourJetCorrelations::RetrieveEventObjects()
{
//
// AliAnalysisTaskFlavourJetCorrelations::RetrieveEventObjects
//

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects()) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlavourJetCorrelations::IsEventSelected()
{
//
// AliAnalysisTaskFlavourJetCorrelations::IsEventSelected
//

  if (!AliAnalysisTaskEmcal::IsEventSelected()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlavourJetCorrelations::FillHistograms()
{
//
// AliAnalysisTaskFlavourJetCorrelations::FillHistograms
//

  if (!AliAnalysisTaskEmcal::FillHistograms()) return kFALSE;

  if (!FillControlHistograms()) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlavourJetCorrelations::FillGeneralHistograms()
{
//
// AliAnalysisTaskFlavourJetCorrelations::FillGeneralHistograms
//

  if (!AliAnalysisTaskEmcal::FillGeneralHistograms()) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::ExecOnce()
{
//
// AliAnalysisTaskFlavourJetCorrelations::ExecOnce
//

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (!fUsedDzeros) fUsedDzeros = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AnaUsedDzero"));
  if (!fUsedD0bars) fUsedD0bars = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AnaUsedD0bar"));
  if (!fUsedDstars) fUsedDstars = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AnaUsedDstar"));

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::CreateDzeroHistograms()
{
//
// AliAnalysisTaskFlavourJetCorrelations::CreateDzeroHistograms
//

  if (!fListAnDzeroHistos) return;
  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nms = kDzeroMatchType;
  const TString sMatched[nms] = { "MatchConeCandi", "MatchAreaCandi", "MatchConeProng", "MatchAreaProng" };
  const Double_t dMass = TDatabasePDG::Instance()->GetParticle(421)->Mass();

  const Int_t    nMassBin = 200;
  const Double_t dMassMin = dMass - 0.15;
  const Double_t dMassMax = dMass + 0.15;

  const Int_t    nJetPtBin   = 7;
  const Double_t aJetPtBin[] = { 0., 5., 10., 20., 30., 50., 80., 120. };

  const Int_t    nFFsBin   =   21;
  const Double_t dFFsBin[] = { 0.0030, 0.0040, 0.0052, 0.0069, 0.0092, 0.0121, 0.0160,
                               0.0212, 0.0280, 0.0371, 0.0490, 0.0648, 0.0856, 0.1132,
                               0.1497, 0.1980, 0.2618, 0.3461, 0.4576, 0.6050, 0.8000, 1. };

  const Int_t    nRelBin   =   12;
  const Double_t dRelBin[] = { 0., 0.06, 0.12, 0.2084, 0.2802, 0.3769, 0.5069, 0.6817, 0.9169, 1.2331, 1.6585, 2.2306, 3. };

  const Int_t    nRhoBin = 20;
  const Double_t dRhoMin = 0.;
  const Double_t dRhoMax = 1.;
//=============================================================================

  TH1D *h1 = 0;
  TH2D *h2 = 0;
  for (Int_t i=0; i<nms; i++) {
    h1 = new TH1D(Form("hDzero_UsedJetPt_%s",sMatched[i].Data()), "", nJetPtBin, aJetPtBin);
    h1->Sumw2(); fListAnDzeroHistos->Add(h1); h1=0;

    h2 = new TH2D(Form("hDzero_InvM_JetPt_%s",sMatched[i].Data()), "", nMassBin, dMassMin, dMassMax, nJetPtBin, aJetPtBin);
    h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

    h2 = new TH2D(Form("hD0bar_InvM_JetPt_%s",sMatched[i].Data()), "", nMassBin, dMassMin, dMassMax, nJetPtBin, aJetPtBin);
    h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

    TString sName;
    for (Int_t j=0; j<nJetPtBin; j++) {
      sName = "hDzero_Reco_FFs_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nFFsBin, dFFsBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hDzero_Reco_Rel_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRelBin, dRelBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hDzero_Reco_Rho_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRhoBin, dRhoMin, dRhoMax, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hDzero_Trks_FFs_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nFFsBin, dFFsBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hDzero_Trks_Rel_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRelBin, dRelBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hDzero_Trks_Rho_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRhoBin, dRhoMin, dRhoMax, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hDzero_2ProngsRho";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRhoBin, dRhoMin, dRhoMax, nRhoBin, dRhoMin, dRhoMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;
//=============================================================================

      sName = "hD0bar_Reco_FFs_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nFFsBin, dFFsBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hD0bar_Reco_Rel_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRelBin, dRelBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hD0bar_Reco_Rho_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRhoBin, dRhoMin, dRhoMax, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hD0bar_Trks_FFs_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nFFsBin, dFFsBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hD0bar_Trks_Rel_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRelBin, dRelBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;

      sName = "hD0bar_Trks_Rho_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRhoBin, dRhoMin, dRhoMax, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDzeroHistos->Add(h2); h2=0;
    }
  }

  TH1::AddDirectory(bStatusTmpH);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::CreateDstarHistograms()
{
//
// AliAnalysisTaskFlavourJetCorrelations::CreateDstarHistograms
//

  if (!fListAnDstarHistos) return;
  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nms = kDzeroMatchType;
  const TString sMatched[nms] = { "MatchConeCandi", "MatchAreaCandi", "MatchConeProng", "MatchAreaProng" };
  Double_t dMassDzero = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t dMassDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t dMassDelta = dMassDstar - dMassDzero;

  const Int_t    nMassBin = 200;
  const Double_t dMassMin = dMassDelta - 0.015;
  const Double_t dMassMax = dMassDelta + 0.015;

  const Int_t    nJetPtBin   = 7;
  const Double_t aJetPtBin[] = { 0., 5., 10., 20., 30., 50., 80., 120. };

  const Int_t    nFFsBin   =   21;
  const Double_t dFFsBin[] = { 0.0030, 0.0040, 0.0052, 0.0069, 0.0092, 0.0121, 0.0160,
                               0.0212, 0.0280, 0.0371, 0.0490, 0.0648, 0.0856, 0.1132,
                               0.1497, 0.1980, 0.2618, 0.3461, 0.4576, 0.6050, 0.8000, 1. };

  const Int_t    nRelBin   =   12;
  const Double_t dRelBin[] = { 0., 0.06, 0.12, 0.2084, 0.2802, 0.3769, 0.5069, 0.6817, 0.9169, 1.2331, 1.6585, 2.2306, 3. };

  const Int_t    nRhoBin = 20;
  const Double_t dRhoMin = 0.;
  const Double_t dRhoMax = 1.;
//=============================================================================

  TH1D *h1 = 0;
  TH2D *h2 = 0;
  for (Int_t i=0; i<nms; i++) {
    h1 = new TH1D(Form("hDstar_UsedJetPt_%s",sMatched[i].Data()), "", nJetPtBin, aJetPtBin);
    h1->Sumw2(); fListAnDstarHistos->Add(h1); h1=0;

    h2 = new TH2D(Form("hDstar_InvM_JetPt_%s",sMatched[i].Data()), "", nMassBin, dMassMin, dMassMax, nJetPtBin, aJetPtBin);
    h2->Sumw2(); fListAnDstarHistos->Add(h2); h2=0;

    TString sName;
    for (Int_t j=0; j<nJetPtBin; j++) {
      sName = "hDstar_Reco_FFs_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nFFsBin, dFFsBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDstarHistos->Add(h2); h2=0;

      sName = "hDstar_Reco_Rel_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRelBin, dRelBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDstarHistos->Add(h2); h2=0;

      sName = "hDstar_Reco_Rho_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRhoBin, dRhoMin, dRhoMax, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDstarHistos->Add(h2); h2=0;

      sName = "hDstar_Trks_FFs_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nFFsBin, dFFsBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDstarHistos->Add(h2); h2=0;

      sName = "hDstar_Trks_Rel_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRelBin, dRelBin, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDstarHistos->Add(h2); h2=0;

      sName = "hDstar_Trks_Rho_InvM";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRhoBin, dRhoMin, dRhoMax, nMassBin, dMassMin, dMassMax);
      h2->Sumw2(); fListAnDstarHistos->Add(h2); h2=0;

      sName = "hDstar_2ProngsRho";
      h2 = new TH2D(Form("%s_%s_JetPtBin%d",sName.Data(),sMatched[i].Data(),j), "", nRhoBin, dRhoMin, dRhoMax, nRhoBin, dRhoMin, dRhoMax);
      h2->Sumw2(); fListAnDstarHistos->Add(h2); h2=0;
    }
  }

  TH1::AddDirectory(bStatusTmpH);
  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlavourJetCorrelations::FillControlHistograms()
{
//
// AliAnalysisTaskFlavourJetCorrelations::FillControlHistograms
//

  if (!fListControlHistos) return kFALSE;

  if (fTracks) {
    AliVTrack *pTrk = 0;
    for (Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++) {
      pTrk = static_cast<AliVTrack*>(fTracks->At(itrk)); if (!pTrk) continue;
      if (!AcceptTrack(pTrk)) { pTrk = 0; continue; }

      ((TH1D*)fListControlHistos->FindObject("hTrksPt"))->Fill(pTrk->Pt()); pTrk=0; 
    }
  }

  if (fCaloClusters) {
    AliVCluster *pCls = 0;
    TLorentzVector vLorentz;
    for (Int_t iClus=0; iClus<fCaloClusters->GetEntriesFast(); iClus++) {
      pCls = static_cast<AliVCluster*>(fCaloClusters->At(iClus)); if (!pCls) continue;

      pCls->GetMomentum(vLorentz, fVertex);
      ((TH1D*)fListControlHistos->FindObject("hClusPt"))->Fill(vLorentz.Pt()); pCls=0;
    }
  }

  if (fJets) {
    AliEmcalJet* pJet = 0;
    static Int_t sortedJets[9999] = { -1 };
    if (GetSortedArray(sortedJets,fJets)) if (sortedJets[0]>=0) {
      pJet = static_cast<AliEmcalJet*>(fJets->At(sortedJets[0])); 
      if (pJet) ((TH1D*)fListControlHistos->FindObject("hLeadingJets_Pt"))->Fill(pJet->Pt());
    } pJet = 0;

    for (Int_t iJet=0; iJet<fJets->GetEntriesFast(); iJet++) {
      pJet = static_cast<AliEmcalJet*>(fJets->At(iJet)); if (!pJet) continue;
      if (!AcceptJet(pJet)) continue;

      Double_t dPt   = pJet->Pt();
      Double_t dArea = pJet->Area();
      ((TH2D*)fListControlHistos->FindObject("hJets_Pt_Area"))->Fill(dPt,dArea);
      ((TH2D*)fListControlHistos->FindObject("hJets_Eta_Phi"))->Fill(pJet->Eta(),pJet->Phi());
      ((TH2D*)fListControlHistos->FindObject("hJets_LeadingPt"))->Fill(dPt,GetLeadingHadronPt(pJet));

      if (fRho) {
        ((TH2D*)fListControlHistos->FindObject("hJets_CorrPt_Area"))->Fill(dPt-fRhoVal*dArea,dArea);
      } pJet = 0;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::MakeControlHistograms()
{
//
// AliAnalysisTaskFlavourJetCorrelations::MakeControlHistograms
//

  if (!fListControlHistos) return;
  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1D *h1 = 0;
  TH2D *h2 = 0;

  if (!fTracksName.IsNull()) {
    h1 = new TH1D("hTrksPt", "", fNbins/2, fMinBinPt, fMaxBinPt/2.);
    h1->Sumw2(); fListControlHistos->Add(h1); h1=0;
  }

  if (!fCaloName.IsNull()) {
    h1 = new TH1D("hClusPt", "", fNbins/2, fMinBinPt, fMaxBinPt/2.);
    h1->Sumw2(); fListControlHistos->Add(h1); h1=0;
  }

  if (!fJetsName.IsNull()) {
    h1 = new TH1D("hLeadingJets_Pt", "", fNbins, fMinBinPt, fMaxBinPt);
    h1->Sumw2(); fListControlHistos->Add(h1); h1=0;

    h2 = new TH2D("hJets_Pt_Area", "", fNbins, fMinBinPt, fMaxBinPt, 30, 0., 3.);
    h2->Sumw2();  fListControlHistos->Add(h2); h2=0;

    h2 = new TH2D("hJets_Eta_Phi", "", 50, -1., 1., 101, 0., 2.*TMath::Pi() + TMath::Pi()/200.);
    h2->Sumw2(); fListControlHistos->Add(h2); h2=0;

    h2 = new TH2D("hJets_LeadingPt", "", fNbins, fMinBinPt, fMaxBinPt, fNbins/2, fMinBinPt, fMaxBinPt/2.);
    h2->Sumw2(); fListControlHistos->Add(h2); h2=0;

    if (!fRhoName.IsNull()) {
      h2 = new TH2D("hJets_CorrPt_Area", "", 2*fNbins, -1.*fMaxBinPt, fMaxBinPt, 30, 0., 3.);
    }
  }

  TH1::AddDirectory(bStatusTmpH);
  return;
}
