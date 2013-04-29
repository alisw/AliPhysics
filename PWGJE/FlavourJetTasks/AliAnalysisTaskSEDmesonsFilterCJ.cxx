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

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>

#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliAnalysisTaskSEDmesonsFilterCJ.h"

ClassImp(AliAnalysisTaskSEDmesonsFilterCJ)

//_____________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::AliAnalysisTaskSEDmesonsFilterCJ() :
AliAnalysisTaskSE(),
fEventAOD(0),
fCutDzero(0),
fCutDstar(0),
fDzeroClArr(0),
fDstarClArr(0),
fUsedDzeros(0),
fUsedD0bars(0),
fUsedDstars(0),
fIsNotExecOnce(kTRUE),
fListCutsDmesons(0),
fListDzeroHistos(0),
fListDstarHistos(0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::AliAnalysisTaskSEDmesonsFilterCJ(const char *name) :
AliAnalysisTaskSE(name),
fEventAOD(0),
fCutDzero(0),
fCutDstar(0),
fDzeroClArr(0),
fDstarClArr(0),
fUsedDzeros(0),
fUsedD0bars(0),
fUsedDstars(0),
fIsNotExecOnce(kTRUE),
fListCutsDmesons(0),
fListDzeroHistos(0),
fListDstarHistos(0)
{
//
// Constructor
//

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::~AliAnalysisTaskSEDmesonsFilterCJ()
{
//
// Default destructor
//

  if (fEventAOD)   { delete fEventAOD;   fEventAOD  =0; }

  if (fCutDzero)   { delete fCutDzero;   fCutDzero  =0; }
  if (fCutDstar)   { delete fCutDstar;   fCutDstar  =0; }

  if (fDzeroClArr) { delete fDzeroClArr; fDzeroClArr=0; }
  if (fDstarClArr) { delete fDstarClArr; fDstarClArr=0; }

  if (fUsedDzeros) { delete fUsedDzeros; fUsedDzeros=0; }
  if (fUsedD0bars) { delete fUsedD0bars; fUsedD0bars=0; }
  if (fUsedDstars) { delete fUsedDstars; fUsedDstars=0; }

  if (fListCutsDmesons) { delete fListCutsDmesons; fListCutsDmesons=0; }
  if (fListDzeroHistos) { delete fListDzeroHistos; fListDzeroHistos=0; }
  if (fListDstarHistos) { delete fListDstarHistos; fListDstarHistos=0; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::Init()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::Init()
//

  if (fListCutsDmesons) return;
  fListCutsDmesons = new TList(); fListCutsDmesons->SetOwner();

  if (fCutDzero) {
    AliRDHFCutsD0toKpi *cutDzero = new AliRDHFCutsD0toKpi(*fCutDzero);
    cutDzero->SetName("AnalysisCutsDzero"); 
    fListCutsDmesons->Add(cutDzero);
  }

  if (fCutDstar) {
    AliRDHFCutsDStartoKpipi *cutDstar = new AliRDHFCutsDStartoKpipi(*fCutDstar);
    cutDstar->SetName("AnalysisCutsDstar"); 
    fListCutsDmesons->Add(cutDstar);
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::UserCreateOutputObjects()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::UserCreateOutputObjects
//

  if (fCutDzero) { fUsedDzeros = new TClonesArray("AliAODRecoDecayHF2Prong",0); fUsedDzeros->SetName("AnaUsedDzero");}
  if (fCutDzero) { fUsedD0bars = new TClonesArray("AliAODRecoDecayHF2Prong",0); fUsedD0bars->SetName("AnaUsedD0bar");}
  if (fCutDstar) { fUsedDstars = new TClonesArray("AliAODRecoCascadeHF",    0); fUsedDstars->SetName("AnaUsedDstar"); }

  if (!fListDzeroHistos) fListDzeroHistos = new TList(); fListDzeroHistos->SetOwner();
  if (!fListDstarHistos) fListDstarHistos = new TList(); fListDstarHistos->SetOwner();

  if (fCutDzero) CreateDzeroHistograms();
  if (fCutDstar) CreateDstarHistograms();

  PostData(1, fListCutsDmesons);
  PostData(2, fListDzeroHistos);
  PostData(3, fListDstarHistos);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::UserExec(Option_t *)
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::UserExec
//

  if (fIsNotExecOnce) {
    fIsNotExecOnce=ExecOnce();
    if (fIsNotExecOnce) return;
  }

  if (ExecEach()) ExecAnas();

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::Terminate(Option_t *)
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::Terminate
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::NotifyRun()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::NotifyRun
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ExecAnas()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::ExecAnas
//

  if (fUsedDzeros) fUsedDzeros->Delete();
  if (fUsedD0bars) fUsedD0bars->Delete();
  if (fUsedDstars) fUsedDstars->Delete();

  if (fDzeroClArr) ExecAnasDzero();
  if (fDstarClArr) ExecAnasDstar();
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ExecAnasDzero()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::ExecAnasDzero
//

  if (!fCutDzero->IsEventSelected(fEventAOD)) return;
  const Int_t nCands = fDzeroClArr->GetEntriesFast(); if (nCands==0) return;

  ((TH1D*)fListDzeroHistos->FindObject("hDzero_candi"))->Fill(nCands);

  Int_t countN = 0;
  AliAODRecoDecayHF2Prong *candi = 0;
  for (Int_t iCand=0; iCand<nCands; iCand++) {
    candi = dynamic_cast<AliAODRecoDecayHF2Prong*>(fDzeroClArr->At(iCand)); if (!candi) continue;

    Double_t dPt = candi->Pt();
    if (!fCutDzero->IsInFiducialAcceptance(dPt,candi->YD0())) continue;
    Int_t mask = fCutDzero->IsSelected(candi,AliRDHFCuts::kAll,fEventAOD); if (mask==0) continue;

    Double_t dMass = 0.; 
    if ((mask==1) || (mask==3)) dMass = candi->InvMassD0();
    if  (mask==2)               dMass = candi->InvMassD0bar();
    ((TH2D*)fListDzeroHistos->FindObject("hDzero_Mass_Pt"))->Fill(dMass,dPt);

    if ((mask==1) || (mask==3)) new((*fUsedDzeros)[countN++]) AliAODRecoDecayHF2Prong(*candi);
    if  (mask==2)               new((*fUsedD0bars)[countN++]) AliAODRecoDecayHF2Prong(*candi);
    candi = 0;
  }

  ((TH1D*)fListDzeroHistos->FindObject("hDzero_seled"))->Fill(countN);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ExecAnasDstar()
{
//
//  AliAnalysisTaskSEDmesonsFilterCJ::ExecAnasDstar
//

  if (!fCutDstar->IsEventSelected(fEventAOD)) return;
  const Int_t nCands = fDstarClArr->GetEntriesFast(); if (nCands==0) return;

  Double_t dMassDzero = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t dMassDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t dMassDelta = dMassDstar - dMassDzero;
  Double_t dSigmaD0Pt[13]={ 0.012, 0.012, 0.012, 0.015, 0.015, 0.018, 0.018, 0.020, 0.020, 0.030, 0.030, 0.037, 0.040 };

  ((TH1D*)fListDstarHistos->FindObject("hDstar_candi"))->Fill(nCands);

  Int_t countN = 0;
  AliAODRecoCascadeHF *candi = 0;
  for (Int_t iCand=0; iCand<nCands; iCand++) {
    candi = dynamic_cast<AliAODRecoCascadeHF*>(fDstarClArr->At(iCand)); if (!candi) continue;

    Double_t dPt = candi->Pt();
    if (!fCutDstar->IsInFiducialAcceptance(dPt,candi->YDstar()))   continue;
    if (!fCutDstar->IsSelected(candi,AliRDHFCuts::kAll)) continue;

    Int_t iBin = fCutDzero->PtBin(dPt);
    if ((iBin<0) || (iBin>=fCutDzero->GetNPtBins())) { AliWarning(Form("pT(D*)=%f out of bounds!!!",dPt)); continue; }

    Double_t dMassResD0 = TMath::Abs(candi->InvMassD0() - dMassDzero); if (dMassResD0>3.*dSigmaD0Pt[iBin]) continue;

    Double_t deltaM = candi->DeltaInvMass();
    ((TH2D*)fListDstarHistos->FindObject("hDstar_deltaMass_Pt"))->Fill(deltaM,dPt);

    if (TMath::Abs(deltaM-dMassDelta)) ((TH1D*)fListDstarHistos->FindObject("hDstar_SoftPiPt"))->Fill(((AliAODTrack*)candi->GetBachelor())->Pt());

    new((*fUsedDstars)[countN++]) AliAODRecoCascadeHF(*candi);
    candi = 0;
  }

  ((TH1D*)fListDstarHistos->FindObject("hDstar_seled"))->Fill(countN);
  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::ExecOnce()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::ExecOnce
//

  if (!InputEvent())  return kTRUE;
  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());

  if (fEventAOD) {
    if (fCutDzero) fDzeroClArr = dynamic_cast<TClonesArray*>(fEventAOD->FindListObject("D0toKpi"));
    if (fCutDstar) fDstarClArr = dynamic_cast<TClonesArray*>(fEventAOD->FindListObject("Dstar"));
  } else if (AODEvent() && IsStandardAOD()) {
    fEventAOD = dynamic_cast<AliAODEvent*>(AODEvent());
    AliAODHandler* aodH = dynamic_cast<AliAODHandler*>((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodH->GetExtensions()) {
      AliAODExtension    *aodExt = dynamic_cast<AliAODExtension*>(aodH->GetExtensions()->FindObject("AliAOD.VertexingHF.root"));
      if (fCutDzero) fDzeroClArr = dynamic_cast<TClonesArray*>(aodExt->GetAOD()->FindListObject("D0toKpi"));
      if (fCutDstar) fDstarClArr = dynamic_cast<TClonesArray*>(aodExt->GetAOD()->FindListObject("Dstar"));
    }
  }

  if (!fEventAOD)                { AliError("No AOD Event!!!");    return kTRUE; }
  if (!fDzeroClArr && fCutDzero) { AliError("No Dzero Branch!!!"); return kTRUE; }
  if (!fDstarClArr && fCutDstar) { AliError("No Dstar Branch!!!"); return kTRUE; }

  if (fDzeroClArr && fUsedDzeros) {
    fUsedDzeros->Delete();
    if (!(InputEvent()->FindListObject("AnaUsedDzero"))) InputEvent()->AddObject(fUsedDzeros);
  }

  if (fDzeroClArr && fUsedD0bars) {
    fUsedD0bars->Delete();
    if (!(InputEvent()->FindListObject("AnaUsedD0bar"))) InputEvent()->AddObject(fUsedD0bars);
  }

  if (fDstarClArr && fUsedDstars) {
    fUsedDstars->Delete();
    if (!(InputEvent()->FindListObject("AnaUsedDstar"))) InputEvent()->AddObject(fUsedDstars);
  }

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::ExecEach()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::ExecEach
//

  if (!fEventAOD) {
    AliWarning("No Input Event, Skip Event!!!");
    return kFALSE;
  }

  if ((!fDzeroClArr) && (!fDstarClArr)) {
    AliWarning("No D meson Branch, Skip Event!!!");
    return kFALSE;
  }

  if ((!fEventAOD->GetPrimaryVertex()) || (TMath::Abs(fEventAOD->GetMagneticField())<1e-3)) return kFALSE;

  Bool_t bIsDmeson = kTRUE;
  if (fDzeroClArr) bIsDmeson = ((fDzeroClArr->GetEntriesFast()==0) && bIsDmeson);
  if (fDstarClArr) bIsDmeson = ((fDstarClArr->GetEntriesFast()==0) && bIsDmeson);
  if (bIsDmeson)   return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::CreateDzeroHistograms()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::CreateDzeroHistograms
//

  if (!fListDzeroHistos) return;
  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1D *h1D = 0;
  h1D = new TH1D("hDzero_candi", "", 100, -0.5, 99.5); h1D->Sumw2(); fListDzeroHistos->Add(h1D); h1D=0;
  h1D = new TH1D("hDzero_seled", "", 100, -0.5, 99.5); h1D->Sumw2(); fListDzeroHistos->Add(h1D); h1D=0;

  TH2D *h2D = 0;
  Double_t dMass = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  h2D = new TH2D("hDzero_Mass_Pt", "", 200, dMass-0.15, dMass+0.15, 100, 0., 50.);
  h2D->Sumw2(); fListDzeroHistos->Add(h2D); h2D=0;

  TH1::AddDirectory(bStatusTmpH);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::CreateDstarHistograms()
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::CreateDstarHistograms
//

  if (!fListDstarHistos) return;
  Bool_t bStatusTmpH = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1D *h1D = 0;
  h1D = new TH1D("hDstar_candi",    "", 100, -0.5, 99.5); h1D->Sumw2(); fListDstarHistos->Add(h1D); h1D=0;
  h1D = new TH1D("hDstar_seled",    "", 100, -0.5, 99.5); h1D->Sumw2(); fListDstarHistos->Add(h1D); h1D=0;
  h1D = new TH1D("hDstar_SoftPiPt", "", 500,  0.,  10.0); h1D->Sumw2(); fListDstarHistos->Add(h1D); h1D=0;

  TH2D *h2D = 0;
  Double_t dMassDzero = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t dMassDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t dMassDelta = dMassDstar - dMassDzero;
  h2D = new TH2D("hDstar_deltaMass_Pt", "", 200, dMassDelta-0.015, dMassDelta+0.015, 100, 0., 50.);
  h2D->Sumw2(); fListDstarHistos->Add(h2D); h2D=0;

  TH1::AddDirectory(bStatusTmpH);
  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::SetCutDzero(AliRDHFCutsD0toKpi *cut)
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::SetCutDzero
//

  if (!cut) return kTRUE;

  if (fCutDzero) { delete fCutDzero; fCutDzero=0; }
  fCutDzero = new AliRDHFCutsD0toKpi(*cut);
  if (!fCutDzero) { AliError("No Dzero Cuts!!!"); return kTRUE; }

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::SetCutDstar(AliRDHFCutsDStartoKpipi *cut)
{
//
// AliAnalysisTaskSEDmesonsFilterCJ::SetCutDstar
//

  if (!cut) return kTRUE;

  if (fCutDstar) { delete fCutDstar; fCutDstar=0; }
  fCutDstar = new AliRDHFCutsDStartoKpipi(*cut);
  if (!fCutDstar)  { AliError("No Dstar Cuts!!!"); return kTRUE; }

  return kFALSE;
}
