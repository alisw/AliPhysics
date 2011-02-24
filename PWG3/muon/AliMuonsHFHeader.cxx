/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// class used to extract and store info at event level
//
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TH1.h>
#include <TH1D.h>
#include <TList.h>

#include "AliVEvent.h"
#include "AliMuonInfoStoreRD.h"
#include "AliMuonInfoStoreMC.h"
#include "AliDimuInfoStoreRD.h"
#include "AliDimuInfoStoreMC.h"
#include "AliMuonsHFHeader.h"

class TNamed;
class AliVVertex;

ClassImp(AliMuonsHFHeader)

const TString AliMuonsHFHeader::fgkStdBranchName("MuEvsH");
Int_t         AliMuonsHFHeader::fgAnaMode = 0;
Bool_t        AliMuonsHFHeader::fgIsMC    = kFALSE;
Double_t      AliMuonsHFHeader::fgCuts[5] = { -999999., 999999., 999999., -999999., 999999. };

//_____________________________________________________________________________
AliMuonsHFHeader::AliMuonsHFHeader() :
TNamed(),
fSelMask(AliVEvent::kAny),
fIsMB(kFALSE),
fIsMU(kFALSE),
fVtxContrsN(0),
fFiredTriggerClass(),
fCentrality(0.)
{
  //
  // default constructor
  //
  for (Int_t i=3; i--;) fVtx[i] = 0.;
}

//_____________________________________________________________________________
AliMuonsHFHeader::AliMuonsHFHeader(const AliMuonsHFHeader &src) :
TNamed(),
fSelMask(src.fSelMask),
fIsMB(src.fIsMB),
fIsMU(src.fIsMU),
fVtxContrsN(src.fVtxContrsN),
fFiredTriggerClass(src.fFiredTriggerClass),
fCentrality(src.fCentrality)
{
  //
  // copy constructor
  //
  for (Int_t i=3; i--;) fVtx[i] = src.fVtx[i];
}

//_____________________________________________________________________________
AliMuonsHFHeader& AliMuonsHFHeader::operator=(const AliMuonsHFHeader &src)
{
  //
  // assignment constructor
  //

  fSelMask           = src.fSelMask;
  fIsMB              = src.fIsMB;
  fIsMU              = src.fIsMU;
  fVtxContrsN        = src.fVtxContrsN;
  fFiredTriggerClass = src.fFiredTriggerClass;
  fCentrality        = src.fCentrality;
  for (Int_t i=3; i--;) fVtx[i] = src.fVtx[i];

  return *this;
}

//_____________________________________________________________________________
AliMuonsHFHeader::~AliMuonsHFHeader()
{
  //
  // default destructor
  //
}

//_____________________________________________________________________________
void AliMuonsHFHeader::SetVertex(AliVVertex *vertex)
{
  // extract event info from AOD event

  vertex->GetXYZ(fVtx);
  fVtxContrsN = vertex->GetNContributors();
  this->SetTitle(vertex->GetTitle());
  return;
}

void AliMuonsHFHeader::SetFiredTriggerClass(TString trigger)
{
  fFiredTriggerClass = trigger;
  fIsMB=(fFiredTriggerClass.Contains("CINT1B-ABCE-NOPF-ALL")  ||   // p-p   min. bias trigger before Aug/2010
         fFiredTriggerClass.Contains("CINT1-B-NOPF-ALLNOTRD") ||   // p-p   min. bias trigger from   Aug/2010
         fFiredTriggerClass.Contains("CMBS2A-B-NOPF-ALL")     ||   // Pb-Pb min. bias trigger 2-out-of-3
         fFiredTriggerClass.Contains("CMBS2C-B-NOPF-ALL")     ||   // Pb-Pb min. bias trigger 2-out-of-3
         fFiredTriggerClass.Contains("CMBAC-B-NOPF-ALL")      ||   // Pb-Pb min. bias trigger 2-out-of-3
         fFiredTriggerClass.Contains("CMBACS2-B-NOPF-ALL")    ||   // Pb-Pb min. bias trigger 3-out-of-3 (early)
         fFiredTriggerClass.Contains("CMBACS2-B-NOPF-ALLNOTRD"));  // Pb-Pb min. bias trigger 3-out-of-3 (late 2010)
  fIsMU=(fFiredTriggerClass.Contains("CMUS1B-ABCE-NOPF-MUON") ||   // p-p MUON trigger before Aug/2010
         fFiredTriggerClass.Contains("CMUS1-B-NOPF-ALLNOTRD"));    // p-p MUON trigger from   Aug/2010 
  return;
}

//_____________________________________________________________________________
Bool_t AliMuonsHFHeader::IsSelected()
{
  // select event according to the event selection cuts
  if (this->VtxContrsN()<fgCuts[0])       return kFALSE;
  if (TMath::Abs(this->Vz())>fgCuts[1])   return kFALSE;
  if (this->Vt()>fgCuts[2])               return kFALSE;

  // centrality selection
  Float_t centr = this->Centrality();
  if (centr<fgCuts[3] || centr>fgCuts[4]) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistograms(TList *list)
{
  // create output histos of muon analysis according to the analysis mode & MC flag

  if (fgIsMC) {
    this->CreateHistosEvnH(list);
    if (fgAnaMode!=2) {
      TString sName[7] = { "Unidentified", "Hadron", "SecondaryMu", "PrimaryMu", "CharmMu", "BottomMu", "" };
      for (Int_t i=7; i--;) this->CreateHistosMuon(list, sName[i]);
    }
    if (fgAnaMode!=1) {
      TString sName[7] = { "Uncorr", "Resonance", "DDsame", "DDdiff", "BBsame", "BBdiff", "" };
      for (Int_t i=7; i--;) this->CreateHistosDimu(list, sName[i]);
    }
    return;
  }

  this->CreateHistosEvnH(list,"MB"); this->CreateHistosEvnH(list,"MU");
  if (fgAnaMode!=2) { this->CreateHistosMuon(list,"MB"); this->CreateHistosMuon(list,"MU"); }
  if (fgAnaMode!=1) { this->CreateHistosDimu(list,"MB"); this->CreateHistosDimu(list,"MU"); }
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosEvnH(TList *list, TString sName)
{
  // create histograms at event level

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nhs    = 3;
  TString tName[nhs] = {  "Vz",  "Vt",  "VtxNcontr" };
  Int_t   nbins[nhs] = {  800 ,   40 ,        202   };
  Double_t xlow[nhs] = {  -40.,    0.,         -2.5 };
  Double_t  xup[nhs] = {   40.,    4.,        199.5 };

  TH1D *histo = 0;
  for (Int_t i=0; i<nhs; i++) {
    char *hName = Form("h%s_%s", sName.Data(), tName[i].Data());
    histo = new TH1D(hName, hName, nbins[i], xlow[i], xup[i]);
    histo->Sumw2(); list->Add(histo); histo = 0;
  }

  TH1::AddDirectory(oldStatus);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosMuon(TList *list, TString sName)
{
  // create histograms for single muon

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nhs    = 7;
  TString tName[nhs] = {   "P",  "Pt",  "Eta",  "DCA",  "TrM",  "Charge", "Rabs" };
  Int_t   nbins[nhs] = { 1500 ,  300 ,    15 ,  1000 ,    4  ,       3  ,    48  };
  Double_t xlow[nhs] = {    0.,    0.,   -4.0,     0.,   -0.5,      -1.5,   17.6 };
  Double_t  xup[nhs] = {  150.,   30.,   -2.5,   500.,    3.5,       1.5,   80.0 };

  TH1D *histo = 0;
  for (Int_t i=0; i<nhs; i++) {
    char *hName = Form("h%s_%s", sName.Data(), tName[i].Data());
    histo = new TH1D(hName, hName, nbins[i], xlow[i], xup[i]);
    histo->Sumw2(); list->Add(histo); histo = 0;
  }

  TH1::AddDirectory(oldStatus);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosDimu(TList *list, TString sName)
{
  // create histograms for dimuon

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1D *histo = 0;
  const Int_t nhs    = 3;
  TString tName[nhs] = {   "P",  "Pt",  "InvM"   };
  Int_t   nbins[nhs] = { 1500 ,  300 ,    300    };
  Double_t xlow[nhs] = {    0.,    0.,      0.   };
  Double_t  xup[nhs] = {  150.,   30.,     30.   };
  TString dimuName[3] = { "DimuNN", "DimuNP", "DimuPP" };
  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<nhs; j++) {
      char *hName = Form("h%s_%s_%s", sName.Data(), dimuName[i].Data(), tName[j].Data());
      histo = new TH1D(hName, hName, nbins[j], xlow[j], xup[j]);
      histo->Sumw2(); list->Add(histo); histo = 0;
    }
  }

  TH1::AddDirectory(oldStatus);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosEvnH(TList *list)
{
  // fill histograms at event level according to event selection cuts

  if (!list)                   return;
  if (!this->IsSelected()) return;

  const Int_t nhs    = 3;
  TString tName[nhs] = {       "Vz",       "Vt",        "VtxNcontr" };
  Double_t dist[nhs] = { this->Vz(), this->Vt(), this->VtxContrsN() };
  if (fgIsMC && (fSelMask & AliVEvent::kAny)) {
    for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h_%s",tName[i].Data())))->Fill(dist[i]);
  } else {
    if (fIsMB && (fSelMask & AliVEvent::kMB))   { for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s","MB",tName[i].Data())))->Fill(dist[i]); }
    if (fIsMU && (fSelMask & AliVEvent::kMUON)) { for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s","MU",tName[i].Data())))->Fill(dist[i]); }
  }
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosMuon(TList *list, AliMuonInfoStoreRD* const infoStore, Int_t s)
{
  // fill histograms for single muon according to event & muon track selection cuts

  if (!list)                    return;
  if (!this->IsSelected())      return;
  if (!infoStore->IsSelected()) return;

  const Int_t nhs    = 7;
  TString tName[nhs] = { "P", "Pt", "Eta", "DCA", "TrM", "Charge", "Rabs" };
  Double_t dist[nhs] = { infoStore->Momentum().Mag(),
                             infoStore->Momentum().Pt(),
                             infoStore->Momentum().Eta(),
                             infoStore->DCA(),
                             infoStore->MatchTrigger(),
                             infoStore->Charge(),
                             infoStore->RabsEnd() };

  if (fgIsMC && (fSelMask & AliVEvent::kAny)) {
    TString sName[7] = { "BottomMu", "CharmMu", "PrimaryMu", "SecondaryMu", "Hadron", "Unidentified", "" };
    for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s",sName[6].Data(),tName[i].Data())))->Fill(dist[i]);
    for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s",sName[s].Data(),tName[i].Data())))->Fill(dist[i]);
  } else {
    if (fIsMB && (fSelMask & AliVEvent::kMB))   { for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s","MB",tName[i].Data())))->Fill(dist[i]); }
    if (fIsMU && (fSelMask & AliVEvent::kMUON)) { for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s","MU",tName[i].Data())))->Fill(dist[i]); }
  }

  return; 
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosDimu(TList *list, AliDimuInfoStoreRD* const infoStore, Int_t s)
{
  // fill histograms for dimuon according to evnet & dimuon candidates selection cuts

  if (!list)                    return;
  if (!this->IsSelected())      return;
  if (!infoStore->IsSelected()) return;

  TString dimuName = "DimuNN";
  if (infoStore->Charge()==0)     dimuName = "DimuNP";
  else if (infoStore->Charge()>0) dimuName = "DimuPP";

  const Int_t nhs    = 3;
  TString tName[nhs] = { "P", "Pt", "InvM" };
  Double_t dist[nhs] = { infoStore->Momentum().Mag(),
                             infoStore->Momentum().Pt(),
                             infoStore->InvM() };

  if (fgIsMC && (fSelMask & AliVEvent::kAny)) {
    TString sName[7] = { "BBdiff", "BBsame", "DDdiff", "DDsame", "Resonance", "Uncorr", "" };
    for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s_%s",sName[6].Data(),dimuName.Data(),tName[i].Data())))->Fill(dist[i]);
    for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s_%s",sName[s].Data(),dimuName.Data(),tName[i].Data())))->Fill(dist[i]);
  } else {
    if (fIsMB && (fSelMask & AliVEvent::kMB)) {
      for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s_%s","MB",dimuName.Data(),tName[i].Data())))->Fill(dist[i]);
    }
    if (fIsMU && (fSelMask & AliVEvent::kMUON)) {
      for (Int_t i=nhs; i--;) ((TH1D*)list->FindObject(Form("h%s_%s_%s","MU",dimuName.Data(),tName[i].Data())))->Fill(dist[i]);
    }
  }

  return;
}
