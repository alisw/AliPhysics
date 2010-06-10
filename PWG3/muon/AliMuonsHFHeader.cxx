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

/////////////////////////////////////////////////////////////
//
// class used to extract and store info at event level
//
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TList.h>

#include "AliMuonInfoStoreRD.h"
#include "AliMuonInfoStoreMC.h"
#include "AliDimuInfoStoreRD.h"
#include "AliDimuInfoStoreMC.h"
#include "AliMuonsHFHeader.h"

class TNamed;
class AliESDVertex;

ClassImp(AliMuonsHFHeader)

const TString AliMuonsHFHeader::fgkStdBranchName("MuEvsH");
Int_t         AliMuonsHFHeader::fgAnaMode = 0;
Bool_t        AliMuonsHFHeader::fgIsMC    = kFALSE;
Double_t      AliMuonsHFHeader::fgCuts[3] = { -999999., 999999., 999999.};

//_____________________________________________________________________________
AliMuonsHFHeader::AliMuonsHFHeader() :
TNamed(),
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
void AliMuonsHFHeader::SetEvent(AliVVertex *vertex)
{
  // extract event info from AOD event

  vertex->GetXYZ(fVtx);
  fVtxContrsN = vertex->GetNContributors();
  this->SetTitle(vertex->GetTitle());
  return;
}

//_____________________________________________________________________________
Bool_t AliMuonsHFHeader::EventSelection()
{
  // select event according to the event selection cuts
  if (this->VtxContrsN()<fgCuts[0])     return kFALSE;
  if (TMath::Abs(this->Vz())>fgCuts[1]) return kFALSE;
  if (this->Vt()>fgCuts[2])             return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistograms(TList *list)
{
  // create output histos of muon analysis according to the analysis mode & MC flag

  this->CreateHistosEvnH(list);

  if (fgAnaMode!=2) {
   this->CreateHistosMuon(list);
   if (fgIsMC) {
      TString sName[6] = { "Unidentified", "Hadron", "SecondaryMu", "PrimaryMu", "CharmMu", "BottomMu" };
      for (Int_t i=AliMuonInfoStoreMC::SourcesN(); i--;) this->CreateHistosMuon(list, sName[i]);
    }
  }

  if (fgAnaMode!=1) {
    this->CreateHistosDimu(list);
    if (fgIsMC) {
      TString sName[6] = { "Uncorr", "Resonance", "DDsame", "DDdiff", "BBsame", "BBdiff" };
      for (Int_t i=AliDimuInfoStoreMC::SourcesN(); i--;) this->CreateHistosDimu(list, sName[i]);
    }
  }

  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosEvnH(TList *list)
{
  // create histograms at event level

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nHistos = 3;
  TString tname[nHistos] = { "hVz", "hVt", "hVtxNcontr" };
  Int_t   nbins[nHistos] = {  800 ,   40 ,        202   };
  Double_t xlow[nHistos] = {  -40.,    0.,         -2.5 };
  Double_t  xup[nHistos] = {   40.,    4.,        199.5 };

  TH1F *histo = 0;
  for (Int_t i=0; i<nHistos; i++) {
    histo = new TH1F(tname[i], tname[i].Data(), nbins[i], xlow[i], xup[i]);
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
  
  const Int_t nHistos = 8;
  TString tName[nHistos] = {   "P",  "Pt",  "Eta",  "DCA",  "TrM",  "Charge", "Rabs", "UnfVtx" };
  Int_t   nbins[nHistos] = { 1500 ,  300 ,   100 ,   500 ,    4  ,       3  ,   300 ,      80  };
  Double_t xlow[nHistos] = {    0.,    0.,   -10.,     0.,   -0.5,      -1.5,     0.,     -40. };
  Double_t  xup[nHistos] = {  150.,   30.,     0.,   500.,    3.5,       1.5,   150.,      40. };

  TH1F *histo = 0;
  for (Int_t i=0; i<nHistos; i++) {
    char *hName = Form("h%s_%s", sName.Data(), tName[i].Data());
    histo = new TH1F(hName, hName, nbins[i], xlow[i], xup[i]);
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

  TH1F *histo = 0;
  const Int_t nHistos = 3;
  TString tName[nHistos] = {   "P",  "Pt",  "InvM"   };
  Int_t   nbins[nHistos] = { 1500 ,  300 ,    300    };
  Double_t xlow[nHistos] = {    0.,    0.,      0.   };
  Double_t  xup[nHistos] = {  150.,   30.,     30.   };
  TString dimuName[3] = { "DimuNN", "DimuNP", "DimuPP" };
  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<nHistos; j++) {
      char *hName = Form("h%s_%s_%s", sName.Data(), dimuName[i].Data(), tName[j].Data());
      histo = new TH1F(hName, hName, nbins[j], xlow[j], xup[j]);
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
  if (!this->EventSelection()) return;

  const Int_t    nHistos = 3;
  TString tname[nHistos] = {      "hVz",      "hVt",       "hVtxNcontr" };
  Double_t dist[nHistos] = { this->Vz(), this->Vt(), this->VtxContrsN() };
  for (Int_t i=nHistos; i--;) ((TH1F*)list->FindObject(tname[i].Data()))->Fill(dist[i]);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosMuon(TList *list, AliMuonInfoStoreRD* const infoStore, Int_t src)
{
  // fill histograms for single muon according to event & muon track selection cuts

  if (!list)                       return;
  if (!this->EventSelection())     return;
  if (!infoStore->MuonSelection()) return;

  const Int_t nHistos    = 8;
  TString tName[nHistos] = { "P", "Pt", "Eta", "DCA", "TrM", "Charge", "Rabs", "UnfVtx" };
  Double_t dist[nHistos] = { infoStore->Momentum().Mag(),
                             infoStore->Momentum().Pt(),
                             infoStore->Momentum().Eta(),
                             infoStore->DCA(),
                             infoStore->MatchTrigger(),
                             infoStore->Charge(),
                             infoStore->RabsEnd(),
                             this->Vz() };
  for (Int_t i=nHistos; i--;) ((TH1F*)list->FindObject(Form("h%s_%s","",tName[i].Data())))->Fill(dist[i]);

  if (fgIsMC && src>=0) {
    TString sName[6] = { "BottomMu", "CharmMu", "PrimaryMu", "SecondaryMu", "Hadron", "Unidentified" };
    for (Int_t i=nHistos; i--;) ((TH1F*)list->FindObject(Form("h%s_%s",sName[src].Data(),tName[i].Data())))->Fill(dist[i]);
  }

  return; 
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosDimu(TList *list, AliDimuInfoStoreRD* const infoStore, Int_t src)
{
  // fill histograms for dimuon according to evnet & dimuon candidates selection cuts

  if (!list)                       return;
  if (!this->EventSelection())     return;
  if (!infoStore->DimuSelection()) return;

  TString dimuName = "DimuNN";
  if (infoStore->Charge()==0)     dimuName = "DimuNP";
  else if (infoStore->Charge()>0) dimuName = "DimuPP";

  const Int_t nHistos    = 3;
  TString tName[nHistos] = { "P", "Pt", "InvM" };
  Double_t dist[nHistos] = { infoStore->Momentum().Mag(),
                             infoStore->Momentum().Pt(),
                             infoStore->InvM() };
  for (Int_t i=nHistos; i--;) ((TH1F*)list->FindObject(Form("h%s_%s_%s","",dimuName.Data(),tName[i].Data())))->Fill(dist[i]);

  if (fgIsMC && src>=0) {
    TString sName[6] = { "BBdiff", "BBsame", "DDdiff", "DDsame", "Resonance", "Uncorr" };
    for (Int_t i=nHistos; i--;) ((TH1F*)list->FindObject(Form("h%s_%s_%s",sName[src].Data(),dimuName.Data(),tName[i].Data())))->Fill(dist[i]);
  }

  return;
}
