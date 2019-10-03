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

#include "AliVEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVVertex.h"
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliAnalysisUtils.h"
#include "AliMultiplicity.h"
#include "AliMuonTrackCuts.h"
#include "AliMuonInfoStoreRD.h"
#include "AliMuonInfoStoreMC.h"
#include "AliDimuInfoStoreRD.h"
#include "AliDimuInfoStoreMC.h"
#include "AliMuonsHFHeader.h"

#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

ClassImp(AliMuonsHFHeader)

const TString AliMuonsHFHeader::fgkStdBranchName("MuEvsH");
Double_t      AliMuonsHFHeader::fgCuts[5] = { -999999., 999999., 999999., -999999., 999999. };

//_____________________________________________________________________________
AliMuonsHFHeader::AliMuonsHFHeader() :
TNamed(),
fAnaMode(0),
fIsMC(kFALSE),
fSelMask(AliVEvent::kAny),
fIsMB(kFALSE),
fIsMU(kFALSE),
fVtxContrsN(0),
fFiredTriggerClass(),
fCentQA(-1),
fEventPlane(0.),
fIsEvtInChunk(kFALSE),
fIsVtxSeled2013pA(kFALSE),
fNumOfTrklets(-1),
fTrgInpts(0),
fImpParam(-1.),
fCentralityV0M(-1.),
fCentralityV0A(-1.),
fCentralityV0C(-1.),
fCentralityCL1(-1.),
fCentralityZNA(-1.),
fCentralityZNC(-1.),
fPUMask(0)
{
  //
  // default constructor
  //
  for (Int_t i=3; i--;) fVtx[i] = 0.;
  for (Int_t i=3; i--;) fVMC[i] = 0.;
}

//_____________________________________________________________________________
AliMuonsHFHeader::AliMuonsHFHeader(const AliMuonsHFHeader &src) :
TNamed(),
fAnaMode(src.fAnaMode),
fIsMC(src.fIsMC),
fSelMask(src.fSelMask),
fIsMB(src.fIsMB),
fIsMU(src.fIsMU),
fVtxContrsN(src.fVtxContrsN),
fFiredTriggerClass(src.fFiredTriggerClass),
fCentQA(src.fCentQA),
fEventPlane(src.fEventPlane),
fIsEvtInChunk(src.fIsEvtInChunk),
fIsVtxSeled2013pA(src.fIsVtxSeled2013pA),
fNumOfTrklets(src.fNumOfTrklets),
fTrgInpts(src.fTrgInpts),
fImpParam(src.fImpParam),
fCentralityV0M(src.fCentralityV0M),
fCentralityV0A(src.fCentralityV0A),
fCentralityV0C(src.fCentralityV0C),
fCentralityCL1(src.fCentralityCL1),
fCentralityZNA(src.fCentralityZNA),
fCentralityZNC(src.fCentralityZNC),
fPUMask(src.fPUMask)
{
  //
  // copy constructor
  //
  for (Int_t i=3; i--;) fVtx[i] = src.fVtx[i];
  for (Int_t i=3; i--;) fVMC[i] = src.fVMC[i];
}

//_____________________________________________________________________________
AliMuonsHFHeader& AliMuonsHFHeader::operator=(const AliMuonsHFHeader &src)
{
  //
  // assignment constructor
  //

  if(&src==this) return *this;

  fAnaMode           = src.fAnaMode;
  fIsMC              = src.fIsMC;
  fSelMask           = src.fSelMask;
  fIsMB              = src.fIsMB;
  fIsMU              = src.fIsMU;
  fVtxContrsN        = src.fVtxContrsN;
  fFiredTriggerClass = src.fFiredTriggerClass;
  fCentQA            = src.fCentQA;
  fEventPlane        = src.fEventPlane;
  fIsEvtInChunk      = src.fIsEvtInChunk;
  fIsVtxSeled2013pA  = src.fIsVtxSeled2013pA;
  fNumOfTrklets      = src.fNumOfTrklets;
  fTrgInpts          = src.fTrgInpts;
  fImpParam          = src.fImpParam;
  fCentralityV0M     = src.fCentralityV0M;
  fCentralityV0A     = src.fCentralityV0A;
  fCentralityV0C     = src.fCentralityV0C; 
  fCentralityCL1     = src.fCentralityCL1; 
  fCentralityZNA     = src.fCentralityZNA;
  fCentralityZNC     = src.fCentralityZNC; 
  fPUMask            = src.fPUMask;
  for (Int_t i=3; i--;) fVtx[i] = src.fVtx[i];
  for (Int_t i=3; i--;) fVMC[i] = src.fVMC[i];

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
void AliMuonsHFHeader::SetEventInfo(AliVEventHandler* const handler)
{
  // fill info at event level

  AliVEvent   *event = handler->GetEvent();
  AliAODEvent *aod   = dynamic_cast<AliAODEvent*>(event);
  AliESDEvent *esd   = dynamic_cast<AliESDEvent*>(event);

  fSelMask = handler->IsEventSelected();
  if (aod) { fFiredTriggerClass = aod->GetFiredTriggerClasses(); fTrgInpts = aod->GetHeader()->GetL0TriggerInputs(); }
  if (esd) { fFiredTriggerClass = esd->GetFiredTriggerClasses(); fTrgInpts = esd->GetHeader()->GetL0TriggerInputs(); }
  fIsMB = fSelMask & AliVEvent::kMB;
  fIsMU = fSelMask & AliVEvent::kMUON;

  const AliVVertex *vertex = event->GetPrimaryVertex();
  vertex->GetXYZ(fVtx);
  fVtxContrsN = vertex->GetNContributors();
  if (fIsMC) {
    AliMCEvent *mcEvent = handler->MCEvent();
    AliGenEventHeader *genHeader = mcEvent->GenEventHeader();
    if (!genHeader) {  AliError("Header not found. Nothing done!"); return; }
    AliGenDPMjetEventHeader *dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);
    if (dpmHeader) fImpParam = dpmHeader->ImpactParameter();

    if (esd)   mcEvent->GetPrimaryVertex()->GetXYZ(fVMC);
    if (aod) ((AliAODMCHeader*)aod->FindListObject(AliAODMCHeader::StdBranchName()))->GetVertex(fVMC);
  } this->SetTitle(vertex->GetTitle());
  //(aod && !aod->GetTracklets()) ? event->IsPileupFromSPD(3,0.8,3.,2.,5.) : event->IsPileupFromSPDInMultBins();

  AliAnalysisUtils *anaUtils = new AliAnalysisUtils();
  if (esd) {
    fIsEvtInChunk     = anaUtils->IsFirstEventInChunk(esd);
    fIsVtxSeled2013pA = anaUtils->IsVertexSelected2013pA(esd);
    fNumOfTrklets     = esd->GetMultiplicity()->GetNumberOfTracklets();
  } 
  if (aod) {
    fIsEvtInChunk     = anaUtils->IsFirstEventInChunk(aod);
    fIsVtxSeled2013pA = anaUtils->IsVertexSelected2013pA(aod);
    fNumOfTrklets     = aod->GetTracklets()->GetNumberOfTracklets();
  }

  fPUMask = this->CollectPUMask(event);

  AliEventplane *evnP = event->GetEventplane();
  if (evnP) fEventPlane = evnP->GetEventplane("Q");
//if (evnP) fEventPlane = evnP->GetEventplane("V0A");

   AliCentrality *cent = event->GetCentrality(); if (cent) {
    fCentQA        = cent->GetQuality();
    fCentralityV0M = cent->GetCentralityPercentileUnchecked("V0M");
    fCentralityV0A = cent->GetCentralityPercentileUnchecked("V0A");
    fCentralityV0C = cent->GetCentralityPercentileUnchecked("V0C");
    fCentralityCL1 = cent->GetCentralityPercentileUnchecked("CL1");
    fCentralityZNA = cent->GetCentralityPercentileUnchecked("ZNA");
    fCentralityZNC = cent->GetCentralityPercentileUnchecked("ZNC");
  }

  aod = 0; esd = 0;
  delete anaUtils; anaUtils = 0;
  
  return;
}

//_____________________________________________________________________________
Bool_t AliMuonsHFHeader::IsSelected()
{
  // select event according to the event selection cuts
  if (this->VtxContrsN()<fgCuts[0])     return kFALSE;
  if (TMath::Abs(this->Vz())>fgCuts[1]) return kFALSE;
  if (this->Vt()>fgCuts[2])             return kFALSE;

  // centrality selection
  Float_t centrV0M = this->Centrality(kV0M);
  Float_t centrV0A = this->Centrality(kV0A);
  Float_t centrV0C = this->Centrality(kV0C); 
  Float_t centrCL1 = this->Centrality(kCL1); 
  Float_t centrZNA = this->Centrality(kZNA);
  Float_t centrZNC = this->Centrality(kZNC); 

  if (centrV0M<fgCuts[3] || centrV0M>fgCuts[4]) return kFALSE;
  if (centrV0A<fgCuts[3] || centrV0A>fgCuts[4]) return kFALSE;
  if (centrV0C<fgCuts[3] || centrV0C>fgCuts[4]) return kFALSE;
  if (centrCL1<fgCuts[3] || centrCL1>fgCuts[4]) return kFALSE;
  if (centrZNA<fgCuts[3] || centrZNA>fgCuts[4]) return kFALSE;
  if (centrZNC<fgCuts[3] || centrZNC>fgCuts[4]) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistograms(TList *list)
{
  // create output histos of muon analysis according to the analysis mode & MC flag

  if (fIsMC) {
    this->CreateHistosEvnH(list);
    if (fAnaMode!=2) {
      TString sName[7] = { "Unidentified", "Hadron", "SecondaryMu", "PrimaryMu", "CharmMu", "BottomMu", "" };
      for (Int_t i=7; i--;) this->CreateHistosMuon(list, sName[i]);
    }
    if (fAnaMode!=1) {
      TString sName[7] = { "Uncorr", "Resonance", "DDsame", "DDdiff", "BBsame", "BBdiff", "" };
      for (Int_t i=7; i--;) this->CreateHistosDimu(list, sName[i]);
    }
    return;
  }

  this->CreateHistosEvnH(list,"MB"); this->CreateHistosEvnH(list,"MU");
  if (fAnaMode!=2) { this->CreateHistosMuon(list,"MB"); this->CreateHistosMuon(list,"MU"); }
  if (fAnaMode!=1) { this->CreateHistosDimu(list,"MB"); this->CreateHistosDimu(list,"MU"); }

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

  if (!list)               return;
  if (!this->IsSelected()) return;

  const Int_t nhs    = 3;
  TString tName[nhs] = {       "Vz",       "Vt",        "VtxNcontr" };
  Double_t dist[nhs] = { this->Vz(), this->Vt(), static_cast<Double_t>(this->VtxContrsN()) };
  if (fIsMC && (fSelMask & AliVEvent::kAny)) {
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

  if (!list)                     return;
  if (!this->IsSelected())       return;
  if (!infoStore->IsSelected(0)) return;

  const Int_t nhs    = 7;
  TString tName[nhs] = { "P", "Pt", "Eta", "DCA", "TrM", "Charge", "Rabs" };
  Double_t dist[nhs] = { infoStore->MomentumAtVtx().Mag(),
                         infoStore->MomentumAtVtx().Pt(),
                         infoStore->MomentumAtVtx().Eta(),
                         infoStore->DCA(),
                         static_cast<Double_t>(infoStore->MatchTrigger()),
                         static_cast<Double_t>(infoStore->Charge()),
                         infoStore->RabsEnd() };

  if (fIsMC && (fSelMask & AliVEvent::kAny)) {
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

  if (!list)                     return;
  if (!this->IsSelected())       return;
  if (!infoStore->IsSelected(0)) return;

  TString dimuName = "DimuNN";
  if (infoStore->Charge()==0)     dimuName = "DimuNP";
  else if (infoStore->Charge()>0) dimuName = "DimuPP";

  const Int_t nhs    = 3;
  TString tName[nhs] = { "P", "Pt", "InvM" };
  Double_t dist[nhs] = { infoStore->Momentum().Mag(),
                         infoStore->Momentum().Pt(),
                         infoStore->InvM() };

  if (fIsMC && (fSelMask & AliVEvent::kAny)) {
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

//_____________________________________________________________________________
Float_t AliMuonsHFHeader::Centrality(Int_t centrality)
{
  // obtain centrality via selected estimators

  if (centrality==kV0M)      return fCentralityV0M;
  else if (centrality==kV0A) return fCentralityV0A; 
  else if (centrality==kV0C) return fCentralityV0C;
  else if (centrality==kCL1) return fCentralityCL1;
  else if (centrality==kZNA) return fCentralityZNA;
  else if (centrality==kZNC) return fCentralityZNC;
  else return -1.;
}

//_____________________________________________________________________________
UInt_t AliMuonsHFHeader::CollectPUMask(AliVEvent *event)
{
  // collect the mask for different combination of the parameters;
  // used to tag pile-up events (IsPileupFromSPD, no IsPileupFromSPDInMultBins at the moment) 

  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(event);
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);

  UInt_t collectMask = 0;
  Bool_t ISc1z1 = kFALSE, ISc1z2 = kFALSE, ISc1z3 = kFALSE, ISc1z4 = kFALSE,
         ISc2z1 = kFALSE, ISc2z2 = kFALSE, ISc2z3 = kFALSE, ISc2z4 = kFALSE,
         ISc3z1 = kFALSE, ISc3z2 = kFALSE, ISc3z3 = kFALSE, ISc3z4 = kFALSE,
         ISc4z1 = kFALSE, ISc4z2 = kFALSE, ISc4z3 = kFALSE, ISc4z4 = kFALSE;

  ISc1z1 = (aod) ? aod->IsPileupFromSPD(3,0.5,3.,2.,5.) : esd->IsPileupFromSPD(3,0.5,3.,2.,5.);
  ISc1z2 = (aod) ? aod->IsPileupFromSPD(3,0.6,3.,2.,5.) : esd->IsPileupFromSPD(3,0.6,3.,2.,5.);
  ISc1z3 = (aod) ? aod->IsPileupFromSPD(3,0.8,3.,2.,5.) : esd->IsPileupFromSPD(3,0.8,3.,2.,5.);
  ISc1z4 = (aod) ? aod->IsPileupFromSPD(3,0.9,3.,2.,5.) : esd->IsPileupFromSPD(3,0.9,3.,2.,5.);
  ISc2z1 = (aod) ? aod->IsPileupFromSPD(4,0.5,3.,2.,5.) : esd->IsPileupFromSPD(4,0.5,3.,2.,5.);
  ISc2z2 = (aod) ? aod->IsPileupFromSPD(4,0.6,3.,2.,5.) : esd->IsPileupFromSPD(4,0.6,3.,2.,5.);
  ISc2z3 = (aod) ? aod->IsPileupFromSPD(4,0.8,3.,2.,5.) : esd->IsPileupFromSPD(4,0.8,3.,2.,5.);
  ISc2z4 = (aod) ? aod->IsPileupFromSPD(4,0.9,3.,2.,5.) : esd->IsPileupFromSPD(4,0.9,3.,2.,5.);
  ISc3z1 = (aod) ? aod->IsPileupFromSPD(5,0.5,3.,2.,5.) : esd->IsPileupFromSPD(5,0.5,3.,2.,5.);
  ISc3z2 = (aod) ? aod->IsPileupFromSPD(5,0.6,3.,2.,5.) : esd->IsPileupFromSPD(5,0.6,3.,2.,5.);
  ISc3z3 = (aod) ? aod->IsPileupFromSPD(5,0.8,3.,2.,5.) : esd->IsPileupFromSPD(5,0.8,3.,2.,5.);
  ISc3z4 = (aod) ? aod->IsPileupFromSPD(5,0.9,3.,2.,5.) : esd->IsPileupFromSPD(5,0.9,3.,2.,5.);
  ISc4z1 = (aod) ? aod->IsPileupFromSPD(6,0.5,3.,2.,5.) : esd->IsPileupFromSPD(6,0.5,3.,2.,5.);
  ISc4z2 = (aod) ? aod->IsPileupFromSPD(6,0.6,3.,2.,5.) : esd->IsPileupFromSPD(6,0.6,3.,2.,5.);
  ISc4z3 = (aod) ? aod->IsPileupFromSPD(6,0.8,3.,2.,5.) : esd->IsPileupFromSPD(6,0.8,3.,2.,5.);
  ISc4z4 = (aod) ? aod->IsPileupFromSPD(6,0.9,3.,2.,5.) : esd->IsPileupFromSPD(6,0.9,3.,2.,5.);
  if (ISc1z1) collectMask |= kPUc1z1;
  if (ISc1z2) collectMask |= kPUc1z2;
  if (ISc1z3) collectMask |= kPUc1z3;
  if (ISc1z4) collectMask |= kPUc1z4;
  if (ISc2z1) collectMask |= kPUc2z1;
  if (ISc2z2) collectMask |= kPUc2z2;
  if (ISc2z3) collectMask |= kPUc2z3;
  if (ISc2z4) collectMask |= kPUc2z4;
  if (ISc3z1) collectMask |= kPUc3z1;
  if (ISc3z2) collectMask |= kPUc3z2;
  if (ISc3z3) collectMask |= kPUc3z3;
  if (ISc3z4) collectMask |= kPUc3z4;
  if (ISc4z1) collectMask |= kPUc4z1;
  if (ISc4z2) collectMask |= kPUc4z2;
  if (ISc4z3) collectMask |= kPUc4z3;
  if (ISc4z4) collectMask |= kPUc4z4;
 
  aod = 0; esd = 0;

  return collectMask;
}
