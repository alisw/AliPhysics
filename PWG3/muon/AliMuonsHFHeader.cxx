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
#include <TH2F.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>

#include "AliTriggerAnalysis.h"
#include "AliBackgroundSelection.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliMuonInfoStoreRD.h"
#include "AliMuonInfoStoreMC.h"
#include "AliDimuInfoStoreRD.h"
#include "AliDimuInfoStoreMC.h"
#include "AliMuonsHFHeader.h"

class TNamed;
class AliESDVertex;

ClassImp(AliMuonsHFHeader)

const TString AliMuonsHFHeader::fgkStdBranchName("MuEvsH");
Int_t         AliMuonsHFHeader::fgAnaMode         = 0;
Bool_t        AliMuonsHFHeader::fgIsMC            = kFALSE;
Bool_t        AliMuonsHFHeader::fgIsEventSelected = kFALSE;
Double_t      AliMuonsHFHeader::fgCuts[3] = { -999999., 999999., 999999.};

//_____________________________________________________________________________
AliMuonsHFHeader::AliMuonsHFHeader() :
TNamed(),
fTriggerMask(0),
fFiredTrigger(0),
fNFiredTrigger(0),
fIsPhysicsTriggered(kFALSE),
fIsPhysicsAccepted(kFALSE),
fEventType(0),
fUnrecoVertex(kFALSE),
fNContributors(0),
fUnrecoVtxSPD(kFALSE),
fNContributorsSPD(0),
fNTrackletsSPD(0),
fCentrality(0.)
{
  //
  // default constructor
  //
  for (Int_t i=3; i--;) fVtx[i] = 0.;
  for (Int_t i=3; i--;) fVtxSPD[i] = 0.;
}

//_____________________________________________________________________________
AliMuonsHFHeader::AliMuonsHFHeader(const AliMuonsHFHeader &src) :
TNamed(),
fTriggerMask(src.fTriggerMask),
fFiredTrigger(src.fFiredTrigger),
fNFiredTrigger(src.fNFiredTrigger),
fIsPhysicsTriggered(src.fIsPhysicsTriggered),
fIsPhysicsAccepted(src.fIsPhysicsAccepted),
fEventType(src.fEventType),
fUnrecoVertex(src.fUnrecoVertex),
fNContributors(src.fNContributors),
fUnrecoVtxSPD(src.fUnrecoVtxSPD),
fNContributorsSPD(src.fNContributorsSPD),
fNTrackletsSPD(src.fNTrackletsSPD),
fCentrality(src.fCentrality)
{
  //
  // copy constructor
  //
  for (Int_t i=3; i--;) fVtx[i]    = src.fVtx[i];
  for (Int_t i=3; i--;) fVtxSPD[i] = src.fVtxSPD[i];
}

//_____________________________________________________________________________
AliMuonsHFHeader& AliMuonsHFHeader::operator=(const AliMuonsHFHeader &src)
{
  //
  // assignment constructor
  //
  fTriggerMask        = src.fTriggerMask;
  fFiredTrigger       = src.fFiredTrigger;
  fNFiredTrigger      = src.fNFiredTrigger; 
  fIsPhysicsTriggered = src.fIsPhysicsTriggered;
  fIsPhysicsAccepted  = src.fIsPhysicsAccepted;
  fEventType          = src.fEventType;
  fUnrecoVertex       = src.fUnrecoVertex;
  fNContributors      = src.fNContributors;
  fUnrecoVtxSPD       = src.fUnrecoVtxSPD;
  fNContributorsSPD   = src.fNContributorsSPD;
  fNTrackletsSPD      = src.fNTrackletsSPD;
  fCentrality         = src.fCentrality;

  for (Int_t i=3; i--;) fVtx[i]      = src.fVtx[i];
  for (Int_t i=3; i--;) fVtxSPD[i]   = src.fVtxSPD[i];

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
void AliMuonsHFHeader::SetEvent(AliAODEvent *event)
{
  // extract event info from AOD event

  fTriggerMask = event->GetTriggerMask();

  AliAODVertex *vertex = event->GetPrimaryVertex(); 
  vertex->GetXYZ(fVtx);
  fNContributors = vertex->GetNContributors();
  fUnrecoVertex = (TMath::Abs(fVtx[0])<1e-6 && TMath::Abs(fVtx[1])<1e-6 &&
                   TMath::Abs(fVtx[2])<1e-6);
  this->SetTitle(vertex->GetTitle());

  this->SetFiredTrigger(event->GetFiredTriggerClasses());
  this->EventSelection();
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::SetEvent(AliESDEvent *event)
{
  // extract event info from ESD event
  // SPD vertex and Physics event selection are implimented

  fTriggerMask = event->GetTriggerMask();

  const AliESDVertex *vertex = event->GetPrimaryVertex(); 
  vertex->GetXYZ(fVtx);
  fNContributors = vertex->GetNContributors();
  fUnrecoVertex = (TMath::Abs(fVtx[0])<1e-6 && TMath::Abs(fVtx[1])<1e-6 &&
                   TMath::Abs(fVtx[2])<1e-6);
  this->SetTitle(vertex->GetTitle());

  const AliESDVertex *vtxSPD = event->GetPrimaryVertexSPD();
  vtxSPD->GetXYZ(fVtxSPD);
  fNContributorsSPD = vtxSPD->GetNContributors();
  fUnrecoVtxSPD = (TMath::Abs(fVtxSPD[0])<1e-6 && TMath::Abs(fVtxSPD[1])<1e-6 &&
                   TMath::Abs(fVtxSPD[2])<1e-6);
  fNTrackletsSPD = event->GetMultiplicity()->GetNumberOfTracklets();

  this->SetFiredTrigger(event->GetFiredTriggerClasses());
  this->PhysicsTriggerAna(event);
  this->EventSelection();
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::PhysicsTriggerAna(const AliESDEvent *esd)
{
  // ESD event trigger condition analysis
  // according to the method in $ALICE_ROOT/ANALYSIS/AliPhysicsSelection.cxx

  fEventType = esd->GetHeader()->GetEventType();
  fIsPhysicsTriggered = kFALSE;
  fIsPhysicsAccepted  = kFALSE;

  AliTriggerAnalysis *triggerAna = new AliTriggerAnalysis();
  triggerAna->SetAnalyzeMC(fgIsMC);
  triggerAna->SetSPDGFOThreshhold(1);

  Int_t  triggerHW  = triggerAna->SPDFiredChips(esd, 1);
  Bool_t isFiredv0A = triggerAna->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0A);
  Bool_t isFiredv0C = triggerAna->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0C);
  if (triggerHW==0 && !isFiredv0A && !isFiredv0C) {
    delete triggerAna;
    triggerAna = 0;
    return;
  }

  Bool_t triggerBG = (triggerAna->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0ABG) ||
                      triggerAna->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0CBG));
  Int_t  isFiredSPD = triggerAna->SPDFiredChips(esd, 0);
  Bool_t triggerFD = ((isFiredSPD>1) || (isFiredSPD>0 && (isFiredv0A || isFiredv0C)) || (isFiredv0A && isFiredv0C));
  if ((!triggerBG) && triggerFD) fIsPhysicsTriggered = kTRUE; 
  delete triggerAna;
  triggerAna = 0;

  if (fIsPhysicsTriggered) {
    AliBackgroundSelection *bkgId = new AliBackgroundSelection();
    fIsPhysicsAccepted = bkgId->IsSelected(const_cast<AliESDEvent*>(esd));
    delete bkgId; bkgId=0;
  }
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::SetFiredTrigger(TString trigger)
{
  // get info of fired trigger classes of event

  fFiredTrigger = trigger;
  TObjArray *tokens = trigger.Tokenize(" ");
  fNFiredTrigger = tokens->GetEntries();
  return;
}

//_____________________________________________________________________________
Bool_t AliMuonsHFHeader::IsTriggerFired(TString trigger)
{
  // check whether the trigger class "trigger" is fired in this event

  if (fNFiredTrigger<=0) return kFALSE;

  TObjArray *tokens = fFiredTrigger.Tokenize(" "); 
  for (Int_t i=fNFiredTrigger; i--;) {
    TString fired = ((TObjString*)tokens->At(i))->String(); 
    if (fired.CompareTo(trigger)==0) return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::EventSelection(TString triggerName)
{
  // select event according to the "triggerName" & event selection cuts

  fgIsEventSelected = kFALSE;
  if (!this->IsPhysicsAccepted()) return;
  if (!fgIsMC && !this->IsTriggerFired(triggerName)) return;
  this->EventSelection();
  return; 
}

//_____________________________________________________________________________
void AliMuonsHFHeader::EventSelection()
{
  // select event according to the event selection cuts

  fgIsEventSelected = kFALSE;
  if (this->NVtxContributorsSPD()<fgCuts[0]) return;
  if (TMath::Abs(this->VzSPD())>fgCuts[1]) return;
  if (this->VtSPD()>fgCuts[2]) return;
  fgIsEventSelected = kTRUE;
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistograms(TList *listEvent, TList *listMuon, TList *listDimu)
{
  // create output histos of muon analysis according to the analysis mode & MC flag

  this->CreateHistosEventH(listEvent);
  if (fgIsMC) {
    if (fgAnaMode!=2) this->CreateHistosMuonMC(listMuon);
    if (fgAnaMode!=1) this->CreateHistosDimuMC(listDimu);
  } else {
    if (fgAnaMode!=2) this->CreateHistosMuonRD(listMuon);
    if (fgAnaMode!=1) this->CreateHistosDimuRD(listDimu);
  }
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosEventH(TList *list)
{
  // create histograms at event level

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nHistos = 5;
  char    *name[nHistos] = { "hVztSPD", "hVzxSPD", "hVzySPD", "hVzNcontr", "hVtNcontr"  };
  Int_t  nbinsX[nHistos] = {      800 ,      800 ,      800 ,       800  ,        400   };
  Double_t xlow[nHistos] = {      -40.,      -40.,      -40.,       -40. ,          0.  };
  Double_t  xup[nHistos] = {       40.,       40.,       40.,        40. ,          4.  };
  Int_t  nbinsY[nHistos] = {      400 ,      600 ,      600 ,       202  ,        202   };
  Double_t ylow[nHistos] = {        0.,       -3.,       -3.,        -2.5,         -2.5 };
  Double_t  yup[nHistos] = {        4.,        3.,        3.,       199.5,        199.5 };

  TH2F *histo = 0;
  for (Int_t i=0; i<nHistos; i++) {
    histo = new TH2F(name[i], name[i], nbinsX[i], xlow[i], xup[i], nbinsY[i], ylow[i], yup[i]);
    histo->Sumw2(); list->AddAt(histo, i); histo = 0;
  }

  TH1::AddDirectory(oldStatus);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosMuonRD(TList *list)
{
  // create histograms for single muon

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  const Int_t nHistos = 8; 
  char    *name[nHistos] = {  "hP", "hPt", "hEta", "hDCA", "hTrg", "hCharge", "hUnfVtx" , "hEtaTimesDCA"};
  Int_t   nbins[nHistos] = { 1500 ,  300 ,   100 ,   500 ,    4  ,       3  ,       80  ,         10000 };
  Double_t xlow[nHistos] = {    0.,    0.,   -10.,     0.,   -0.5,      -1.5,      -40. ,             0.};
  Double_t  xup[nHistos] = {  150.,   30.,     0.,   500.,    3.5,       1.5,       40. ,          1000.};

  TH1F *histo = 0;
  for (Int_t i=0; i<nHistos; i++) {
    histo = new TH1F(name[i], name[i], nbins[i], xlow[i], xup[i]);
    histo->Sumw2(); list->AddAt(histo, i); histo = 0;
  }

  TH1::AddDirectory(oldStatus);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosDimuRD(TList *list)
{
  // create histograms for dimuon

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1F *histo = 0;
  const Int_t nHistos = 3;
  char    *name[nHistos] = {  "hP", "hPt", "hInvM" };
  Int_t   nbins[nHistos] = { 1500 ,  300 ,    300  };
  Double_t xlow[nHistos] = {    0.,    0.,      0. };
  Double_t  xup[nHistos] = {  150.,   30.,     30. };
  char *dimuName[3] = {"DimuNN", "DimuNP", "DimuPP"};
  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<nHistos; j++) {
      histo = new TH1F(Form("%s_%s",name[j],dimuName[i]), name[j], nbins[j], xlow[j], xup[j]);
      histo->Sumw2(); list->AddAt(histo,i*nHistos+j); histo = 0;
    }
  }

  TH1::AddDirectory(oldStatus);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosMuonMC(TList *list)
{
  // create histograms for single muon with MC info

  if (!list) list = new TList[AliMuonInfoStoreMC::NSources()];
  for (Int_t i=AliMuonInfoStoreMC::NSources(); i--;) this->CreateHistosMuonRD(&list[i]);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::CreateHistosDimuMC(TList *list)
{
  // create histograms for dimuon with MC info

  if (!list) list = new TList[AliDimuInfoStoreMC::NSources()];
  for (Int_t i=AliDimuInfoStoreMC::NSources(); i--;) this->CreateHistosDimuRD(&list[i]);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosEventH(TList *list)
{
  // fill histograms at event level according to event selection cuts

  if (!list) return;
  if (!AliMuonsHFHeader::IsEventSelected()) return;

  const Int_t nHistos = 5;
  Double_t vz = this->VzSPD();
  Double_t vt = this->VtSPD();
  Int_t    nc = this->NVtxContributorsSPD();
  Double_t distX[nHistos] = { vz,            vz,            vz, vz, vt };
  Double_t distY[nHistos] = { vt, this->VxSPD(), this->VySPD(), nc, nc };
  for (Int_t i=nHistos; i--;) ((TH2F*)list->At(i))->Fill(distX[i], distY[i]);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosMuonRD(TList *list, AliMuonInfoStoreRD* const muonStoreRD)
{
  // fill histograms for single muon according to event & muon track selection cuts

  if (!list) return;
  if (!AliMuonsHFHeader::IsEventSelected()) return;
  if (!muonStoreRD->MuonSelection()) return;

  const Int_t nHistos = 8;
  Double_t dist[nHistos] = {muonStoreRD->Momentum().Mag(),
                            muonStoreRD->Momentum().Pt(),
                            muonStoreRD->Momentum().Eta(),
                            muonStoreRD->DCA(),
                            muonStoreRD->MatchTrigger(),
                            muonStoreRD->Charge(),
                            this->VzSPD(),
                            muonStoreRD->DCA()*TMath::Exp(-2.*muonStoreRD->Momentum().Eta())/1000.};
  for (Int_t i=nHistos; i--;) ((TH1F*)list->At(i))->Fill(dist[i]);
  return; 
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosDimuRD(TList *list, AliDimuInfoStoreRD* const dimuStoreRD)
{
  // fill histograms for dimuon according to evnet & dimuon candidates selection cuts

  if (!list) return;
  if (!AliMuonsHFHeader::IsEventSelected()) return;
  if (!dimuStoreRD->DimuSelection()) return;

  Int_t theDimu = 0;
  if (dimuStoreRD->Charge()==0)     theDimu = 1;
  else if (dimuStoreRD->Charge()>0) theDimu = 2;

  const Int_t nHistos = 3;
  Double_t dist[nHistos] = {dimuStoreRD->Momentum().Mag(),
                            dimuStoreRD->Momentum().Pt(),
                            dimuStoreRD->InvM()};
  for (Int_t i=nHistos; i--;) ((TH1F*)list->At(theDimu*nHistos+i))->Fill(dist[i]);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosMuonMC(TList *list, AliMuonInfoStoreMC* const muonStoreMC)
{
  // fill histograms for single muon with MC info

  Int_t srcMuon = muonStoreMC->MuonSource(); 
  this->FillHistosMuonRD(&list[6], (AliMuonInfoStoreRD*)muonStoreMC);
  if (srcMuon<0) this->FillHistosMuonRD(&list[5], (AliMuonInfoStoreRD*)muonStoreMC);
  else this->FillHistosMuonRD(&list[srcMuon], (AliMuonInfoStoreRD*)muonStoreMC);
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::FillHistosDimuMC(TList *list, AliDimuInfoStoreMC* const dimuStoreMC)
{
  // fill histograms for dimuon with MC info

  Int_t srcDimu = dimuStoreMC->DimuSource();
  this->FillHistosDimuRD(&list[6], (AliDimuInfoStoreRD*)dimuStoreMC);
  this->FillHistosDimuRD(&list[srcDimu], (AliDimuInfoStoreRD*)dimuStoreMC);
  return;
}
