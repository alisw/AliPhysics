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

///////////////////////////////////////////////////////////////////////////
//
// class used to extract ,store info and fill histograms at event level
//
// Author: H-S. Zhu, hongsheng.zhu@cern.ch
//                   hszhu@iopp.ccnu.edu.cn
///////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TList.h>
#include <TClonesArray.h>

#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVCaloCells.h"
#include "AliPHOSGeoUtils.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"

#include "AliCaloClusterInfo.h"
#include "AliPHOSpPbPi0Header.h"

class TNamed;

ClassImp(AliPHOSpPbPi0Header)

Bool_t   AliPHOSpPbPi0Header::fgIsMC           = kFALSE;
Bool_t   AliPHOSpPbPi0Header::fgIspARun        = kFALSE;
Bool_t   AliPHOSpPbPi0Header::fgUseFiducialCut = kFALSE;
Int_t    AliPHOSpPbPi0Header::fgNCent          = 10;
Double_t AliPHOSpPbPi0Header::fgCuts[4]        = { 1., 10., 0., 100. };

//_____________________________________________________________________________
AliPHOSpPbPi0Header::AliPHOSpPbPi0Header() :
TNamed(),
fFiredTriggerClass(),
fSelMask(0),
fVtxContrsN(0),
fIsVertexOK(kFALSE),
fIsPileupSPD(kFALSE),
fCentrality(0.)
{
  //
  // default constructor
  //
  for (Int_t i=3; i--;) fVtx[i] = 0.;
}

//_____________________________________________________________________________
AliPHOSpPbPi0Header::AliPHOSpPbPi0Header(const AliPHOSpPbPi0Header &src) :
TNamed(),
fFiredTriggerClass(src.fFiredTriggerClass),
fSelMask(src.fSelMask),
fVtxContrsN(src.fVtxContrsN),
fIsVertexOK(src.fIsVertexOK),
fIsPileupSPD(src.fIsPileupSPD),
fCentrality(src.fCentrality)
{
  //
  // copy constructor
  //
  for (Int_t i=3; i--;) fVtx[i] = src.fVtx[i];
}

//_____________________________________________________________________________
AliPHOSpPbPi0Header& AliPHOSpPbPi0Header::operator=(const AliPHOSpPbPi0Header &src)
{
  //
  // assignment constructor
  //

  if(&src==this) return *this;

  fFiredTriggerClass            = src.fFiredTriggerClass;
  fSelMask                      = src.fSelMask;
  fVtxContrsN                   = src.fVtxContrsN;
  fIsVertexOK                   = src.fIsVertexOK;
  fIsPileupSPD                  = src.fIsPileupSPD;
  fCentrality                   = src.fCentrality;
  for (Int_t i=3; i--;) fVtx[i] = src.fVtx[i];

  return *this;
}

//_____________________________________________________________________________
AliPHOSpPbPi0Header::~AliPHOSpPbPi0Header()
{
  //
  // default destructor
  //
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::SetEventInfo(AliInputEventHandler* const handler)
{
  // fill info at event level

  AliVEvent *event = handler->GetEvent();
  AliAODEvent *aod = 0x0;   aod = dynamic_cast<AliAODEvent*>(event);
  AliESDEvent *esd = 0x0;   esd = dynamic_cast<AliESDEvent*>(event);

  fIsVertexOK   = CheckEventVertex(aod, esd);
  fSelMask      = handler->IsEventSelected();
  fIsPileupSPD  = (aod && !aod->GetTracklets()) ? event->IsPileupFromSPD(3,0.8,3.,2.,5.) :
                                                  event->IsPileupFromSPDInMultBins(); //TODO not sure!!!

  AliCentrality *cent = event->GetCentrality();
  if (cent) fCentrality = cent->GetCentralityPercentile("V0M");

//this->SetTitle(vertex->GetTitle());
  return;
}

//_____________________________________________________________________________
Bool_t AliPHOSpPbPi0Header::IsSelected()
{
  // select event according to the event selection cuts
  if (fVtxContrsN<(Int_t)fgCuts[0])                    return kFALSE; // num. of vtx contributors cut
  if (TMath::Abs(fVtx[2])>fgCuts[1])                   return kFALSE; // Vtxz cut
  if (fCentrality<fgCuts[2] || fCentrality>fgCuts[3])  return kFALSE; // centrality selection
  if (fgIspARun && !fIsVertexOK)                       return kFALSE; // pA vertex cut

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPHOSpPbPi0Header::CheckEventVertex(AliAODEvent* const aod, AliESDEvent* const esd)
{
  // check event vertex

  Bool_t   isAOD = kFALSE;
  if (aod) isAOD = kTRUE;
  if (esd) isAOD = kFALSE;

  // set event basic info
  fFiredTriggerClass       = isAOD ? aod->GetFiredTriggerClasses() : esd->GetFiredTriggerClasses();
  const AliVVertex* trkVtx = isAOD ? dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertex()) :
                                     dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertex()) ;
  trkVtx->GetXYZ(fVtx);   fVtxContrsN  = trkVtx->GetNContributors();

  if (!trkVtx || fVtxContrsN<1)                                                           return kFALSE;
  if (!fgIspARun)                                                                         return kTRUE;  // Following cuts are for pA vertex
  if (!((TString)trkVtx->GetTitle()).Contains("VertexerTracks"))                          return kFALSE;

  const AliVVertex* spdVtx = isAOD ? dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertexSPD()) :
                                     dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertexSPD()) ;
  Double_t cov[6]={0.};  spdVtx->GetCovarianceMatrix(cov);

  if (spdVtx->GetNContributors()<1)                                                       return kFALSE;
  if (((TString)spdVtx->GetTitle()).Contains("vertexer:Z") && (TMath::Sqrt(cov[5])>0.25)) return kFALSE; // Double_t zRes = TMath::Sqrt(cov[5]);
  if (TMath::Abs(spdVtx->GetZ() - fVtx[2])>0.5)                                           return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistograms(TList *listQA, TList *listRD, TList *listMC)
{
  // create output histos of pi0 analysis according to the MC flag

  this->CreateHistosEvent(listQA);
  this->CreateHistosCaloCellsQA(listQA);
  this->CreateHistosCaloCluster(listQA);
  this->CreateHistosPi0(listRD);
  this->CreateHistosMixPi0(listRD);
  if (fgIsMC) this->CreateHistosMC(listMC);

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosEvent(TList *list)
{
  // create histograms at event level

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nhs    = 4;
  TString tName[nhs] = { "VtxNcontr", "Centrality", "Vz", "Pileup" };
  Int_t   nbins[nhs] = {       202  ,         110 , 500 ,     500  };
  Double_t xlow[nhs] = {        -2.5,         -10., -25.,     -25. };
  Double_t  xup[nhs] = {       199.5,         100.,  25.,      25. };

  TString hName;
  TH1D   *histo = 0;
  for (Int_t i=0; i<nhs; i++) {
    hName = Form("hEvent_%s", tName[i].Data());
    histo = new TH1D(hName.Data(), hName.Data(), nbins[i], xlow[i], xup[i]);
    histo->Sumw2(); list->Add(histo); histo = 0;

    hName = Form("hEventSel_%s", tName[i].Data());
    histo = new TH1D(hName.Data(), hName.Data(), nbins[i], xlow[i], xup[i]);
    histo->Sumw2(); list->Add(histo); histo = 0;
  }
  TH1::AddDirectory(oldStatus);

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosCaloCellsQA(TList *list)
{
  // create QA histograms for CaloCells

  if (!list) list = new TList();
  list->SetOwner();
  Bool_t oldStatus1 = TH1::AddDirectoryStatus(); TH1::AddDirectory(kFALSE);
  Bool_t oldStatus2 = TH2::AddDirectoryStatus(); TH2::AddDirectory(kFALSE);

  TH1D *histo1       = 0;
  TH2D *histo2       = 0;
  const Int_t nhs    = 3;
  const Int_t nMod   = 3;
  TString hName = "hCaloCellsQA_E";
  histo1 = new TH1D(hName, hName, 300, 0., 30.);
  histo1->Sumw2(); list->Add(histo1); histo1 = 0;
  for (Int_t i=0; i<nMod; i++) {
    hName = Form("hCaloCellsQA_E_Mod%d", i+1);
    histo1 = new TH1D(hName.Data(), hName.Data(), 300, 0., 30.);
    histo1->Sumw2(); list->Add(histo1); histo1 = 0;
  }

  TString tName[nhs] = { "CentNCells", "Nxz", "Exz" };
  Int_t  nbinsx[nhs] = {          10 ,   64 ,  64   };
  Double_t xlow[nhs] = {           0.,   0.5,   0.5 };
  Double_t  xup[nhs] = {         100.,  64.5,  64.5 };
  Int_t  nbinsy[nhs] = {        1000 ,   56 ,  56   };
  Double_t ylow[nhs] = {           0.,   0.5,   0.5 };
  Double_t  yup[nhs] = {        1000.,  56.5,  56.5 };

  for (Int_t i=0; i<nhs; i++) {
    if (i == 0) {
      hName = Form("hCaloCellsQA_%s", tName[i].Data());
      histo2 = new TH2D(hName.Data(), hName.Data(), nbinsx[i], xlow[i], xup[i], nbinsy[i], ylow[i], yup[i]);
      histo2->Sumw2(); list->Add(histo2); histo2 = 0;
    } else {
      for (Int_t j=0; j<nMod; j++) {
        hName = Form("hCaloCellsQA_%s_Mod%d", tName[i].Data(), j+1);
        histo2 = new TH2D(hName.Data(), hName.Data(), nbinsx[i], xlow[i], xup[i], nbinsy[i], ylow[i], yup[i]);
        histo2->Sumw2(); list->Add(histo2); histo2 = 0;
      }
    }
  }

  TH1::AddDirectory(oldStatus1);
  TH2::AddDirectory(oldStatus2);
  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosCaloCluster(TList *list)
{
  // create histograms for CaloCluster

  if (!list) list = new TList();
  list->SetOwner();

  Bool_t oldStatus1 = TH1::AddDirectoryStatus();   TH1::AddDirectory(kFALSE);
  Bool_t oldStatus2 = TH2::AddDirectoryStatus();   TH2::AddDirectory(kFALSE);

  const Int_t nMods = 3;
  const Int_t nTRUs = 8;
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs };   // PID
  TString namePID[kPIDs] =  { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"           };
                         // { kPtClu, kEtaClu, kPhiClu, kM02Clu, kM20Clu, kTOFClu, kNCellsClu, kNClustersClu, kVarsClu };   // clusters
  TString tName[kVarsClu] = {   "Pt",   "Eta",   "Phi",   "M02",   "M20",   "TOF",   "NCells",   "NClusters"           };
  Int_t    bins[kVarsClu] = {   500 ,  300   ,   120  ,    200 ,    200 ,    400 ,      1000 ,          20             };
  Double_t xMin[kVarsClu] = {     0.,   -0.15,     4.4,      0.,      0.,   -2e-6,         0.,           0.            };
  Double_t xMax[kVarsClu] = {    50.,    0.15,     5.6,      2.,      2.,    2e-6,      1000.,          20.            };

  TH1I *hncells = 0;
  TH1D *histo1  = 0;
  TH2D *histo2  = 0;
  TString hName;

  hName   = Form("hCaloCluster_%s", tName[kNCellsClu].Data());
  hncells = new TH1I(hName.Data(), hName.Data(), bins[kNCellsClu], (Int_t)xMin[kNCellsClu], (Int_t)xMax[kNCellsClu]);
  hncells->Sumw2(); list->Add(hncells); hncells = 0;

  hName   = Form("hCaloCluster_%s", tName[kNClustersClu].Data());
  hncells = new TH1I(hName.Data(), hName.Data(), bins[kNClustersClu], (Int_t)xMin[kNClustersClu], (Int_t)xMax[kNClustersClu]);
  hncells->Sumw2(); list->Add(hncells); hncells = 0;

  for (Int_t iPID=0; iPID<kPIDs; iPID++) {
    hName  = Form("hCaloCluster_%s_%s", tName[kPtClu].Data(), namePID[iPID].Data());
    histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtClu], xMin[kPtClu], xMax[kPtClu]);
    histo1->Sumw2(); list->Add(histo1); histo1 = 0;

    for (Int_t icent=0; icent<fgNCent; icent++) {
      hName  = Form("hCaloCluster_%s_%s_cent%d", tName[kPtClu].Data(), namePID[iPID].Data(), icent);
      histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtClu], xMin[kPtClu], xMax[kPtClu]);
      histo1->Sumw2(); list->Add(histo1); histo1 = 0;
    }

    hName = Form("hCaloCluster_%s%s_%s", tName[kEtaClu].Data(), tName[kPhiClu].Data(), namePID[iPID].Data());
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[kEtaClu], xMin[kEtaClu], xMax[kEtaClu], bins[kPhiClu], xMin[kPhiClu], xMax[kPhiClu]);
    histo2->Sumw2(); list->Add(histo2); histo2 = 0;
    for (Int_t i=0; i<kVarsClu-3; i++) {
      hName = Form("hCaloCluster_%s%s_%s", tName[0].Data(), tName[i+1].Data(), namePID[iPID].Data());
      histo2 = new TH2D(hName.Data(), hName.Data(), bins[0], xMin[0], xMax[0], bins[i+1], xMin[i+1], xMax[i+1]);
      histo2->Sumw2(); list->Add(histo2); histo2 = 0;
    }
  }

  for (Int_t iMod=0; iMod<nMods; iMod++) {
    hName   = Form("hCaloCluster_%s_Mod%d", tName[kNCellsClu].Data(), iMod+1);
    hncells = new TH1I(hName.Data(), hName.Data(), bins[kNCellsClu], (Int_t)xMin[kNCellsClu], (Int_t)xMax[kNCellsClu]);
    hncells->Sumw2(); list->Add(hncells); hncells = 0;

    hName   = Form("hCaloCluster_%s_Mod%d", tName[kNClustersClu].Data(), iMod+1);
    hncells = new TH1I(hName.Data(), hName.Data(), bins[kNClustersClu], (Int_t)xMin[kNClustersClu], (Int_t)xMax[kNClustersClu]);
    hncells->Sumw2(); list->Add(hncells); hncells = 0;

    for (Int_t iTRU=0; iTRU<nTRUs; iTRU++) {
      hName  = Form("hCaloCluster_%s_Mod%d_TRU%d", tName[kPtClu].Data(), iMod+1, iTRU+1);
      histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtClu], xMin[kPtClu], xMax[kPtClu]);
      histo1->Sumw2(); list->Add(histo1); histo1 = 0;
    }

    for (Int_t iPID=0; iPID<kPIDs; iPID++) {
      hName  = Form("hCaloCluster_%s_Mod%d_%s", tName[kPtClu].Data(), iMod+1, namePID[iPID].Data());
      histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtClu], xMin[kPtClu], xMax[kPtClu]);
      histo1->Sumw2(); list->Add(histo1); histo1 = 0;
    }
  }

  hName = Form("hCaloCluster_%s%s", tName[kEtaClu].Data(), tName[kPhiClu].Data());
  histo2 = new TH2D(hName.Data(), hName.Data(), bins[kEtaClu], xMin[kEtaClu], xMax[kEtaClu], bins[kPhiClu], xMin[kPhiClu], xMax[kPhiClu]);
  histo2->Sumw2(); list->Add(histo2); histo2 = 0;
  for (Int_t i=0; i<kVarsClu-3; i++) {
    hName = Form("hCaloCluster_%s%s", tName[0].Data(), tName[i+1].Data());
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[0], xMin[0], xMax[0], bins[i+1], xMin[i+1], xMax[i+1]);
    histo2->Sumw2(); list->Add(histo2); histo2 = 0;
  }

  TH1::AddDirectory(oldStatus1);
  TH2::AddDirectory(oldStatus2);
  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosPi0(TList *list)
{
  // create histograms for Pi0

  if (!list) list = new TList();
  list->SetOwner();

  Bool_t oldStatus = TH2::AddDirectoryStatus();   TH2::AddDirectory(kFALSE);

  const Int_t nComb = 2;
  TString srcMod[nComb]  =  { "Mod11", "Mod33" }; 
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs      };  // PID
  TString namePID[kPIDs] =  { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"                };
                         // { kPtPi0, kEtaPi0, kPhiPi0, kAsyPi0,   kAnglePi0, kInvMassPi0, kVarsPi0 };  // pi0
  TString tName[kVarsPi0] = {   "Pt",   "Eta",   "Phi",   "Asy",     "Angle",   "InvMass"           };
  Int_t    bins[kVarsPi0] = {   500 ,  300   ,   120  ,    100 ,        180 ,       1000            };
  Double_t xMin[kVarsPi0] = {     0.,   -0.15,     4.4,      0.,          0.,          0.           };
  Double_t xMax[kVarsPi0] = {    50.,    0.15,     5.6,      1., TMath::Pi(),          1.           };

  TH2D   *histo = 0;
  TString hName;

  hName  = Form("hPi0_%s%s", tName[kEtaPi0].Data(), tName[kPhiPi0].Data());
  histo = new TH2D(hName.Data(), hName.Data(), bins[kEtaPi0], xMin[kEtaPi0], xMax[kEtaPi0], bins[kPhiPi0], xMin[kPhiPi0], xMax[kPhiPi0]);
  histo->Sumw2(); list->Add(histo); histo = 0;
  for (Int_t i=0; i<kVarsPi0-2; i++) {
    hName  = Form("hPi0_%s%s", tName[0].Data(), tName[i+1].Data());
    histo = new TH2D(hName.Data(), hName.Data(), bins[0], xMin[0], xMax[0], bins[i+1], xMin[i+1], xMax[i+1]);
    histo->Sumw2(); list->Add(histo); histo = 0;
  }

  for (Int_t iPID=0; iPID<kPIDs; iPID++) {
    hName  = Form("hPi0_%s%s_%s", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(), namePID[iPID].Data());
    histo = new TH2D(hName.Data(), hName.Data(), bins[kPtPi0], xMin[kPtPi0], xMax[kPtPi0], bins[kInvMassPi0], xMin[kInvMassPi0], xMax[kInvMassPi0]);
    histo->Sumw2(); list->Add(histo); histo = 0;

    hName  = Form("hPi0_%s%s_%s", tName[kEtaPi0].Data(), tName[kPhiPi0].Data(), namePID[iPID].Data());
    histo = new TH2D(hName.Data(), hName.Data(), bins[kEtaPi0], xMin[kEtaPi0], xMax[kEtaPi0], bins[kPhiPi0], xMin[kPhiPi0], xMax[kPhiPi0]);
    histo->Sumw2(); list->Add(histo); histo = 0;
    for (Int_t icent=0; icent<fgNCent; icent++) {
      hName  = Form("hPi0_%s%s_%s_cent%d", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(), namePID[iPID].Data(), icent);
      histo = new TH2D(hName.Data(), hName.Data(), bins[kPtPi0], xMin[kPtPi0], xMax[kPtPi0], bins[kInvMassPi0], xMin[kInvMassPi0], xMax[kInvMassPi0]);
      histo->Sumw2(); list->Add(histo); histo = 0;
    }

    for (Int_t i=0; i<kVarsPi0-2; i++) {
      hName  = Form("hPi0_%s%s_%s", tName[0].Data(), tName[i+1].Data(), namePID[iPID].Data());
      histo = new TH2D(hName.Data(), hName.Data(), bins[0], xMin[0], xMax[0], bins[i+1], xMin[i+1], xMax[i+1]);
      histo->Sumw2(); list->Add(histo); histo = 0;
    }
    for (Int_t iComb=0; iComb<nComb; iComb++) {
      hName  = Form("hPi0_%s%s_%s_%s", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(), namePID[iPID].Data(), srcMod[iComb].Data());
      histo = new TH2D(hName.Data(), hName.Data(), bins[kPtPi0], xMin[kPtPi0], xMax[kPtPi0], bins[kInvMassPi0], xMin[kInvMassPi0], xMax[kInvMassPi0]);
      histo->Sumw2(); list->Add(histo); histo = 0;
    }
  }

  TH2::AddDirectory(oldStatus);
  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosMixPi0(TList *list)
{
  // create Histograms for Mixed Pi0

  if (!list) list = new TList();
  list->SetOwner();

  Bool_t oldStatus = TH2::AddDirectoryStatus(); TH2::AddDirectory(kFALSE);
  const Int_t nComb = 2;
  TString  srcMod[nComb] =  { "Mod11", "Mod33" };
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs       };  // PID
  TString namePID[kPIDs] =  { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"                 };
                             //{ kPtMixPi0, kEtaMixPi0, kPhiMixPi0, kInvMassMixPi0, kVarsMixPi0      };  // Mix Pi0
  TString tName[kVarsMixPi0] = {      "Pt",      "Eta",      "Phi",      "InvMass"                   };
  Int_t    bins[kVarsMixPi0] = {      500 ,     300   ,      120  ,          1000                    };
  Double_t xMin[kVarsMixPi0] = {        0.,      -0.15,        4.4,             0.                   };
  Double_t xMax[kVarsMixPi0] = {       50.,       0.15,        5.6,             1.                   };

  TH2D   *histo = 0;
  TString hName;

  hName  = Form("hMixPi0_%s%s", tName[kEtaMixPi0].Data(), tName[kPhiMixPi0].Data());
  histo = new TH2D(hName.Data(), hName.Data(), bins[kEtaMixPi0], xMin[kEtaMixPi0], xMax[kEtaMixPi0],
                                               bins[kPhiMixPi0], xMin[kPhiMixPi0], xMax[kPhiMixPi0]);
  histo->Sumw2(); list->Add(histo); histo = 0;
  
  for (Int_t iPID=0; iPID<kPIDs; iPID++) {
    hName  = Form("hMixPi0_%s%s_%s", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(), namePID[iPID].Data());
    histo = new TH2D(hName.Data(), hName.Data(), bins[kPtMixPi0],      xMin[kPtMixPi0],      xMax[kPtMixPi0],
                                                 bins[kInvMassMixPi0], xMin[kInvMassMixPi0], xMax[kInvMassMixPi0]);
    histo->Sumw2(); list->Add(histo); histo = 0;
    hName  = Form("hMixPi0_%s%s_%s", tName[kEtaMixPi0].Data(), tName[kPhiMixPi0].Data(), namePID[iPID].Data());
    histo = new TH2D(hName.Data(), hName.Data(), bins[kEtaMixPi0], xMin[kEtaMixPi0], xMax[kEtaMixPi0], 
                                                 bins[kPhiMixPi0], xMin[kPhiMixPi0], xMax[kPhiMixPi0]);
    histo->Sumw2(); list->Add(histo); histo = 0;
    for (Int_t icent=0; icent<fgNCent; icent++) {
      hName  = Form("hMixPi0_%s%s_%s_cent%d", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(), namePID[iPID].Data(), icent);
      histo = new TH2D(hName.Data(), hName.Data(), bins[kPtMixPi0],      xMin[kPtMixPi0],      xMax[kPtMixPi0],
                                                   bins[kInvMassMixPi0], xMin[kInvMassMixPi0], xMax[kInvMassMixPi0]);
      histo->Sumw2(); list->Add(histo); histo = 0;
    } 
    for (Int_t iComb=0; iComb<nComb; iComb++) {
      hName  = Form("hMixPi0_%s%s_%s_%s", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(), namePID[iPID].Data(), srcMod[iComb].Data());
      histo = new TH2D(hName.Data(), hName.Data(), bins[kPtMixPi0],      xMin[kPtMixPi0],      xMax[kPtMixPi0],
                                                   bins[kInvMassMixPi0], xMin[kInvMassMixPi0], xMax[kInvMassMixPi0]);
      histo->Sumw2(); list->Add(histo); histo = 0;
    }
  }

  TH2::AddDirectory(oldStatus);
  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosMC(TList *list)
{
  // create Histograms for MC
  if (!list) list = new TList();
  list->SetOwner();

  Bool_t oldStatus1 = TH2::AddDirectoryStatus(); TH1::AddDirectory(kFALSE);
  Bool_t oldStatus2 = TH2::AddDirectoryStatus(); TH2::AddDirectory(kFALSE);
  const Int_t method = 2;
  const char *mName[method] = {"IsPrimary", "RCut"};
  const Int_t cut    = 3;
  const char *cName[cut]    = {"WideY", "NarrY", "Acc"};
  const Int_t type   = 3;
  const char *pName[type]   = {"Pi0", "Eta", "Gamma"};
                       //  { kVertexMC, kPtMC, kRapidityMC,          kPhiMC, kWeightMC, kVarsMC };   // MC
  TString tName[kVarsMC] = {  "Vertex",  "Pt",  "Rapidity",           "Phi",  "Weight"          };
  Int_t    bins[kVarsMC] = {     1000 ,  500 ,        200 ,            360 ,     100            };
  Double_t xMin[kVarsMC] = {        0.,    0.,         -1.,              0.,       0.           };
  Double_t xMax[kVarsMC] = {      500.,   50.,          1.,  TMath::TwoPi(),      10.           };

  TH1D   *histo1 = 0;
  TH2D   *histo2 = 0;
  TString hName;

  for (Int_t iType=0; iType<type; iType++) {
    hName  = Form("hMC%s_%s%s", pName[iType], tName[kPtMC].Data(), tName[kVertexMC].Data());
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC], bins[kVertexMC], xMin[kVertexMC], xMax[kVertexMC]);
    histo2->Sumw2(); list->Add(histo2); histo2 = 0;
    hName  = Form("hMC%s_%s%s", pName[iType], tName[kPtMC].Data(), tName[kRapidityMC].Data());
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC], bins[kRapidityMC], xMin[kRapidityMC], xMax[kRapidityMC]);
    histo2->Sumw2(); list->Add(histo2); histo2 = 0;
    hName  = Form("hMC%s_%s%s_%s", pName[iType], tName[kPtMC].Data(), tName[kVertexMC].Data(), mName[0]);
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC], bins[kVertexMC], xMin[kVertexMC], xMax[kVertexMC]);
    histo2->Sumw2(); list->Add(histo2); histo2 = 0;
    for (Int_t iMethod=0; iMethod<method; iMethod++) {
      for (Int_t iCent=0; iCent<fgNCent; iCent++) {
        for (Int_t i=0; i<kVarsMC-2; i++) {
          hName  = Form("hMC%s_%s_cent%d_%s", pName[iType], tName[i+1].Data(), iCent, mName[iMethod]);
          histo1 = new TH1D(hName.Data(), hName.Data(), bins[i+1], xMin[i+1], xMax[i+1] );
          histo1->Sumw2(); list->Add(histo1); histo1 = 0;

        }
        for (Int_t iCut=0; iCut<cut; iCut++) {
          hName  = Form("hMC%s_%s_cent%d_%s_%s", pName[iType], tName[kPtMC].Data(), iCent, cName[iCut], mName[iMethod]);
          histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC]);
          histo1->Sumw2(); list->Add(histo1); histo1 = 0;
        }
      }
    }
  }


  TH1::AddDirectory(oldStatus1);
  TH2::AddDirectory(oldStatus2);
  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosEvent(TList *list)
{
  // fill histograms at event level according to event selection cuts

  if (!list) return;

  const Int_t nhs      = 3;
  TString tName[nhs+1] = { "VtxNcontr", "Centrality",    "Vz", "Pileup" };
  Double_t dist[nhs]   = { fVtxContrsN,  fCentrality, fVtx[2]           };
  for (Int_t i=nhs; i--;) {
    ((TH1D*)list->FindObject(Form("hEvent_%s",tName[i].Data())))->Fill(dist[i]);
    if (this->IsSelected()) ((TH1D*)list->FindObject(Form("hEventSel_%s",tName[i].Data())))->Fill(dist[i]);
  }
  if (fIsPileupSPD) {
    ((TH1D*)list->FindObject(Form("hEvent_%s",tName[nhs].Data())))->Fill(dist[nhs-1]);
    if (this->IsSelected())  ((TH1D*)list->FindObject(Form("hEventSel_%s",tName[nhs].Data())))->Fill(dist[nhs-1]);
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosCaloCellsQA(TList *list, AliVCaloCells* const cells, AliPHOSGeoUtils* const phosGeo)
{
  // fill QA histograms for calocells according to evnet selection cuts

  if (!list) return;

  const Int_t nhs    = 3;
  TString tName[nhs] = { "CentNCells", "Nxz", "Exz" };

  ((TH2D*)list->FindObject(Form("hCaloCellsQA_%s", tName[0].Data())))->Fill(fCentrality, cells->GetNumberOfCells());
  for (Int_t iCell=0; iCell<cells->GetNumberOfCells(); iCell++) {
    Int_t relID[4] = {0,0,0,0};
    phosGeo->AbsToRelNumbering(cells->GetCellNumber(iCell),relID); // Int_t mod   = relID[0]; Int_t cellX = relID[2]; Int_t cellZ = relID[3];
    Double_t wight[2] = { 1., cells->GetAmplitude(iCell)};
    ((TH1D*)list->FindObject(Form("hCaloCellsQA_%s", "E")))->Fill(wight[1]);
    ((TH1D*)list->FindObject(Form("hCaloCellsQA_%s_Mod%d", "E", relID[0])))->Fill(wight[1]);
    for (Int_t i=1; i<nhs; i++)
      ((TH2D*)list->FindObject(Form("hCaloCellsQA_%s_Mod%d", tName[i].Data(), relID[0])))->Fill(relID[2], relID[3], wight[i-1]);
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosCaloCluster(TList *list, TClonesArray* const caloClArr, Int_t cent)
{
  // fill histograms for calocluster according to evnet && calocluster candidates selection cuts

  if (!list)  return;
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs       };  // PID
  TString namePID[kPIDs] =  { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"                 };
  UInt_t  PIDBIT[kPIDs]   = { 0x0,  BIT(0), BIT(1), BIT(0)|BIT(1), BIT(2), BIT(3), BIT(2)|BIT(3)     };
                         // { kPtClu, kEtaClu, kPhiClu, kM02Clu, kM20Clu, kTOFClu, kNCellsClu, kNClustersClu, kVarsClu };   // clusters
  TString tName[kVarsClu] = {   "Pt",   "Eta",   "Phi",   "M02",   "M20",   "TOF",   "NCells",   "NClusters"           };

  Int_t entries = caloClArr->GetEntriesFast();
  Int_t iMod = 0, iTRU = 0, nclusters[3] = { 0, 0, 0 };
  UInt_t pidBit = 0;
  AliCaloClusterInfo *caloCluster = 0;
  TLorentzVector gamma;
  for(Int_t i=0; i<entries; i++) { // Loop calo cluster
    caloCluster = (AliCaloClusterInfo*)caloClArr->At(i);
    iMod   = caloCluster->GetModule();   nclusters[iMod-1]++;
    iTRU   = caloCluster->GetTRUNumber();
    pidBit = caloCluster->GetPIDBit();
    gamma  = caloCluster->LorentzVector();
    Double_t vars[kVarsClu-1] = { gamma.E(),
                                  gamma.Eta(),
                                  (gamma.Phi()+TMath::TwoPi()),
                                  caloCluster->GetM02(),
                                  caloCluster->GetM20(),
                                  caloCluster->GetTOF(),
                                  (Double_t)caloCluster->GetNCells() };
 
    ((TH1I*)list->FindObject(Form("hCaloCluster_%s", tName[kNCellsClu].Data())))->Fill((Int_t)vars[kNCellsClu]);
    ((TH1I*)list->FindObject(Form("hCaloCluster_%s_Mod%d", tName[kNCellsClu].Data(), iMod)))->Fill((Int_t)vars[kNCellsClu]);
    ((TH1D*)list->FindObject(Form("hCaloCluster_%s_Mod%d_TRU%d", tName[kPtClu].Data(), iMod, iTRU)))->Fill(vars[kPtClu]);

    for (Int_t iPID=0; iPID<kPIDs; iPID++) {
      if ((pidBit&PIDBIT[iPID])==PIDBIT[iPID]) {
        ((TH1D*)list->FindObject(Form("hCaloCluster_%s_%s", tName[kPtClu].Data(), namePID[iPID].Data())))->Fill(vars[kPtClu]);
        ((TH1D*)list->FindObject(Form("hCaloCluster_%s_Mod%d_%s", tName[kPtClu].Data(), iMod, namePID[iPID].Data())))->Fill(vars[kPtClu]);
        ((TH1D*)list->FindObject(Form("hCaloCluster_%s_%s_cent%d", tName[kPtClu].Data(), namePID[iPID].Data(), cent)))->Fill(vars[kPtClu]);
        ((TH2D*)list->FindObject(Form("hCaloCluster_%s%s_%s", tName[kEtaClu].Data(), tName[kPhiClu].Data(),
                                                              namePID[iPID].Data())))->Fill(vars[kEtaClu], vars[kPhiClu]);
        for (Int_t j=0; j<kVarsClu-3; j++)
          ((TH2D*)list->FindObject(Form("hCaloCluster_%s%s_%s", tName[0].Data(), tName[j+1].Data(), namePID[iPID].Data())))->Fill(vars[0], vars[j+1]);
      }
    }

    ((TH2D*)list->FindObject(Form("hCaloCluster_%s%s", tName[kEtaClu].Data(), tName[kPhiClu].Data())))->Fill(vars[kEtaClu], vars[kPhiClu]);
    for (Int_t j=0; j<kVarsClu-3; j++) {
      ((TH2D*)list->FindObject(Form("hCaloCluster_%s%s", tName[0].Data(), tName[j+1].Data())))->Fill(vars[0], vars[j+1]);
    }
  } // End loop calo cluster

  ((TH1I*)list->FindObject(Form("hCaloCluster_%s", tName[kNClustersClu].Data())))->Fill(entries);
  for (Int_t jMod=1; jMod<4; jMod++) {
    if (jMod==2) continue;
    ((TH1I*)list->FindObject(Form("hCaloCluster_%s_Mod%d", tName[kNClustersClu].Data(), jMod)))->Fill(nclusters[jMod-1]);
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosPi0(TList *list, TClonesArray* const caloClArr, Int_t cent)
{
  // fill histograms for pi0 according to evnet && pi0 candidates selection cuts

  if (!list) return;
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs     };  // PID
  TString namePID[kPIDs] =  { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"               };
  UInt_t  PIDBIT[kPIDs]   = { 0x0,  BIT(0), BIT(1), BIT(0)|BIT(1), BIT(2), BIT(3), BIT(2)|BIT(3)   };
                         // { kPtPi0, kEtaPi0, kPhiPi0, kAsyPi0, kAnglePi0, kInvMassPi0, kVarsPi0  };  // Pi0
  TString tName[kVarsPi0] = {   "Pt",   "Eta",   "Phi",   "Asy",   "Angle",   "InvMass"            };

  Int_t entries = caloClArr->GetEntriesFast();
  Int_t iMod = 0, jMod = 0, srcMod = 0;
  UInt_t iPIDBit = 0, jPIDBit = 0;
  Double_t vars[kVarsPi0];
  AliCaloClusterInfo *iCaloCluster = 0, *jCaloCluster = 0;
  TLorentzVector iGamma, jGamma;
  for(Int_t i=0; i<entries-1; i++) { // Loop calo cluster i
    iCaloCluster = (AliCaloClusterInfo*)caloClArr->At(i);
    iMod    = iCaloCluster->GetModule();
    iPIDBit = iCaloCluster->GetPIDBit();
    iGamma  = iCaloCluster->LorentzVector();
    
    for (Int_t j=i+1; j<entries; j++) { // Loop calo cluster j
      jCaloCluster = (AliCaloClusterInfo*)caloClArr->At(j);
      jMod    = jCaloCluster->GetModule();
      jPIDBit = jCaloCluster->GetPIDBit();
      jGamma  = jCaloCluster->LorentzVector();
      if (TMath::Abs(iMod-jMod)>1) { jCaloCluster = 0; continue; }

      if (iMod==1 && jMod==1 )      srcMod = 11; // Mod11
      else if (iMod==3 && jMod==3)  srcMod = 33; // Mod33
      else if (iMod==2 && jMod==2)  srcMod = 22; // Mod22
      else if (iMod==1 && jMod==2)  srcMod = 12; // Mod12
      else if (iMod==2 && jMod==3)  srcMod = 23; // Mod23

      vars[kPtPi0]      = (iGamma+jGamma).E();
      vars[kEtaPi0]     = (iGamma+jGamma).Eta();
      vars[kPhiPi0]     = (iGamma+jGamma).Phi() + TMath::TwoPi();
      vars[kInvMassPi0] = (iGamma+jGamma).M();
      vars[kAsyPi0]     = TMath::Abs((iGamma-jGamma).E())/(iGamma+jGamma).E();
      vars[kAnglePi0]   = jGamma.Angle(iGamma.Vect());

      ((TH2D*)list->FindObject(Form("hPi0_%s%s", tName[kEtaPi0].Data(), tName[kPhiPi0].Data())))->Fill(vars[kEtaPi0], vars[kPhiPi0]);
      for (Int_t k=0; k<kVarsPi0-2; k++) {
        ((TH2D*)list->FindObject(Form("hPi0_%s%s", tName[0].Data(), tName[k+1].Data())))->Fill(vars[0], vars[k+1]);
      }

      for (Int_t iPID=0; iPID<kPIDs; iPID++) {
        if ((iPIDBit & jPIDBit & PIDBIT[iPID])==PIDBIT[iPID]) {
          ((TH2D*)list->FindObject(Form("hPi0_%s%s_%s", tName[kEtaPi0].Data(), tName[kPhiPi0].Data(),
                                                        namePID[iPID].Data())))->Fill(vars[kEtaPi0], vars[kPhiPi0]);
          ((TH2D*)list->FindObject(Form("hPi0_%s%s_%s", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(),
                                                        namePID[iPID].Data())))->Fill(vars[kPtPi0], vars[kInvMassPi0]);
          ((TH2D*)list->FindObject(Form("hPi0_%s%s_%s_Mod%d", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(),
                                                              namePID[iPID].Data(), srcMod)))->Fill(vars[kPtPi0], vars[kInvMassPi0]);
          ((TH2D*)list->FindObject(Form("hPi0_%s%s_%s_cent%d", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(), 
                                                               namePID[iPID].Data(), cent)))->Fill(vars[kPtPi0], vars[kInvMassPi0]);
          for (Int_t k=0; k<kVarsPi0-2; k++) {
            ((TH2D*)list->FindObject(Form("hPi0_%s%s_%s", tName[0].Data(), tName[k+1].Data(), namePID[iPID].Data())))->Fill(vars[0], vars[k+1]);
          }
        }
      }
      jCaloCluster = 0;
    }// End loop j
    iCaloCluster = 0;
  } // End loop i

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosMixPi0(TList *list, TClonesArray* const iCaloClArr, TList* const eventList, Int_t cent)
{
  // fill histograms for mixed pi0 according to evnet && pi0 candidates selection cuts

  if (!list) return;
                        // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs     };  // PID
  TString namePID[kPIDs] = { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"               };
  UInt_t  PIDBIT[kPIDs]  = { 0x0,  BIT(0), BIT(1), BIT(0)|BIT(1), BIT(2), BIT(3), BIT(2)|BIT(3)   };
                            // { kPtMixPi0, kEtaMixPi0, kPhiMixPi0, kInvMassMixPi0, kVarsMixPi0   };   // Mix pi0
  TString tName[kVarsMixPi0] = {      "Pt",      "Eta",      "Phi",      "InvMass"                };

  Int_t iEntries = iCaloClArr->GetEntriesFast();
  Int_t iMod = 0, jMod = 0, srcMod = 0;
  UInt_t iPIDBit = 0, jPIDBit = 0;
  Double_t vars[kVarsMixPi0];
  AliCaloClusterInfo *iCaloCluster = 0, *jCaloCluster = 0;
  TLorentzVector iGamma, jGamma;
  TClonesArray *jCaloClArr = 0;
  for (Int_t i=0; i<iEntries; i++) { // Loop calo cluster i
    iCaloCluster=(AliCaloClusterInfo*)iCaloClArr->At(i);
    iMod    = iCaloCluster->GetModule();
    iPIDBit = iCaloCluster->GetPIDBit();
    iGamma  = iCaloCluster->LorentzVector();

    Int_t nev = eventList->GetSize();
    for(Int_t ev=0; ev<nev; ev++){ // Loop events for mixing
      jCaloClArr     = static_cast<TClonesArray*>(eventList->At(ev));
      Int_t jEntries = jCaloClArr->GetEntriesFast();
      for(Int_t j=0; j<jEntries; j++){ // Loop calo cluster j
        jCaloCluster=(AliCaloClusterInfo*)jCaloClArr->At(j);
        jMod    = jCaloCluster->GetModule();   if (TMath::Abs(iMod-jMod)>1) { jCaloCluster = 0; continue; }
        jPIDBit = jCaloCluster->GetPIDBit();
        jGamma  = jCaloCluster->LorentzVector();

        if (iMod==1 && jMod==1 )                                srcMod = 11; // Mod11
        else if (iMod==3 && jMod==3)                            srcMod = 33; // Mod33
        else if (iMod==2 && jMod==2)                            srcMod = 22; // Mod22
        else if ((iMod==1 && jMod==2) || (iMod==2 && jMod==1))  srcMod = 12; // Mod12
        else if ((iMod==2 && jMod==3) || (iMod==3 && jMod==2))  srcMod = 23; // Mod23

        vars[kPtMixPi0]      = (iGamma+jGamma).E();
        vars[kEtaMixPi0]     = (iGamma+jGamma).Eta();
        vars[kPhiMixPi0]     = (iGamma+jGamma).Phi() + TMath::TwoPi();
        vars[kInvMassMixPi0] = (iGamma+jGamma).M();

        ((TH2D*)list->FindObject(Form("hMixPi0_%s%s", tName[kEtaMixPi0].Data(), tName[kPhiMixPi0].Data())))->Fill(vars[kEtaMixPi0], vars[kPhiMixPi0]);
        for (Int_t iPID=0; iPID<kPIDs; iPID++) {
          if ((iPIDBit & jPIDBit & PIDBIT[iPID])==PIDBIT[iPID]) {
            ((TH2D*)list->FindObject(Form("hMixPi0_%s%s_%s", tName[kEtaMixPi0].Data(), tName[kPhiMixPi0].Data(),
                                                             namePID[iPID].Data())))->Fill(vars[kEtaMixPi0], vars[kPhiMixPi0]);
            ((TH2D*)list->FindObject(Form("hMixPi0_%s%s_%s", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(),
                                                             namePID[iPID].Data())))->Fill(vars[kPtMixPi0], vars[kInvMassMixPi0]);
            ((TH2D*)list->FindObject(Form("hMixPi0_%s%s_%s_Mod%d", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(),
                                                                   namePID[iPID].Data(), srcMod)))->Fill(vars[kPtMixPi0], vars[kInvMassMixPi0]);
            ((TH2D*)list->FindObject(Form("hMixPi0_%s%s_%s_cent%d", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(),
                                                                    namePID[iPID].Data(), cent)))->Fill(vars[kPtMixPi0], vars[kInvMassMixPi0]);
          }
        }
        jCaloCluster = 0;
      } // End loop j
      jCaloClArr = 0;
    } // End loop events
    iCaloCluster = 0;
  } // End loop i

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosMC(TList *list, AliMCEvent* const mcEvent, AliPHOSGeoUtils* const phosGeo, Int_t cent)
{
  // fill histograms for AOD MC
  if (!list) return;

  TString pName;
  const Int_t ntrks = mcEvent->GetNumberOfTracks();
  for(Int_t i=0; i<ntrks; i++) {
    AliAODMCParticle *pMC = (AliAODMCParticle*)mcEvent->GetTrack(i);
    Int_t pdg = pMC->GetPdgCode();
    if (pdg == 111)      pName = "Pi0";   // Pi0
    else if (pdg == 221) pName = "Eta";   // Eta
    else if (pdg == 22)  pName = "Gamma"; // Gamma
    else continue;

    //Primary particle
    Double_t pt     = pMC->Pt();
    Double_t r      = TMath::Sqrt(TMath::Power(pMC->Xv(), 2) + TMath::Power(pMC->Yv(), 2));
    Double_t y      = pMC->Y();
    Double_t phi    = pMC->Phi();   while (phi<0.) phi += TMath::TwoPi(); while (phi>TMath::TwoPi()) phi -= TMath::TwoPi();
//  Double_t weight = PrimaryParticleWeight(pdg, pt);

    ((TH2D*)list->FindObject(Form("hMC%s_PtVertex", pName.Data())))->Fill(pt, r);
    ((TH2D*)list->FindObject(Form("hMC%s_PtRapidity", pName.Data())))->Fill(pt, y);

    Int_t mod1 = 0, mod2 = 0;
    AliAODMCParticle *gamma1=0x0, *gamma2=0x0;
    if (r<1.) {
      ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_RCut", pName.Data(), cent)))->Fill(pt);
      ((TH1D*)list->FindObject(Form("hMC%s_Rapidity_cent%d_RCut", pName.Data(), cent)))->Fill(y);
      ((TH1D*)list->FindObject(Form("hMC%s_Phi_cent%d_RCut", pName.Data(), cent)))->Fill(phi);
      if (y<1.)   ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_WideY_RCut", pName.Data(), cent)))->Fill(pt); 
      if (y<0.15) ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_NarrY_RCut", pName.Data(), cent)))->Fill(pt); 
      if (pMC->GetNDaughters()==2) { // Do not account Dalitz decays
        gamma1 = (AliAODMCParticle*)mcEvent->GetTrack(pMC->GetDaughter(0));
        gamma2 = (AliAODMCParticle*)mcEvent->GetTrack(pMC->GetDaughter(1));
        //Number of pi0s decayed into acceptance
        mod1 = HitPHOSModule(gamma1, phosGeo);
        mod2 = HitPHOSModule(gamma2, phosGeo);
        if (mod1!=-1 && mod1!=2 && mod2!=-1 && mod2!=2 && TMath::Abs(mod1-mod2)<2) // !remove module 2
          ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_Acc_RCut", pName.Data(), cent)))->Fill(pt);
      }
    }
    if (pMC->IsPhysicalPrimary()) {
      ((TH2D*)list->FindObject(Form("hMC%s_PtVertex_IsPrimary", pName.Data())))->Fill(pt, r);
      ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_IsPrimary", pName.Data(), cent)))->Fill(pt);
      ((TH1D*)list->FindObject(Form("hMC%s_Rapidity_cent%d_IsPrimary", pName.Data(), cent)))->Fill(y);
      ((TH1D*)list->FindObject(Form("hMC%s_Phi_cent%d_IsPrimary", pName.Data(), cent)))->Fill(phi);
      if (y<1.)   ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_WideY_IsPrimary", pName.Data(), cent)))->Fill(pt); 
      if (y<0.15) ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_NarrY_IsPrimary", pName.Data(), cent)))->Fill(pt); 
      if (pMC->GetNDaughters()==2) { // Do not account Dalitz decays
        gamma1 = (AliAODMCParticle*)mcEvent->GetTrack(pMC->GetDaughter(0));
        gamma2 = (AliAODMCParticle*)mcEvent->GetTrack(pMC->GetDaughter(1));
        //Number of pi0s decayed into acceptance
        mod1 = HitPHOSModule(gamma1, phosGeo);
        mod2 = HitPHOSModule(gamma2, phosGeo);
        if (mod1!=-1 && mod1!=2 && mod2!=-1 && mod2!=2 && TMath::Abs(mod1-mod2)<2) // !remove module 2
          ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_Acc_IsPrimary", pName.Data(), cent)))->Fill(pt);
      }
    }
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosMC(TList *list, AliStack* const stack, AliPHOSGeoUtils* const phosGeo, Int_t cent)
{
  // fill histograms for ESD MC
  if (!list)  return;
  if (!stack) return;

  TString pName;
  const Int_t ntrks = stack->GetNtrack();
  for(Int_t i=0; i<ntrks; i++) {
    TParticle* pMC = stack->Particle(i);
    Int_t pdg = pMC->GetPdgCode();
    if (pdg == 111)      pName = "Pi0";   // Pi0
    else if (pdg == 221) pName = "Eta";   // Eta
    else if (pdg == 22)  pName = "Gamma"; // Gamma
    else continue;

    //Primary particle
    Double_t pt     = pMC->Pt();
    Double_t r      = pMC->R();
    Double_t y      = pMC->Y();
    Double_t phi    = pMC->Phi();   while (phi<0.) phi += TMath::TwoPi(); while (phi>TMath::TwoPi()) phi -= TMath::TwoPi();
//  Double_t weight = PrimaryParticleWeight(pdg, pt);

    ((TH2D*)list->FindObject(Form("hMC%s_PtVertex", pName.Data())))->Fill(pt, r);
    ((TH2D*)list->FindObject(Form("hMC%s_PtRapidity", pName.Data())))->Fill(pt, y);

    Int_t mod1 = 0, mod2 = 0;
    TParticle *gamma1 = 0x0, *gamma2 = 0x0;
    if (r<1.) { // R_Cut
      ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_RCut", pName.Data(), cent)))->Fill(pt);
      ((TH1D*)list->FindObject(Form("hMC%s_Rapidity_cent%d_RCut", pName.Data(), cent)))->Fill(y);
      ((TH1D*)list->FindObject(Form("hMC%s_Phi_cent%d_RCut", pName.Data(), cent)))->Fill(phi);
      if (y<0.5)  ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_WideY_RCut", pName.Data(), cent)))->Fill(pt);
      if (y<0.15) ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_NarrY_RCut", pName.Data(), cent)))->Fill(pt);
      if (pMC->GetNDaughters()==2) { // Do not account Dalitz decays
        gamma1 = stack->Particle(pMC->GetFirstDaughter());
        gamma2 = stack->Particle(pMC->GetLastDaughter());
        //Number of pi0s decayed into acceptance
        mod1 = HitPHOSModule(gamma1, phosGeo);
        mod2 = HitPHOSModule(gamma2, phosGeo);
        if (mod1!=-1 && mod1!=2 && mod2!=-1 && mod2!=2 && TMath::Abs(mod1-mod2)<2) // !remove module 2
          ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_Acc_RCut", pName.Data(), cent)))->Fill(pt);
      }
    }
    if (stack->IsPhysicalPrimary(i)) { // IsPrimary
      ((TH2D*)list->FindObject(Form("hMC%s_PtVertex_IsPrimary", pName.Data())))->Fill(pt, r);
      ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_IsPrimary", pName.Data(), cent)))->Fill(pt);
      ((TH1D*)list->FindObject(Form("hMC%s_Rapidity_cent%d_IsPrimary", pName.Data(), cent)))->Fill(y);
      ((TH1D*)list->FindObject(Form("hMC%s_Phi_cent%d_IsPrimary", pName.Data(), cent)))->Fill(phi);
      if (y<0.5)  ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_WideY_IsPrimary", pName.Data(), cent)))->Fill(pt);
      if (y<0.15) ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_NarrY_IsPrimary", pName.Data(), cent)))->Fill(pt);
      if (pMC->GetNDaughters()==2) { // Do not account Dalitz decays
        gamma1 = stack->Particle(pMC->GetFirstDaughter());
        gamma2 = stack->Particle(pMC->GetLastDaughter());
        //Number of pi0s decayed into acceptance
        mod1 = HitPHOSModule(gamma1, phosGeo);
        mod2 = HitPHOSModule(gamma2, phosGeo);
        if (mod1!=-1 && mod1!=2 && mod2!=-1 && mod2!=2 && TMath::Abs(mod1-mod2)<2) // !remove module 2
          ((TH1D*)list->FindObject(Form("hMC%s_Pt_cent%d_Acc_IsPrimary", pName.Data(), cent)))->Fill(pt);
      }
    }
  }

  return;
}

//________________________________________________________________________
Double_t AliPHOSpPbPi0Header::PrimaryParticleWeight(Int_t pdg, Double_t pt)
{

  Int_t type=0 ;
  if (pdg == 111 || TMath::Abs(pdg)==211) type = 1; // Pi0/Eta 
  else if (TMath::Abs(pdg)<1000)          type = 2; // Kaon-like
  else                                    type = 3; // baryons

  if (type==1) {
    if (fCentrality<5.)       // 0-5
      return (1.662990+1.140890*pt-0.192088*pt*pt)/(1.-0.806630*pt+0.304771*pt*pt)+0.141690*pt;
    else if (fCentrality<10.) // 5-10
      return (1.474351+0.791492*pt-0.066369*pt*pt)/(1.-0.839338*pt+0.317312*pt*pt)+0.093289*pt;
    else if (fCentrality<20.) // 10-20
      return (1.174728+0.959681*pt-0.137695*pt*pt)/(1.-0.788873*pt+0.299538*pt*pt)+0.128759*pt; 
    else if (fCentrality<40.) // 20-40
      return (0.927335+0.475349*pt+0.004364*pt*pt)/(1.-0.817966*pt+0.309787*pt*pt)+0.086899*pt; 
    else if (fCentrality<60.) // 40-60
      return (0.676878+0.190680*pt+0.077031*pt*pt)/(1.-0.790623*pt+0.305183*pt*pt)+0.064510*pt; 
    else if (fCentrality<80.) // 60-80
      return (0.684726-0.606262*pt+0.409963*pt*pt)/(1.-1.080061*pt+0.456933*pt*pt)+0.005151*pt; 
  }
  if (type==2) {
    if (fCentrality<5.)       // 0-5
      return (-0.417131+2.253936*pt-0.337731*pt*pt)/(1.-0.909892*pt+0.316820*pt*pt)+0.157312*pt;
    else if (fCentrality<10.) // 5-10
      return (-0.352275+1.844466*pt-0.248598*pt*pt)/(1.-0.897048*pt+0.316462*pt*pt)+0.132461*pt; 
    else if (fCentrality<20.) // 10-20
      return (-0.475481+1.975108*pt-0.336013*pt*pt)/(1.-0.801028*pt+0.276705*pt*pt)+0.188164*pt; 
    else if (fCentrality<40.) // 20-40
      return (-0.198954+1.068789*pt-0.103540*pt*pt)/(1.-0.848354*pt+0.299209*pt*pt)+0.112939*pt; 
    else if (fCentrality<60.) // 40-60
      return (-0.111052+0.664041*pt-0.019717*pt*pt)/(1.-0.804916*pt+0.300779*pt*pt)+0.095784*pt;
    else if (fCentrality<80.) // 60-80
      return (0.202788-0.439832*pt+0.564585*pt*pt)/(1.-1.254029*pt+0.679444*pt*pt)+0.016235*pt;
  }
  if (type==3) {
    if (fCentrality<5.)       // 0-5
      return (-1.312732+2.743568*pt-0.375775*pt*pt)/(1.-0.717533*pt+0.164694*pt*pt)+0.164445*pt;
    else if (fCentrality<10.) // 5-10
      return (-1.229425+2.585889*pt-0.330164*pt*pt)/(1.-0.715892*pt+0.167386*pt*pt)+0.133085*pt;
    else if (fCentrality<20.) // 10-20
      return (-1.135677+2.397489*pt-0.320355*pt*pt)/(1.-0.709312*pt+0.164350*pt*pt)+0.146095*pt;
    else if (fCentrality<40.) // 20-40
      return (-0.889993+1.928263*pt-0.220785*pt*pt)/(1.-0.715991*pt+0.174729*pt*pt)+0.095098*pt;
    else if (fCentrality<60.) // 40-60
      return (-0.539237+1.329118*pt-0.115439*pt*pt)/(1.-0.722906*pt+0.186832*pt*pt)+0.059267*pt;
    else if (fCentrality<80.) // 60-80
      return (-0.518126+1.327628*pt-0.130881*pt*pt)/(1.-0.665649*pt+0.184300*pt*pt)+0.081701*pt;
  }

  return 1.;
}

//________________________________________________________________________
Int_t AliPHOSpPbPi0Header::HitPHOSModule(AliAODMCParticle* const pMC, AliPHOSGeoUtils* const phosGeo)
{ 

  Int_t mod=0, relID[4]={0,0,0,0}; 
  Double_t x=0., z=0.;
  Double_t vtx[3]; pMC->XvYvZv(vtx);
  Double_t theta = pMC->Theta();
  Double_t phi   = pMC->Phi();
    
  if (!phosGeo->ImpactOnEmc(vtx, theta, phi, mod, z, x)) return -1;
  phosGeo->RelPosToRelId(mod, x, z, relID);
  if (fgUseFiducialCut) {
    const Int_t edgeX = 2;
    const Int_t edgeZ = 2;
    if (relID[2] >edgeX && relID[2]<(65-edgeX) && relID[3]>edgeZ && relID[3] <(57-edgeZ)) return relID[0];
    else return -1;
  }

  return relID[0];
}

//________________________________________________________________________
Int_t AliPHOSpPbPi0Header::HitPHOSModule(TParticle* const pMC, AliPHOSGeoUtils* const phosGeo)
{

  Int_t mod=0, relID[4]={0,0,0,0}; 
  Double_t x=0., z=0.;
  Double_t vtx[3] = { pMC->Vx(), pMC->Vy(), pMC->Vz() };
  Double_t theta = pMC->Theta();
  Double_t phi   = pMC->Phi();

  if (!phosGeo->ImpactOnEmc(vtx, theta, phi, mod, z, x)) return -1;
  phosGeo->RelPosToRelId(mod, x, z, relID);
  if (fgUseFiducialCut) {
    const Int_t edgeX = 2;
    const Int_t edgeZ = 2;
    if (relID[2] >edgeX && relID[2]<(65-edgeX) && relID[3]>edgeZ && relID[3] <(57-edgeZ)) return relID[0];
    else return -1;
  } 

  return relID[0];
}
