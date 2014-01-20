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
#include "AliAnalysisUtils.h"
#include "AliVCaloCells.h"
#include "AliPHOSGeoUtils.h"
#include "AliCentrality.h"

#include "AliCaloClusterInfo.h"
#include "AliPHOSpPbPi0Header.h"

class TNamed;

ClassImp(AliPHOSpPbPi0Header)

Bool_t   AliPHOSpPbPi0Header::fgIsMC           = kFALSE;
Bool_t   AliPHOSpPbPi0Header::fgUseFiducialCut = kFALSE;
Int_t    AliPHOSpPbPi0Header::fgNCent          = 10;
Double_t AliPHOSpPbPi0Header::fgCuts[3]        = { 10., 0., 100. };

//_____________________________________________________________________________
AliPHOSpPbPi0Header::AliPHOSpPbPi0Header() :
  TNamed(), fFiredTriggerClass(), fSelMask(0), fIsVertexOK(kFALSE), fIsPileup(kFALSE), fCentrality(0.)
{
  //
  // default constructor
  //
  for (Int_t i=3; i--;) fVtx[i] = 0.;
}

//_____________________________________________________________________________
AliPHOSpPbPi0Header::AliPHOSpPbPi0Header(const AliPHOSpPbPi0Header &src) :
  TNamed(), fFiredTriggerClass(src.fFiredTriggerClass), fSelMask(src.fSelMask), fIsVertexOK(src.fIsVertexOK),
  fIsPileup(src.fIsPileup), fCentrality(src.fCentrality)
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
  fIsVertexOK                   = src.fIsVertexOK;
  fIsPileup                     = src.fIsPileup;
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

  fSelMask      = handler->IsEventSelected();
  fIsVertexOK   = CheckEventVertex(event);

  AliCentrality *cent = event->GetCentrality();
  if (cent) fCentrality = cent->GetCentralityPercentile("V0A");
  else      fCentrality = -1.;

//this->SetTitle(vertex->GetTitle());
  return;
}

//_____________________________________________________________________________
Bool_t AliPHOSpPbPi0Header::IsSelected()
{
  // select event according to the event selection cuts
  if (TMath::Abs(fVtx[2])>fgCuts[0])                   return kFALSE; // Vtxz cut
  if (fCentrality<fgCuts[1] || fCentrality>fgCuts[2])  return kFALSE; // centrality selection
  if (!fIsVertexOK)                                    return kFALSE; // pA vertex cut

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPHOSpPbPi0Header::CheckEventVertex(AliVEvent* const event)
{
  // check event vertex

  // set event basic info
  fFiredTriggerClass       = event->GetFiredTriggerClasses();
  const AliVVertex* trkVtx = event->GetPrimaryVertex();
  if (!trkVtx)                               return kFALSE;
  trkVtx->GetXYZ(fVtx);

  AliAnalysisUtils *utils = new AliAnalysisUtils();
  if (!utils->IsVertexSelected2013pA(event)) return kFALSE;
  fIsPileup = utils->IsPileUpEvent(event);

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
void AliPHOSpPbPi0Header::CreateHistosEvent(TList *listQA)
{
  // create histograms at event level

  const Int_t nhs    = 3;
  TString tName[nhs] = { "Centrality", "Vz", "Pileup" };
  Int_t   nbins[nhs] = {         110 , 500 ,     500  };
  Double_t xlow[nhs] = {         -10., -25.,     -25. };
  Double_t  xup[nhs] = {         100.,  25.,      25. };

  TString hName;
  TH1D   *histo = 0;
  for (Int_t i=0; i<nhs; i++) {
    hName = Form("hEvent_%s", tName[i].Data());
    histo = new TH1D(hName.Data(), hName.Data(), nbins[i], xlow[i], xup[i]);
    histo->Sumw2(); listQA->Add(histo); histo = 0;

    hName = Form("hEventSel_%s", tName[i].Data());
    histo = new TH1D(hName.Data(), hName.Data(), nbins[i], xlow[i], xup[i]);
    histo->Sumw2(); listQA->Add(histo); histo = 0;
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosCaloCellsQA(TList *listQA)
{
  // create QA histograms for CaloCells

  TH1D *histo1       = 0;
  TH2D *histo2       = 0;
  const Int_t nhs    = 3;
  const Int_t nMod   = 3;
  TString hName = "hCaloCellsQA_E";
  histo1 = new TH1D(hName, hName, 300, 0., 30.);
  histo1->Sumw2(); listQA->Add(histo1); histo1 = 0;
  for (Int_t i=0; i<nMod; i++) {
    hName  = Form("hCaloCellsQA_E_Mod%d", i+1);
    histo1 = new TH1D(hName.Data(), hName.Data(), 300, 0., 30.);
    histo1->Sumw2(); listQA->Add(histo1); histo1 = 0;
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
      hName  = Form("hCaloCellsQA_%s", tName[i].Data());
      histo2 = new TH2D(hName.Data(), hName.Data(), nbinsx[i], xlow[i], xup[i], nbinsy[i], ylow[i], yup[i]);
      histo2->Sumw2(); listQA->Add(histo2); histo2 = 0;
    } else {
      for (Int_t j=0; j<nMod; j++) {
        hName  = Form("hCaloCellsQA_%s_Mod%d", tName[i].Data(), j+1);
        histo2 = new TH2D(hName.Data(), hName.Data(), nbinsx[i], xlow[i], xup[i], nbinsy[i], ylow[i], yup[i]);
        histo2->Sumw2(); listQA->Add(histo2); histo2 = 0;
      }
    }
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosCaloCluster(TList *listQA)
{
  // create histograms for CaloCluster

  const Int_t nMods = 3;
  const Int_t nTRUs = 8;
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs };   // PID
  TString namePID[kPIDs] =  { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"           };
                         // { kPtClu, kEtaClu, kPhiClu, kM02Clu, kM20Clu, kTOFClu, kNCellsClu, kNClustersClu, kVarsClu };   // clusters
  TString tName[kVarsClu] = {   "Pt",   "Eta",   "Phi",   "M02",   "M20",   "TOF",   "NCells",   "NClusters"           };
  Int_t    bins[kVarsClu] = {   500 ,  300   ,   120  ,    200 ,    200 ,    400 ,       200 ,          20             };
  Double_t xMin[kVarsClu] = {     0.,   -0.15,     4.4,      0.,      0.,   -2e-6,         0.,           0.            };
  Double_t xMax[kVarsClu] = {    50.,    0.15,     5.6,      2.,      2.,    2e-6,      1000.,          20.            };

  TH1I *hncells = 0;
  TH1D *histo1  = 0;
  TH2D *histo2  = 0;
  TString hName;

  hName   = Form("hCaloCluster_%s", tName[kNCellsClu].Data());
  hncells = new TH1I(hName.Data(), hName.Data(), bins[kNCellsClu], (Int_t)xMin[kNCellsClu], (Int_t)xMax[kNCellsClu]);
  hncells->Sumw2(); listQA->Add(hncells); hncells = 0;

  hName   = Form("hCaloCluster_%s", tName[kNClustersClu].Data());
  hncells = new TH1I(hName.Data(), hName.Data(), bins[kNClustersClu], (Int_t)xMin[kNClustersClu], (Int_t)xMax[kNClustersClu]);
  hncells->Sumw2(); listQA->Add(hncells); hncells = 0;

  for (Int_t iPID=0; iPID<kPIDs; iPID++) {
    hName  = Form("hCaloCluster_%s_%s", tName[kPtClu].Data(), namePID[iPID].Data());
    histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtClu], xMin[kPtClu], xMax[kPtClu]);
    histo1->Sumw2(); listQA->Add(histo1); histo1 = 0;

    for (Int_t icent=0; icent<fgNCent; icent++) {
      hName  = Form("hCaloCluster_%s_%s_cent%d", tName[kPtClu].Data(), namePID[iPID].Data(), icent);
      histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtClu], xMin[kPtClu], xMax[kPtClu]);
      histo1->Sumw2(); listQA->Add(histo1); histo1 = 0;
    }

    hName  = Form("hCaloCluster_%s%s_%s", tName[kEtaClu].Data(), tName[kPhiClu].Data(), namePID[iPID].Data());
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[kEtaClu], xMin[kEtaClu], xMax[kEtaClu], bins[kPhiClu], xMin[kPhiClu], xMax[kPhiClu]);
    histo2->Sumw2(); listQA->Add(histo2); histo2 = 0;
    for (Int_t i=0; i<kVarsClu-3; i++) {
      hName  = Form("hCaloCluster_%s%s_%s", tName[0].Data(), tName[i+1].Data(), namePID[iPID].Data());
      histo2 = new TH2D(hName.Data(), hName.Data(), bins[0], xMin[0], xMax[0], bins[i+1], xMin[i+1], xMax[i+1]);
      histo2->Sumw2(); listQA->Add(histo2); histo2 = 0;
    }
  }

  for (Int_t iMod=0; iMod<nMods; iMod++) {
    hName   = Form("hCaloCluster_%s_Mod%d", tName[kNCellsClu].Data(), iMod+1);
    hncells = new TH1I(hName.Data(), hName.Data(), bins[kNCellsClu], (Int_t)xMin[kNCellsClu], (Int_t)xMax[kNCellsClu]);
    hncells->Sumw2(); listQA->Add(hncells); hncells = 0;

    hName   = Form("hCaloCluster_%s_Mod%d", tName[kNClustersClu].Data(), iMod+1);
    hncells = new TH1I(hName.Data(), hName.Data(), bins[kNClustersClu], (Int_t)xMin[kNClustersClu], (Int_t)xMax[kNClustersClu]);
    hncells->Sumw2(); listQA->Add(hncells); hncells = 0;

    for (Int_t iTRU=0; iTRU<nTRUs; iTRU++) {
      hName  = Form("hCaloCluster_%s_Mod%d_TRU%d", tName[kPtClu].Data(), iMod+1, iTRU+1);
      histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtClu], xMin[kPtClu], xMax[kPtClu]);
      histo1->Sumw2(); listQA->Add(histo1); histo1 = 0;
    }

    for (Int_t iPID=0; iPID<kPIDs; iPID++) {
      hName  = Form("hCaloCluster_%s_Mod%d_%s", tName[kPtClu].Data(), iMod+1, namePID[iPID].Data());
      histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtClu], xMin[kPtClu], xMax[kPtClu]);
      histo1->Sumw2(); listQA->Add(histo1); histo1 = 0;
    }
  }

  hName  = Form("hCaloCluster_%s%s", tName[kEtaClu].Data(), tName[kPhiClu].Data());
  histo2 = new TH2D(hName.Data(), hName.Data(), bins[kEtaClu], xMin[kEtaClu], xMax[kEtaClu], bins[kPhiClu], xMin[kPhiClu], xMax[kPhiClu]);
  histo2->Sumw2(); listQA->Add(histo2); histo2 = 0;
  for (Int_t i=0; i<kVarsClu-3; i++) {
    hName  = Form("hCaloCluster_%s%s", tName[0].Data(), tName[i+1].Data());
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[0], xMin[0], xMax[0], bins[i+1], xMin[i+1], xMax[i+1]);
    histo2->Sumw2(); listQA->Add(histo2); histo2 = 0;
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosPi0(TList *listRD)
{
  // create histograms for Pi0

  const Int_t nComb = 2;
  TString srcMod[nComb]  =  { "Mod11", "Mod33" }; 
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs      };  // PID
  TString namePID[kPIDs] =  { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"                };
                         // { kPtPi0, kEtaPi0, kPhiPi0, kAsyPi0,   kAnglePi0, kInvMassPi0, kVarsPi0 };  // pi0
  TString tName[kVarsPi0] = {   "Pt",   "Eta",   "Phi",   "Asy",     "Angle",   "InvMass"           };
  Int_t    bins[kVarsPi0] = {   500 ,  300   ,   120  ,    100 ,        180 ,       1000            };
  Double_t xMin[kVarsPi0] = {     0.,   -0.15,     4.4,      0.,          0.,          0.           };
  Double_t xMax[kVarsPi0] = {    50.,    0.15,     5.6,      1., TMath::Pi(),          1.           };

  const Int_t    nPtbins     = 32;
  Double_t ptbins[nPtbins+1] = {  0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,
                                  2.8,  3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.5,  5.0,  5.5,
                                  6.0,  7.0,  8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0,
                                 30.0, 36.0, 44.0 };

  TH2D   *histo = 0;
  TString hName;

  hName = Form("hPi0_%s%s", tName[kEtaPi0].Data(), tName[kPhiPi0].Data());
  histo = new TH2D(hName.Data(), hName.Data(), bins[kEtaPi0], xMin[kEtaPi0], xMax[kEtaPi0], bins[kPhiPi0], xMin[kPhiPi0], xMax[kPhiPi0]);
  histo->Sumw2(); listRD->Add(histo); histo = 0;
  for (Int_t i=0; i<kVarsPi0-2; i++) {
    hName = Form("hPi0_%s%s", tName[0].Data(), tName[i+1].Data());
    histo = new TH2D(hName.Data(), hName.Data(), bins[0], xMin[0], xMax[0], bins[i+1], xMin[i+1], xMax[i+1]);
    histo->Sumw2(); listRD->Add(histo); histo = 0;
  }

  for (Int_t iPID=0; iPID<kPIDs; iPID++) {
    hName = Form("hPi0_%s%s_%s", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(), namePID[iPID].Data());
    histo = new TH2D(hName.Data(), hName.Data(), nPtbins, ptbins, bins[kInvMassPi0], xMin[kInvMassPi0], xMax[kInvMassPi0]);
    histo->Sumw2(); listRD->Add(histo); histo = 0;

    hName = Form("hPi0_%s%s_%s", tName[kEtaPi0].Data(), tName[kPhiPi0].Data(), namePID[iPID].Data());
    histo = new TH2D(hName.Data(), hName.Data(), bins[kEtaPi0], xMin[kEtaPi0], xMax[kEtaPi0], bins[kPhiPi0], xMin[kPhiPi0], xMax[kPhiPi0]);
    histo->Sumw2(); listRD->Add(histo); histo = 0;
    for (Int_t icent=0; icent<fgNCent; icent++) {
      hName = Form("hPi0_%s%s_%s_cent%d", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(), namePID[iPID].Data(), icent);
      histo = new TH2D(hName.Data(), hName.Data(), nPtbins, ptbins, bins[kInvMassPi0], xMin[kInvMassPi0], xMax[kInvMassPi0]);
      histo->Sumw2(); listRD->Add(histo); histo = 0;
    }

    for (Int_t i=0; i<kVarsPi0-2; i++) {
      hName = Form("hPi0_%s%s_%s", tName[0].Data(), tName[i+1].Data(), namePID[iPID].Data());
      histo = new TH2D(hName.Data(), hName.Data(), bins[0], xMin[0], xMax[0], bins[i+1], xMin[i+1], xMax[i+1]);
      histo->Sumw2(); listRD->Add(histo); histo = 0;
    }
    for (Int_t iComb=0; iComb<nComb; iComb++) {
      hName = Form("hPi0_%s%s_%s_%s", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(), namePID[iPID].Data(), srcMod[iComb].Data());
      histo = new TH2D(hName.Data(), hName.Data(), nPtbins, ptbins, bins[kInvMassPi0], xMin[kInvMassPi0], xMax[kInvMassPi0]);
      histo->Sumw2(); listRD->Add(histo); histo = 0;
    }
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosMixPi0(TList *listRD)
{
  // create Histograms for Mixed Pi0

  const Int_t nComb = 2;
  TString  srcMod[nComb]     = { "Mod11", "Mod33" };
                             //{ kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs    };  // PID
  TString namePID[kPIDs]     = { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"              };
                             //{ kPtMixPi0, kEtaMixPi0, kPhiMixPi0, kInvMassMixPi0, kVarsMixPi0      };  // Mix Pi0
  TString tName[kVarsMixPi0] = {      "Pt",      "Eta",      "Phi",      "InvMass"                   };
  Int_t    bins[kVarsMixPi0] = {      500 ,     300   ,      120  ,          1000                    };
  Double_t xMin[kVarsMixPi0] = {        0.,      -0.15,        4.4,             0.                   };
  Double_t xMax[kVarsMixPi0] = {       50.,       0.15,        5.6,             1.                   };

  const Int_t    nPtbins     = 32;
  Double_t ptbins[nPtbins+1] = {  0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,
                                  2.8,  3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.5,  5.0,  5.5,
                                  6.0,  7.0,  8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0,
                                 30.0, 36.0, 44.0 };

  TH2D   *histo = 0;
  TString hName;

  hName = Form("hMixPi0_%s%s", tName[kEtaMixPi0].Data(), tName[kPhiMixPi0].Data());
  histo = new TH2D(hName.Data(), hName.Data(), bins[kEtaMixPi0], xMin[kEtaMixPi0], xMax[kEtaMixPi0],
                                               bins[kPhiMixPi0], xMin[kPhiMixPi0], xMax[kPhiMixPi0]);
  histo->Sumw2(); listRD->Add(histo); histo = 0;
  
  for (Int_t iPID=0; iPID<kPIDs; iPID++) {
    hName = Form("hMixPi0_%s%s_%s", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(), namePID[iPID].Data());
    histo = new TH2D(hName.Data(), hName.Data(), nPtbins, ptbins, bins[kInvMassMixPi0], xMin[kInvMassMixPi0], xMax[kInvMassMixPi0]);
    histo->Sumw2(); listRD->Add(histo); histo = 0;
    hName = Form("hMixPi0_%s%s_%s", tName[kEtaMixPi0].Data(), tName[kPhiMixPi0].Data(), namePID[iPID].Data());
    histo = new TH2D(hName.Data(), hName.Data(), bins[kEtaMixPi0], xMin[kEtaMixPi0], xMax[kEtaMixPi0], 
                                                 bins[kPhiMixPi0], xMin[kPhiMixPi0], xMax[kPhiMixPi0]);
    histo->Sumw2(); listRD->Add(histo); histo = 0;
    for (Int_t icent=0; icent<fgNCent; icent++) {
      hName = Form("hMixPi0_%s%s_%s_cent%d", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(), namePID[iPID].Data(), icent);
      histo = new TH2D(hName.Data(), hName.Data(), nPtbins, ptbins, bins[kInvMassMixPi0], xMin[kInvMassMixPi0], xMax[kInvMassMixPi0]);
      histo->Sumw2(); listRD->Add(histo); histo = 0;
    } 
    for (Int_t iComb=0; iComb<nComb; iComb++) {
      hName = Form("hMixPi0_%s%s_%s_%s", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(), namePID[iPID].Data(), srcMod[iComb].Data());
      histo = new TH2D(hName.Data(), hName.Data(), nPtbins, ptbins, bins[kInvMassMixPi0], xMin[kInvMassMixPi0], xMax[kInvMassMixPi0]);
      histo->Sumw2(); listRD->Add(histo); histo = 0;
    }
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::CreateHistosMC(TList *listMC)
{
  // create Histograms for MC

  const Int_t acc    = 2;
  const char *aName[acc]   = {"Pi0InAcc", "GammaInAcc"};
  const Int_t types  = 5;
  const char *pName[types] = {"InclusivePi0", "PrimaryPi0", "2ndPi0FromK0s", "2ndPi0FromMaterials", "2ndPi0FromOtherDecays"};
                         //  { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs     };  // PID
  TString namePID[kPIDs]   = { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"               };
                         //  { kPtMC, kRadiusMC, kRapidityMC,          kPhiMC, kInvMassMC, kVarsMC  };  // MC
  TString tName[kVarsMC]   = {  "Pt",  "Radius",  "Rapidity",           "Phi",  "InvMass"           };
  Int_t    bins[kVarsMC]   = {  500 ,      1000,        200 ,            360 ,       500            };
  Double_t xMin[kVarsMC]   = {    0.,        0.,         -1.,              0.,         0.           };
  Double_t xMax[kVarsMC]   = {   50.,      500.,          1.,  TMath::TwoPi(),         0.5          };

  const Int_t    nPtbins     = 32;
  Double_t ptbins[nPtbins+1] = {  0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,
                                  2.8,  3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.5,  5.0,  5.5,
                                  6.0,  7.0,  8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0,
                                 30.0, 36.0, 44.0 };

  TH1D   *histo1 = 0;
  TH2D   *histo2 = 0;
  TString hName;

  for (Int_t iType=0; iType<types; iType++) {
    hName  = Form("hMC%s_%s%s", pName[iType], tName[kPtMC].Data(), tName[kRadiusMC].Data());
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC], bins[kRadiusMC], xMin[kRadiusMC], xMax[kRadiusMC]);
    histo2->Sumw2(); listMC->Add(histo2); histo2 = 0;
    hName  = Form("hMC%s_%s%s", pName[iType], tName[kPtMC].Data(), tName[kRapidityMC].Data());
    histo2 = new TH2D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC], bins[kRapidityMC], xMin[kRapidityMC], xMax[kRapidityMC]);
    histo2->Sumw2(); listMC->Add(histo2); histo2 = 0;
    for (Int_t iAcc=0; iAcc<acc; iAcc++) {
      hName  = Form("hMC%s_%s_%s", pName[iType], tName[kPtMC].Data(), aName[iAcc]);
      histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC]);
      histo1->Sumw2(); listMC->Add(histo1); histo1 = 0;
    }
    for (Int_t iPID=0; iPID<kPIDs; iPID++) {
      hName  = Form("hMC%s_%s%s_%s_Reco", pName[iType], tName[kPtMC].Data(), tName[kInvMassMC].Data(), namePID[iPID].Data());
      histo2 = new TH2D(hName.Data(), hName.Data(), nPtbins, ptbins, bins[kInvMassMC], xMin[kInvMassMC], xMax[kInvMassMC]);
      histo2->Sumw2(); listMC->Add(histo2); histo2 = 0;
    }

    for (Int_t iCent=0; iCent<fgNCent; iCent++) {
      hName  = Form("hMC%s_%s%s_cent%d", pName[iType], tName[kPtMC].Data(), tName[kRadiusMC].Data(), iCent);
      histo2 = new TH2D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC], bins[kRadiusMC], xMin[kRadiusMC], xMax[kRadiusMC]);
      histo2->Sumw2(); listMC->Add(histo2); histo2 = 0;
      hName  = Form("hMC%s_%s%s_cent%d", pName[iType], tName[kPtMC].Data(), tName[kRapidityMC].Data(), iCent);
      histo2 = new TH2D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC], bins[kRapidityMC], xMin[kRapidityMC], xMax[kRapidityMC]);
      histo2->Sumw2(); listMC->Add(histo2); histo2 = 0;
      for (Int_t iAcc=0; iAcc<acc; iAcc++) {
        hName  = Form("hMC%s_%s_%s_cent%d", pName[iType], tName[kPtMC].Data(), aName[iAcc], iCent);
        histo1 = new TH1D(hName.Data(), hName.Data(), bins[kPtMC], xMin[kPtMC], xMax[kPtMC]);
        histo1->Sumw2(); listMC->Add(histo1); histo1 = 0;
      }
      for (Int_t iPID=0; iPID<kPIDs; iPID++) {
        hName  = Form("hMC%s_%s%s_%s_Reco_cent%d", pName[iType], tName[kPtMC].Data(), tName[kInvMassMC].Data(), namePID[iPID].Data(), iCent);
        histo2 = new TH2D(hName.Data(), hName.Data(), nPtbins, ptbins, bins[kInvMassMC], xMin[kInvMassMC], xMax[kInvMassMC]);
        histo2->Sumw2(); listMC->Add(histo2); histo2 = 0;
      }
    }
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosEvent(TList *listQA)
{
  // fill histograms at event level according to event selection cuts

  if (!listQA) return;

  const Int_t nhs      = 2;
  TString tName[nhs+1] = { "Centrality",    "Vz", "Pileup" };
  Double_t dist[nhs]   = {  fCentrality, fVtx[2]           };
  for (Int_t i=nhs; i--;) {
    ((TH1D*)listQA->FindObject(Form("hEvent_%s",tName[i].Data())))->Fill(dist[i]);
    if (this->IsSelected()) ((TH1D*)listQA->FindObject(Form("hEventSel_%s",tName[i].Data())))->Fill(dist[i]);
  }
  if (fIsPileup) {
    ((TH1D*)listQA->FindObject(Form("hEvent_%s",tName[nhs].Data())))->Fill(dist[nhs-1]);
    if (this->IsSelected())  ((TH1D*)listQA->FindObject(Form("hEventSel_%s",tName[nhs].Data())))->Fill(dist[nhs-1]);
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosCaloCellsQA(TList *listQA, AliVCaloCells* const cells, AliPHOSGeoUtils* const phosGeo)
{
  // fill QA histograms for calocells according to evnet selection cuts

  if (!listQA) return;

  const Int_t nhs    = 3;
  TString tName[nhs] = { "CentNCells", "Nxz", "Exz" };

  ((TH2D*)listQA->FindObject(Form("hCaloCellsQA_%s", tName[0].Data())))->Fill(fCentrality, cells->GetNumberOfCells());
  for (Int_t iCell=0; iCell<cells->GetNumberOfCells(); iCell++) {
    Int_t relID[4] = {0,0,0,0};
    phosGeo->AbsToRelNumbering(cells->GetCellNumber(iCell),relID); // Int_t mod   = relID[0]; Int_t cellX = relID[2]; Int_t cellZ = relID[3];
    Double_t wight[2] = { 1., cells->GetAmplitude(iCell)};
    ((TH1D*)listQA->FindObject(Form("hCaloCellsQA_%s", "E")))->Fill(wight[1]);
    ((TH1D*)listQA->FindObject(Form("hCaloCellsQA_%s_Mod%d", "E", relID[0])))->Fill(wight[1]);
    for (Int_t i=1; i<nhs; i++)
      ((TH2D*)listQA->FindObject(Form("hCaloCellsQA_%s_Mod%d", tName[i].Data(), relID[0])))->Fill(relID[2], relID[3], wight[i-1]);
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosCaloCluster(TList *listQA, TClonesArray* const caloClArr, Int_t cent)
{
  // fill histograms for calocluster according to evnet && calocluster candidates selection cuts

  if (!listQA)  return;
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs       };  // PID
  TString namePID[kPIDs]  = { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"                 };
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
 
    ((TH1I*)listQA->FindObject(Form("hCaloCluster_%s", tName[kNCellsClu].Data())))->Fill((Int_t)vars[kNCellsClu]);
    ((TH1I*)listQA->FindObject(Form("hCaloCluster_%s_Mod%d", tName[kNCellsClu].Data(), iMod)))->Fill((Int_t)vars[kNCellsClu]);
    ((TH1D*)listQA->FindObject(Form("hCaloCluster_%s_Mod%d_TRU%d", tName[kPtClu].Data(), iMod, iTRU)))->Fill(vars[kPtClu]);

    for (Int_t iPID=0; iPID<kPIDs; iPID++) {
      if ((pidBit&PIDBIT[iPID])==PIDBIT[iPID]) {
        ((TH1D*)listQA->FindObject(Form("hCaloCluster_%s_%s", tName[kPtClu].Data(), namePID[iPID].Data())))->Fill(vars[kPtClu]);
        ((TH1D*)listQA->FindObject(Form("hCaloCluster_%s_Mod%d_%s", tName[kPtClu].Data(), iMod, namePID[iPID].Data())))->Fill(vars[kPtClu]);
        ((TH1D*)listQA->FindObject(Form("hCaloCluster_%s_%s_cent%d", tName[kPtClu].Data(), namePID[iPID].Data(), cent)))->Fill(vars[kPtClu]);
        ((TH2D*)listQA->FindObject(Form("hCaloCluster_%s%s_%s", tName[kEtaClu].Data(), tName[kPhiClu].Data(),
                                                              namePID[iPID].Data())))->Fill(vars[kEtaClu], vars[kPhiClu]);
        for (Int_t j=0; j<kVarsClu-3; j++)
          ((TH2D*)listQA->FindObject(Form("hCaloCluster_%s%s_%s", tName[0].Data(), tName[j+1].Data(), namePID[iPID].Data())))->Fill(vars[0], vars[j+1]);
      }
    }

    ((TH2D*)listQA->FindObject(Form("hCaloCluster_%s%s", tName[kEtaClu].Data(), tName[kPhiClu].Data())))->Fill(vars[kEtaClu], vars[kPhiClu]);
    for (Int_t j=0; j<kVarsClu-3; j++) {
      ((TH2D*)listQA->FindObject(Form("hCaloCluster_%s%s", tName[0].Data(), tName[j+1].Data())))->Fill(vars[0], vars[j+1]);
    }
  } // End loop calo cluster

  ((TH1I*)listQA->FindObject(Form("hCaloCluster_%s", tName[kNClustersClu].Data())))->Fill(entries);
  for (Int_t jMod=1; jMod<4; jMod++) {
    if (jMod==2) continue;
    ((TH1I*)listQA->FindObject(Form("hCaloCluster_%s_Mod%d", tName[kNClustersClu].Data(), jMod)))->Fill(nclusters[jMod-1]);
  }

  return;
}

//_____________________________________________________________________________
void AliPHOSpPbPi0Header::FillHistosPi0(TList *listRD, TClonesArray* const caloClArr, Int_t cent)
{
  // fill histograms for pi0 according to evnet && pi0 candidates selection cuts

  if (!listRD) return;
                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs     };  // PID
  TString namePID[kPIDs]  = { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"               };
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

      if (iMod==1 && jMod==1)       srcMod = 11; // Mod11
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

      ((TH2D*)listRD->FindObject(Form("hPi0_%s%s", tName[kEtaPi0].Data(), tName[kPhiPi0].Data())))->Fill(vars[kEtaPi0], vars[kPhiPi0]);
      for (Int_t k=0; k<kVarsPi0-2; k++) {
        ((TH2D*)listRD->FindObject(Form("hPi0_%s%s", tName[0].Data(), tName[k+1].Data())))->Fill(vars[0], vars[k+1]);
      }

      for (Int_t iPID=0; iPID<kPIDs; iPID++) {
        if ((iPIDBit & jPIDBit & PIDBIT[iPID])==PIDBIT[iPID]) {
          ((TH2D*)listRD->FindObject(Form("hPi0_%s%s_%s", tName[kEtaPi0].Data(), tName[kPhiPi0].Data(),
                                                        namePID[iPID].Data())))->Fill(vars[kEtaPi0], vars[kPhiPi0]);
          ((TH2D*)listRD->FindObject(Form("hPi0_%s%s_%s", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(),
                                                        namePID[iPID].Data())))->Fill(vars[kPtPi0], vars[kInvMassPi0]);
          ((TH2D*)listRD->FindObject(Form("hPi0_%s%s_%s_Mod%d", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(),
                                                              namePID[iPID].Data(), srcMod)))->Fill(vars[kPtPi0], vars[kInvMassPi0]);
          ((TH2D*)listRD->FindObject(Form("hPi0_%s%s_%s_cent%d", tName[kPtPi0].Data(), tName[kInvMassPi0].Data(), 
                                                               namePID[iPID].Data(), cent)))->Fill(vars[kPtPi0], vars[kInvMassPi0]);
          for (Int_t k=0; k<kVarsPi0-2; k++) {
            ((TH2D*)listRD->FindObject(Form("hPi0_%s%s_%s", tName[0].Data(), tName[k+1].Data(), namePID[iPID].Data())))->Fill(vars[0], vars[k+1]);
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
void AliPHOSpPbPi0Header::FillHistosMixPi0(TList *listRD, TClonesArray* const iCaloClArr, TList* const eventList, Int_t cent)
{
  // fill histograms for mixed pi0 according to evnet && pi0 candidates selection cuts

  if (!listRD) return;
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

        ((TH2D*)listRD->FindObject(Form("hMixPi0_%s%s", tName[kEtaMixPi0].Data(), tName[kPhiMixPi0].Data())))->Fill(vars[kEtaMixPi0], vars[kPhiMixPi0]);
        for (Int_t iPID=0; iPID<kPIDs; iPID++) {
          if ((iPIDBit & jPIDBit & PIDBIT[iPID])==PIDBIT[iPID]) {
            ((TH2D*)listRD->FindObject(Form("hMixPi0_%s%s_%s", tName[kEtaMixPi0].Data(), tName[kPhiMixPi0].Data(),
                                                               namePID[iPID].Data())))->Fill(vars[kEtaMixPi0], vars[kPhiMixPi0]);
            ((TH2D*)listRD->FindObject(Form("hMixPi0_%s%s_%s", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(),
                                                               namePID[iPID].Data())))->Fill(vars[kPtMixPi0], vars[kInvMassMixPi0]);
            ((TH2D*)listRD->FindObject(Form("hMixPi0_%s%s_%s_Mod%d", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(),
                                                                     namePID[iPID].Data(), srcMod)))->Fill(vars[kPtMixPi0], vars[kInvMassMixPi0]);
            ((TH2D*)listRD->FindObject(Form("hMixPi0_%s%s_%s_cent%d", tName[kPtMixPi0].Data(), tName[kInvMassMixPi0].Data(),
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
void AliPHOSpPbPi0Header::FillHistosMC(TList *listMC, AliStack* const stack, TClonesArray* const caloClArr, AliPHOSGeoUtils* const phosGeo, Int_t cent)
{
  // fill histograms for ESD MC
  if (!listMC)  return;
  if (!stack)   return;

                         // { kAll,  kCpv,  kDisp,  kBoth,  kCpv2,  kDisp2,  kBoth2,     kPIDs     };  // PID
  TString namePID[kPIDs]  = { "All", "Cpv", "Disp", "Both", "Cpv2", "Disp2", "Both2"               };
  UInt_t  PIDBIT[kPIDs]   = { 0x0,  BIT(0), BIT(1), BIT(0)|BIT(1), BIT(2), BIT(3), BIT(2)|BIT(3)   };

  // MC truth info
  TString pName;
  const Int_t nParticles = stack->GetNtrack();
  Int_t iMod = 0, jMod = 0;
  for(Int_t iParticle=0; iParticle<nParticles; iParticle++) {
    TParticle* pMC  = stack->Particle(iParticle);
    if (pMC->GetPdgCode() != 111) continue;   // Pi0
    pName = ClassifyMCPi0(iParticle, stack);

    //Pi0 Information
    Double_t pt     = pMC->Pt();
    Double_t r      = pMC->R();
    Double_t y      = pMC->Y();
    Double_t phi    = pMC->Phi();   while (phi<0.) phi += TMath::TwoPi(); while (phi>TMath::TwoPi()) phi -= TMath::TwoPi();

    ((TH2D*)listMC->FindObject("hMCInclusivePi0_PtRadius"))->Fill(pt, r);
    ((TH2D*)listMC->FindObject(Form("hMCInclusivePi0_PtRadius_cent%d", cent)))->Fill(pt, r);
    ((TH2D*)listMC->FindObject("hMCInclusivePi0_PtRapidity"))->Fill(pt, y);
    ((TH2D*)listMC->FindObject(Form("hMCInclusivePi0_PtRapidity_cent%d", cent)))->Fill(pt, y);

    ((TH2D*)listMC->FindObject(Form("hMC%s_PtRadius", pName.Data())))->Fill(pt, r);
    ((TH2D*)listMC->FindObject(Form("hMC%s_PtRadius_cent%d", pName.Data(), cent)))->Fill(pt, r);
    ((TH2D*)listMC->FindObject(Form("hMC%s_PtRapidity", pName.Data())))->Fill(pt, y);
    ((TH2D*)listMC->FindObject(Form("hMC%s_PtRapidity_cent%d", pName.Data(), cent)))->Fill(pt, y);

    if ((((phi>260.*TMath::DegToRad()) && (phi<280.*TMath::DegToRad())) || ((phi>300.*TMath::DegToRad()) && (phi<320.*TMath::DegToRad()))) && (TMath::Abs(y)<0.135)) {
      ((TH1D*)listMC->FindObject("hMCInclusivePi0_Pt_Pi0InAcc"))->Fill(pt);
      ((TH1D*)listMC->FindObject(Form("hMCInclusivePi0_Pt_Pi0InAcc_cent%d", cent)))->Fill(pt);

      ((TH1D*)listMC->FindObject(Form("hMC%s_Pt_Pi0InAcc", pName.Data())))->Fill(pt);
      ((TH1D*)listMC->FindObject(Form("hMC%s_Pt_Pi0InAcc_cent%d", pName.Data(), cent)))->Fill(pt);
    }

    TParticle *iGamma = 0x0, *jGamma = 0x0;
    if (pMC->GetNDaughters()==2) { // Do not account Dalitz decays
      iGamma = stack->Particle(pMC->GetFirstDaughter());
      jGamma = stack->Particle(pMC->GetLastDaughter());
      //Pi0's daughter particles in PHOS acceptance
      iMod = HitPHOSModule(iGamma, phosGeo);
      jMod = HitPHOSModule(jGamma, phosGeo);
      if (iMod!=-1 && iMod!=2 && jMod!=-1 && jMod!=2 && TMath::Abs(iMod-jMod)<2) { // !remove module 2
        ((TH1D*)listMC->FindObject("hMCInclusivePi0_Pt_GammaInAcc"))->Fill(pt);
        ((TH1D*)listMC->FindObject(Form("hMCInclusivePi0_Pt_GammaInAcc_cent%d", cent)))->Fill(pt);

        ((TH1D*)listMC->FindObject(Form("hMC%s_Pt_GammaInAcc", pName.Data())))->Fill(pt);
        ((TH1D*)listMC->FindObject(Form("hMC%s_Pt_GammaInAcc_cent%d", pName.Data(), cent)))->Fill(pt);
      }
    }
  }

  // Reconstruction info
  Int_t  entries  = caloClArr->GetEntriesFast();
  Int_t  iPi0Indx =-1, jPi0Indx =-1;
  UInt_t iPIDBit  = 0, jPIDBit  = 0;
  AliCaloClusterInfo *iCaloCluster = 0, *jCaloCluster = 0;
  TLorentzVector iLorentz, jLorentz;
  for(Int_t i=0; i<entries-1; i++) { // Loop calo cluster i
    iCaloCluster = (AliCaloClusterInfo*)caloClArr->At(i);
    if (!iCaloCluster->CheckIsClusterFromPi0(stack, iPi0Indx)) { iCaloCluster = 0; continue; }
    iMod     = iCaloCluster->GetModule();
    iPIDBit  = iCaloCluster->GetPIDBit();
    iLorentz = iCaloCluster->LorentzVector();

    for (Int_t j=i+1; j<entries; j++) { // Loop calo cluster j
      jCaloCluster = (AliCaloClusterInfo*)caloClArr->At(j);
      if (!jCaloCluster->CheckIsClusterFromPi0(stack, jPi0Indx)) { jCaloCluster = 0; continue; }
      if (jPi0Indx != iPi0Indx) { jCaloCluster = 0; continue; }  // coming from the same pi0
      pName    = ClassifyMCPi0(iPi0Indx, stack);
      jMod     = jCaloCluster->GetModule();
      jPIDBit  = jCaloCluster->GetPIDBit();
      jLorentz = jCaloCluster->LorentzVector();
      if (TMath::Abs(iMod-jMod)>1) { jCaloCluster = 0; continue; }

      Double_t pi0Pt      = (iLorentz+jLorentz).E();
      Double_t pi0InvMass = (iLorentz+jLorentz).M();

      for (Int_t iPID=0; iPID<kPIDs; iPID++) {
        if ((iPIDBit & jPIDBit & PIDBIT[iPID])==PIDBIT[iPID]) {
          ((TH2D*)listMC->FindObject(Form("hMCInclusivePi0_PtInvMass_%s_Reco", namePID[iPID].Data())))->Fill(pi0Pt, pi0InvMass);
          ((TH2D*)listMC->FindObject(Form("hMCInclusivePi0_PtInvMass_%s_Reco_cent%d", namePID[iPID].Data(), cent)))->Fill(pi0Pt, pi0InvMass);

          ((TH2D*)listMC->FindObject(Form("hMC%s_PtInvMass_%s_Reco", pName.Data(), namePID[iPID].Data())))->Fill(pi0Pt, pi0InvMass);
          ((TH2D*)listMC->FindObject(Form("hMC%s_PtInvMass_%s_Reco_cent%d", pName.Data(), namePID[iPID].Data(), cent)))->Fill(pi0Pt, pi0InvMass);
        }
      }
    } // End loop calo cluster j 
  } // End loop calo cluster i

  return;
}

//________________________________________________________________________
TString AliPHOSpPbPi0Header::ClassifyMCPi0(Int_t index, AliStack* const stack)
{

    TString  pName;
    TParticle* pMC  = stack->Particle(index);
    Int_t      mpdg = 0;

    if (pMC->GetFirstMother()>-1) {
      TParticle* pMother = stack->Particle(pMC->GetFirstMother());
      mpdg = TMath::Abs(pMother->GetPdgCode());
    }

    // Classify Pi0 sources
    if (index<stack->GetNprimary())                              pName = "PrimaryPi0";
    else if (mpdg==310)                                          pName = "2ndPi0FromK0s";
    else if (mpdg==321 || mpdg==211 || mpdg==2212 || mpdg==2112) pName = "2ndPi0FromMaterials";     // K+-, pi+-, p(pbar), n(nbar)
    else                                                         pName = "2ndPi0FromOtherDecays";   // all other secondary pi0s

    return pName;
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
