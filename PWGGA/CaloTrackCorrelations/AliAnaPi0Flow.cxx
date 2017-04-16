/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// --- ROOT system ---
#include "TH3.h"
#include "TH2F.h"
//#include "Riostream.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TDatabasePDG.h"
#include <TParticle.h>
#include <TRefArray.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TVirtualFFT.h>


//---- AliRoot system ----
#include "AliAnaPi0Flow.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliVEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliNeutralMesonSelection.h"
#include "AliMixedEvent.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisTaskSE.h"
#include "AliRunLoader.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliStack.h"
#include "AliEMCALGeometry.h"
#include "AliEPFlattener.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliKFParticle.h"
#include "AliAODEvent.h"
#include "AliVCluster.h"
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliOADBContainer.h"

// --- Detectors ---
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"

#include <iostream>

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnaPi0Flow);
/// \endcond

//
// Constructor
//
AliAnaPi0Flow::AliAnaPi0Flow() : AliAnaCaloTrackCorrBaseClass(),

isPhosCali(0), isCentFlat(0),
fInputEvent(0x0), fEventPlane(0x0), fCentrality(-999.), fRunNumber(0), fInternalRunNum(0), fFlatContainer(0x0),

fTPCFlat(0x0), fV0AFlat(0x0), fV0CFlat(0x0),
fEPTPC(-999.), fEPTPCResolution(0.),
fEPV0(-999.), fEPV0A(-999.), fEPV0C(-999.),
fEPV0AR(-999.), fEPV0CR(-999.), fEPV0R(-999.),
fEPV0AR4(-999.), fEPV0AR5(-999.), fEPV0AR6(-999.), fEPV0AR7(-999.),
fEPV0CR0(-999.), fEPV0CR1(-999.), fEPV0CR2(-999.), fEPV0CR3(-999.),

fHistStatCentrality(0), fHistStatCentralityCorrected(0), fHistStatRunNum(0),

fHistEPTPC(0), fHistEPTPCResolution(0),
fHistEPV0(0), fHistEPV0A(0), fHistEPV0C(0),
fHistEPV0AR(0), fHistEPV0CR(0), fHistEPV0R(0),
fHistEPV0AR4(0), fHistEPV0AR7(0), fHistEPV0CR0(0), fHistEPV0CR3(0),
fHistEPTPCFlatten(0), fHistEPV0AFlatten(0), fHistEPV0CFlatten(0),
fHistEPDiffV0A_V0CR0(0), fHistEPDiffV0A_V0CR3(0), fHistEPDiffV0CR0_V0CR3(0),
fHistEPDiffV0C_V0AR4(0), fHistEPDiffV0C_V0AR7(0), fHistEPDiffV0AR4_V0AR7(0), fHistEPDiffV0AR_V0CR(0),
fHistClusterEtaPhi(0), fHistClusterEN(0), fHistClusterEtN(0), fHistClusterEM02(0), fHistClusterEtM02(0),

fDataV0(0), fDataV0A(0), fDataV0C(0), fDataTPC(0)
{
  InitParameters();
}


AliAnaPi0Flow::~AliAnaPi0Flow()
{}

//
// Initialize parameters
//
void AliAnaPi0Flow::InitParameters()
{
  SetInputAODName("AliAODPWG4Particle");  
  AddToHistogramsName("AnaPi0Flow_");

  fFlatContainer = new AliOADBContainer("phosFlat");
  fFlatContainer->InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root", "phosFlat");
}

//
// Add cuts (needed by AliAnaCaloTrackCorrMaker)
//
TObjString* AliAnaPi0Flow::GetAnalysisCuts()
{
  TString parList ; //this will be list of parameters used for this analysis.

  //Get parameters set in base class.
  parList += GetBaseParametersList() ;

  return new TObjString(parList) ;
}

//
// Create output (needed by AliAnaCaloTrackCorrMaker) 
//
TList * AliAnaPi0Flow::GetCreateOutputObjects()
{
  TList * outputList = new TList() ;
  outputList->SetName(GetName());
  
  //
  // statistics
  //
  fHistStatRunNum = new TH1I("fHistStatRunNum", "", 45, 0, 45);
  outputList->Add(fHistStatRunNum);

  fHistStatCentrality = new TH1D("fHistStatCentrality", "", 200, 0, 100);
  outputList->Add(fHistStatCentrality);

  fHistStatCentralityCorrected = new TH1D("fHistStatCentralityCorrected", "", 200, 0, 100);
  outputList->Add(fHistStatCentralityCorrected);

  //
  // event plane
  //
  fHistEPTPC = new TH2F("fHistEPTPC","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPTPC);

  fHistEPTPCResolution = new TH2F("fHistEPTPCResolution","",100,0,100,100,0.0,1.0);
  outputList->Add(fHistEPTPCResolution); 

  fHistEPV0 = new TH2F("fHistEPV0","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0);

  fHistEPV0A = new TH2F("fHistEPV0A","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0A);

  fHistEPV0C = new TH2F("fHistEPV0C","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0C);

  fHistEPV0AR = new TH2F("fHistEPV0AR","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0AR);

  fHistEPV0CR = new TH2F("fHistEPV0CR","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0CR);

  fHistEPV0R = new TH2F("fHistEPV0R","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0R);

  fHistEPV0AR4 = new TH2F("fHistEPV0AR4","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0AR4);

  fHistEPV0AR7 = new TH2F("fHistEPV0AR7","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0AR7);

  fHistEPV0CR0 = new TH2F("fHistEPV0CR0","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0CR0);

  fHistEPV0CR3 = new TH2F("fHistEPV0CR3","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0CR3);

  fHistEPV0AFlatten = new TH2F("fHistEPV0AFlatten","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0AFlatten);

  fHistEPV0CFlatten = new TH2F("fHistEPV0CFlatten","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPV0CFlatten);

  fHistEPTPCFlatten = new TH2F("fHistEPTPFlatten","",100,0,100,100,0.0,TMath::Pi());
  outputList->Add(fHistEPTPCFlatten);

  fHistEPDiffV0A_V0CR0 = new TH2F("fHistEPDiffV0A_V0CR0","",100,0,100,100,-1.0,1.0);
  outputList->Add(fHistEPDiffV0A_V0CR0);

  fHistEPDiffV0A_V0CR3 = new TH2F("fHistEPDiffV0A_V0CR3","",100,0,100,100,-1.0,1.0);
  outputList->Add(fHistEPDiffV0A_V0CR3);

  fHistEPDiffV0CR0_V0CR3 = new TH2F("fHistEPDiffV0CR0_V0CR3","",100,0,100,100,-1.0,1.0);
  outputList->Add(fHistEPDiffV0CR0_V0CR3);

  fHistEPDiffV0C_V0AR4 = new TH2F("fHistEPDiffV0C_V0AR4","",100,0,100,100,-1.0,1.0);
  outputList->Add(fHistEPDiffV0C_V0AR4);

  fHistEPDiffV0C_V0AR7 = new TH2F("fHistEPDiffV0C_V0AR7","",100,0,100,100,-1.0,1.0);
  outputList->Add(fHistEPDiffV0C_V0AR7);

  fHistEPDiffV0AR4_V0AR7 = new TH2F("fHistEPDiffV0AR4_V0AR7","",100,0,100,100,-1.0,1.0);
  outputList->Add(fHistEPDiffV0AR4_V0AR7);

  fHistEPDiffV0AR_V0CR = new TH2F("fHistEPDiffV0AR_V0CR","",100,0,100,100,-1.0,1.0);
  outputList->Add(fHistEPDiffV0AR_V0CR);

  fHistClusterEtaPhi = new TH2F("fHistClusterEtaPhi","Cluster #eta vs #phi",100,-1.0,1.0,100,0.0,6.29);
  fHistClusterEtaPhi->GetXaxis()->SetTitle("#eta"); fHistClusterEtaPhi->GetYaxis()->SetTitle("#phi [rad]");
  outputList->Add(fHistClusterEtaPhi);

  fHistClusterEN = new TH2F("fHistClusterEN","N vs E",100,0.0,50.0,30,0.0,30.0);
  fHistClusterEN->GetYaxis()->SetTitle("N"); fHistClusterEN->GetXaxis()->SetTitle("E [GeV]");
  outputList->Add(fHistClusterEN);

  fHistClusterEtN = new TH2F("fHistClusterEtN","N vs Cluster E_{T}",100,0.0,50.0,30,0.0,30.0);
  fHistClusterEtN->GetYaxis()->SetTitle("N"); fHistClusterEtN->GetXaxis()->SetTitle("E_{T} [GeV]");
  outputList->Add(fHistClusterEtN);

  fHistClusterEM02 = new TH2F("fHistClusterEM02","Cluster E vs M02",500,0.0,50.0,100,0.0,10.0);
  fHistClusterEM02->GetYaxis()->SetTitle("M02"); fHistClusterEM02->GetXaxis()->SetTitle("E [GeV]");
  outputList->Add(fHistClusterEM02);

  fHistClusterEtM02 = new TH2F("fHistClusterEtM02","Cluster E_{T} vs M02",100,0.0,50.0,100,0.0,10.0);
  fHistClusterEtM02->GetYaxis()->SetTitle("M02"); fHistClusterEtM02->GetXaxis()->SetTitle("E_{T} [GeV]");
  outputList->Add(fHistClusterEtM02);

  //                   E    Et    M02   Cent     DeltaPhi      Cos[2*DeltaPhi]
  Int_t    bins[6] = {500,  500,  500,  100,          100,     100};
  Double_t min[6]  = {0,    0.0,  0.0,    0,          0.0,    -1.0};
  Double_t max[6]  = {50.0, 50.0, 10.0, 100,  TMath::Pi(),     1.0};
   
  fDataV0 = new THnSparseF("fDataV0","",6,bins,min,max);
  fDataV0->GetAxis(0)->SetTitle("E"); fDataV0->GetAxis(1)->SetTitle("E_{T} [GeV]"); fDataV0->GetAxis(2)->SetTitle("M02");
  fDataV0->GetAxis(3)->SetTitle("Centrality"); fDataV0->GetAxis(4)->SetTitle("#Delta(#phi) [rad]"); fDataV0->GetAxis(5)->SetTitle("Cos[2*#Delta(#phi)]");  
  outputList->Add(fDataV0);
   
  fDataV0A = new THnSparseF("fDataV0A","",6,bins,min,max);
  fDataV0A->GetAxis(0)->SetTitle("E"); fDataV0A->GetAxis(1)->SetTitle("E_{T} [GeV]"); fDataV0A->GetAxis(2)->SetTitle("M02"); 
  fDataV0A->GetAxis(3)->SetTitle("Centrality"); fDataV0A->GetAxis(4)->SetTitle("#Delta(#phi) [rad]"); fDataV0A->GetAxis(5)->SetTitle("Cos[2*#Delta(#phi)]");  
  outputList->Add(fDataV0A);
   
  fDataV0C = new THnSparseF("fDataV0C","",6,bins,min,max);
  fDataV0C->GetAxis(0)->SetTitle("E"); fDataV0C->GetAxis(1)->SetTitle("E_{T} [GeV]"); fDataV0C->GetAxis(2)->SetTitle("M02");
  fDataV0C->GetAxis(3)->SetTitle("Centrality"); fDataV0C->GetAxis(4)->SetTitle("#Delta(#phi) [rad]"); fDataV0C->GetAxis(5)->SetTitle("Cos[2*#Delta(#phi)]");  
  outputList->Add(fDataV0C);
   
  fDataTPC = new THnSparseF("fDataTPC","",6,bins,min,max);
  fDataTPC->GetAxis(0)->SetTitle("E"); fDataTPC->GetAxis(1)->SetTitle("E_{T} [GeV]"); fDataTPC->GetAxis(2)->SetTitle("M02");
  fDataTPC->GetAxis(3)->SetTitle("Centrality"); fDataTPC->GetAxis(4)->SetTitle("#Delta(#phi) [rad]"); fDataTPC->GetAxis(5)->SetTitle("Cos[2*#Delta(#phi)]");  
  outputList->Add(fDataTPC);

  return outputList;
}

//
// Print (needed by AliAnaCaloTrackCorrMaker) 
//
void AliAnaPi0Flow::Print(const Option_t * /*opt*/) const
{
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");  
  printf("------------------------------------------------------\n") ;
}

//
// Main method (needed by AliAnaCaloTrackCorrMaker)
//
void AliAnaPi0Flow::MakeAnalysisFillHistograms()
{  
  if (!GetInputAODBranch()) {
    AliFatal(Form("ERROR: No input in AOD with name branch < %s >", GetInputAODName().Data()));
    return;
  }

  fInputEvent = GetReader()->GetInputEvent();
  fEventPlane = fInputEvent->GetEventplane();
  fCentrality = GetEventCentrality();
  fHistStatCentrality->Fill(fCentrality);
  if (isCentFlat && !IsCentAccepted()) return;
  fHistStatCentralityCorrected->Fill(fCentrality);

  Int_t fRunNumber = fInputEvent->GetRunNumber();
  fInternalRunNum = GetInternalRunNum(fRunNumber);
  fHistStatRunNum->Fill(fInternalRunNum);

  TObjArray *maps = (TObjArray*)fFlatContainer->GetObject(fRunNumber,"phosFlat");
  fTPCFlat = (AliEPFlattener*)maps->At(0);
  fV0AFlat = (AliEPFlattener*)maps->At(1);
  fV0CFlat = (AliEPFlattener*)maps->At(2);
  if (fEventPlane) {
    GetVZEROEventPlane(isPhosCali);
  } else {
    AliFatal("ERROR: Could not retrieve the Event Plane!");
    return;
  }

  //
  // Get me some pions
  //
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if (naod == 0) {
    AliDebug(1, "No cluster found in this event!");
    return;
  }

  for (Int_t i = 0; i < naod; ++i) {
    AliAODPWG4Particle* c = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i));

    Double_t E = c->E();
    Double_t eta = c->Eta();
    Double_t Et = c->Pt();
    Double_t phi = c->Phi();
    Double_t M02 = c->GetM02();
    Double_t N = c->GetNCells();

    Double_t dphiV0  = TVector2::Phi_0_2pi(phi-fEPV0R); if(dphiV0 >TMath::Pi()) {dphiV0 -=TMath::Pi();}
    Double_t dphiV0A = TVector2::Phi_0_2pi(phi-fEPV0A); if(dphiV0A>TMath::Pi()) {dphiV0A-=TMath::Pi();}
    Double_t dphiV0C = TVector2::Phi_0_2pi(phi-fEPV0C); if(dphiV0C>TMath::Pi()) {dphiV0C-=TMath::Pi();}
    Double_t dphiTPC = TVector2::Phi_0_2pi(phi-fEPTPC); if(dphiTPC>TMath::Pi()) {dphiTPC-=TMath::Pi();}

    //
    // Fill flow data
    //
    Double_t dataV0[6];
    dataV0[0] = E;
    dataV0[1] = Et;
    dataV0[2] = M02;
    dataV0[3] = fCentrality;
    dataV0[4] = dphiV0;
    dataV0[5] = TMath::Cos(2.0*(dphiV0));
    fDataV0->Fill(dataV0);

    Double_t dataV0A[6];
    dataV0A[0] = E;
    dataV0A[1] = Et;
    dataV0A[2] = M02;
    dataV0A[3] = fCentrality;
    dataV0A[4] = dphiV0A;
    dataV0A[5] = TMath::Cos(2.0*(dphiV0A));
    fDataV0A->Fill(dataV0A);

    Double_t dataV0C[6];
    dataV0C[0] = E;
    dataV0C[1] = Et;
    dataV0C[2] = M02;
    dataV0C[3] = fCentrality;
    dataV0C[4] = dphiV0C;
    dataV0C[5] = TMath::Cos(2.0*(dphiV0C));
    fDataV0C->Fill(dataV0C);

    Double_t dataTPC[6];
    dataTPC[0] = E;
    dataTPC[1] = Et;
    dataTPC[2] = M02;
    dataTPC[3] = fCentrality;
    dataTPC[4] = dphiTPC;
    dataTPC[5] = TMath::Cos(2.0*(dphiTPC));
    fDataTPC->Fill(dataTPC);

    //
    // Fill some QA hists
    //
    fHistClusterEtaPhi->Fill(eta,phi);
    fHistClusterEN->Fill(E,N);
    fHistClusterEtN->Fill(Et,N);
    fHistClusterEM02->Fill(E,M02);
    fHistClusterEtM02->Fill(Et,M02);
  }

  AliDebug(1, "End fill histograms");
}

//
// Flatten centrality
//
Bool_t AliAnaPi0Flow::IsCentAccepted()
{
  if (fCentrality<=10) {   // 0-10%
    Bool_t isGoodCent = kFALSE;
    TRandom3 *rndm = new TRandom3(0);
    Double_t Nrndm = rndm->Uniform(0.,1.);
    Double_t rejectR[10] = {0.864235, 0.855794, 0.838345, 0.837401, 0.838991, 
                            0.859132, 0.852702, 0.859955, 0.880277, 1}; 

    if (fCentrality<=1) {
      if (Nrndm < rejectR[0]) isGoodCent = kTRUE;
    } else if (1<fCentrality && fCentrality<=2) {
      if (Nrndm < rejectR[1]) isGoodCent = kTRUE;
    } else if (2<fCentrality && fCentrality<=3) {
      if (Nrndm < rejectR[2]) isGoodCent = kTRUE;
    } else if (3<fCentrality && fCentrality<=4) {
      if (Nrndm < rejectR[3]) isGoodCent = kTRUE;
    } else if (4<fCentrality && fCentrality<=5) {
      if (Nrndm < rejectR[4]) isGoodCent = kTRUE;
    } else if (5<fCentrality && fCentrality<=6) {
      if (Nrndm < rejectR[5]) isGoodCent = kTRUE;
    } else if (6<fCentrality && fCentrality<=7) {
      if (Nrndm < rejectR[6]) isGoodCent = kTRUE;
    } else if (7<fCentrality && fCentrality<=8) {
      if (Nrndm < rejectR[7]) isGoodCent = kTRUE;
    } else if (8<fCentrality && fCentrality<=9) {
      if (Nrndm < rejectR[8]) isGoodCent = kTRUE;
    } else if (9<fCentrality && fCentrality<=10) {
      if (Nrndm < rejectR[9]) isGoodCent = kTRUE;
    }

    delete rndm; rndm = NULL;
    return isGoodCent;
  } else if (10<fCentrality && fCentrality<=50) {  // 10-50%
    TString centfired;
    centfired = fInputEvent->GetFiredTriggerClasses();
    if (!centfired.Contains("CVLN_B2-B-NOPF-ALLNOTRD") && 
        !centfired.Contains("CVLN_R1-B-NOPF-ALLNOTRD") && 
        !centfired.Contains("CSEMI_R1-B-NOPF-ALLNOTRD")) return kFALSE;
    else return kTRUE;
  } else return kTRUE; // other%

  return kFALSE;
}

//
// Internal run number
//
Int_t AliAnaPi0Flow::GetInternalRunNum(Int_t runnumber)
{ // LHC11h_AOD_145, train runlist listPCMgood+VaryingVoltages (no phi cut)
  Int_t stdRunNum[45] = {167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 
                         168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168984, 
                         168988, 168992, 169035, 169091, 169094, 169138, 169143, 169144, 169145, 169148, 
                         169156, 169160, 169167, 169238, 169411, 169415, 169417, 169835, 169837, 169846, 
                         169855, 169858, 169859, 169923, 169965};

  for (int i = 0; i < 45; ++i) {
    if (i%10==0) fHistStatRunNum->GetXaxis()->SetBinLabel(i+1, Form("%i",stdRunNum[i]));    
    if (runnumber==stdRunNum[i]) return i;
  }

  return -1;
}


//
// Get different event planes, flatten them, and fill some QA plots
//
void AliAnaPi0Flow::GetVZEROEventPlane(Bool_t flattenEP)
{ // Calculate the V0 Event Plane
  if (fEventPlane->GetQVector()) { 
    // fEPTPC = TVector2::Phi_0_2pi(fEventPlane->GetQVector()->Phi())/2.0; //if(fEPTPC>TMath::Pi()) {fEPTPC-=TMath::Pi();} 
    fEPTPC = fEventPlane->GetEventplane("Q");
  } else fEPTPC = -999;

  if (fEventPlane->GetQsub1() && fEventPlane->GetQsub2()) {
    // fEPTPCResolution = TMath::Cos(2.0*(fEventPlane->GetQsub1()->Phi()/2.0-fEventPlane->GetQsub2()->Phi()/2.0));
    fEPTPCResolution = TMath::Cos(2.*(fEventPlane->GetQsubRes()));
  } else fEPTPCResolution = -1;

  fEPV0  = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0",  fInputEvent)); if (fEPV0>TMath::Pi()) fEPV0-=TMath::Pi();
  fEPV0A = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0A", fInputEvent)); if (fEPV0A>TMath::Pi()) fEPV0A-=TMath::Pi();
  fEPV0C = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0C", fInputEvent)); if (fEPV0C>TMath::Pi()) fEPV0C-=TMath::Pi();

  Double_t qx=0, qy=0, qxr=0, qyr=0;
  fEPV0AR = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 4, 5, 2, qxr, qyr)); if (fEPV0AR>TMath::Pi()) fEPV0AR-=TMath::Pi();
  fEPV0CR = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 2, 3, 2, qx,  qy)); if (fEPV0CR>TMath::Pi()) fEPV0CR-=TMath::Pi();
  qxr += qx; qyr += qy;
  fEPV0R   = TVector2::Phi_0_2pi(TMath::ATan2(qyr,qxr))/2.0; // equals to ring 2-5
  fEPV0AR4 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 4, 2, qx, qy)); if (fEPV0AR4>TMath::Pi()) fEPV0AR4-=TMath::Pi();
  fEPV0AR5 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 5, 2, qx, qy)); if (fEPV0AR5>TMath::Pi()) fEPV0AR5-=TMath::Pi();
  fEPV0AR6 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 6, 2, qx, qy)); if (fEPV0AR6>TMath::Pi()) fEPV0AR6-=TMath::Pi();
  fEPV0AR7 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 7, 2, qx, qy)); if (fEPV0AR7>TMath::Pi()) fEPV0AR7-=TMath::Pi();
  fEPV0CR0 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 0, 2, qx, qy)); if (fEPV0CR0>TMath::Pi()) fEPV0CR0-=TMath::Pi();
  fEPV0CR1 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 1, 2, qx, qy)); if (fEPV0CR1>TMath::Pi()) fEPV0CR1-=TMath::Pi();
  fEPV0CR2 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 2, 2, qx, qy)); if (fEPV0CR2>TMath::Pi()) fEPV0CR2-=TMath::Pi();
  fEPV0CR3 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fInputEvent, 3, 2, qx, qy)); if (fEPV0CR3>TMath::Pi()) fEPV0CR3-=TMath::Pi();

  fHistEPTPC->Fill(fCentrality,  fEPTPC); 
  if (fEPTPCResolution!=-1) fHistEPTPCResolution->Fill(fCentrality, fEPTPCResolution);
  fHistEPV0->Fill(fCentrality, fEPV0);
  fHistEPV0A->Fill(fCentrality, fEPV0A);
  fHistEPV0C->Fill(fCentrality, fEPV0C);
  fHistEPV0AR->Fill(fCentrality, fEPV0AR);
  fHistEPV0CR->Fill(fCentrality, fEPV0CR);
  fHistEPV0R->Fill(fCentrality, fEPV0R);
  fHistEPV0AR4->Fill(fCentrality, fEPV0AR4);
  fHistEPV0AR7->Fill(fCentrality, fEPV0AR7);
  fHistEPV0CR0->Fill(fCentrality, fEPV0CR0);
  fHistEPV0CR3->Fill(fCentrality, fEPV0CR3);

  if (flattenEP) {
    fEPV0A = ApplyFlatteningV0A(fEPV0A, fCentrality);
    fEPV0C = ApplyFlatteningV0C(fEPV0C, fCentrality);
    if (fEPTPC != -999.) fEPTPC = ApplyFlatteningTPC(fEPTPC, fCentrality);
    while (fEPV0A<0.) fEPV0A+=TMath::Pi(); while (fEPV0A>TMath::Pi()) fEPV0A-=TMath::Pi();
    while (fEPV0C<0.) fEPV0C+=TMath::Pi(); while (fEPV0C>TMath::Pi()) fEPV0C-=TMath::Pi();
    while (fEPTPC<0.) fEPTPC+=TMath::Pi(); while (fEPTPC>TMath::Pi()) fEPTPC-=TMath::Pi();

    fHistEPTPCFlatten->Fill(fCentrality, fEPTPC);
    fHistEPV0AFlatten->Fill(fCentrality, fEPV0A);
    fHistEPV0CFlatten->Fill(fCentrality, fEPV0C);
  }

  fHistEPDiffV0A_V0CR0->Fill(fCentrality, TMath::Cos(2.0*(fEPV0A - fEPV0CR0)));
  fHistEPDiffV0A_V0CR3->Fill(fCentrality, TMath::Cos(2.0*(fEPV0A - fEPV0CR3)));
  fHistEPDiffV0CR0_V0CR3->Fill(fCentrality, TMath::Cos(2.0*(fEPV0CR0 - fEPV0CR3)));
  fHistEPDiffV0C_V0AR4->Fill(fCentrality, TMath::Cos(2.0*(fEPV0C - fEPV0AR4)));
  fHistEPDiffV0C_V0AR7->Fill(fCentrality, TMath::Cos(2.0*(fEPV0C - fEPV0AR7)));
  fHistEPDiffV0AR4_V0AR7->Fill(fCentrality, TMath::Cos(2.0*(fEPV0AR4 - fEPV0AR7)));   
  fHistEPDiffV0AR_V0CR->Fill(fCentrality, TMath::Cos(2.0*(fEPV0AR - fEPV0CR)));
}

Double_t AliAnaPi0Flow::ApplyFlatteningV0A(Double_t phi, Double_t c)
{
  if (fV0AFlat) return fV0AFlat->MakeFlat(phi,c);
  return phi;
}

Double_t AliAnaPi0Flow::ApplyFlatteningV0C(Double_t phi, Double_t c)
{ 
  if (fV0CFlat) return fV0CFlat->MakeFlat(phi,c);
  return phi;
}

Double_t AliAnaPi0Flow::ApplyFlatteningTPC(Double_t phi, Double_t c)
{
  if (fTPCFlat) return fTPCFlat->MakeFlat(phi,c);
  return phi;
}

