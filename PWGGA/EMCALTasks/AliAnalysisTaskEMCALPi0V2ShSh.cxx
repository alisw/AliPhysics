// $Id: AliAnalysisTaskEMCALPi0V2ShSh.cxx$

#include "AliAnalysisTaskEMCALPi0V2ShSh.h"

//Root include files 
//#include <Riostream.h>
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

//AliRoot include files 
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



ClassImp(AliAnalysisTaskEMCALPi0V2ShSh)

//________________________________________________________________________
AliAnalysisTaskEMCALPi0V2ShSh::AliAnalysisTaskEMCALPi0V2ShSh() : 
  AliAnalysisTaskSE("default_name"), 
  fESD(0), fAOD(0),
  fCentrality(-999.), fEventPlane(0), fFlatContainer(0),
  fRunNumber(0), fInternalRunNum(0), 
  fTPCFlat(0x0), fV0AFlat(0x0), fV0CFlat(0x0),
  fEMCALClustersListName(""), fESDClusters(0), fAODClusters(0),
  fESDCells(0), fAODCells(0),
  fGeom(0), fGeoName("EMCAL_COMPLETEV1"), fOADBContainer(0),
  fOutputList(0),
  fEPTPC(-999.), fEPTPCResolution(0.),
  fEPV0(-999.), fEPV0A(-999.), fEPV0C(-999.),
  fEPV0AR(-999.), fEPV0CR(-999.), fEPV0R(-999.),
  fEPV0AR4(-999.), fEPV0AR5(-999.), fEPV0AR6(-999.), fEPV0AR7(-999.),
  fEPV0CR0(-999.), fEPV0CR1(-999.), fEPV0CR2(-999.), fEPV0CR3(-999.),
  fHistEPV0MB(0), fHistEPV0RMB(0), fHistEPV0AMB(0), fHistEPV0CMB(0), fHistEPTPCMB(0),
  fHistEPTPC(0), fHistEPTPCResolution(0),
  fHistEPV0(0), fHistEPV0A(0), fHistEPV0C(0),
  fHistEPV0AR(0), fHistEPV0CR(0), fHistEPV0R(0),
  fHistEPV0AR4(0), fHistEPV0AR7(0), fHistEPV0CR0(0), fHistEPV0CR3(0),
  fHistDiffV0A_V0CR0(0), fHistDiffV0A_V0CR3(0), fHistDiffV0CR0_V0CR3(0),
  fHistDiffV0C_V0AR4(0), fHistDiffV0C_V0AR7(0), fHistDiffV0AR4_V0AR7(0), fHistDiffV0AR_V0CR(0),
  fHistClusterEta(0), fHistClusterPhi(0), fHistClusterPhiEta(0),
  fHistClusterE(0), fHistClusterEt(0),
  fHistClusterN(0), fHistClusterM02(0),
  fHistClusterEN(0), fHistClusterEtN(0), 
  fHistClusterEM02(0), fHistClusterEtM02(0),
  fHistClusterdphiV0(0),
  fHistTrackPt(0), fHistTrackEta(0), fHistTrackPhi(0), fHistTrackPhiEta(0),
  fClusterV0(0), fClusterV0A(0), fClusterV0C(0), fClusterTPC(0),
  fHistEvtStat(0), fHistClusterStat(0), fHistRunNumStat(0), fHistCentrality(0), fHistCentralityCorrected(0),
  fHistEPTPCFlatten(0), fHistEPV0AFlatten(0), fHistEPV0CFlatten(0), 
  fHistRBRCosV0A(0), fHistRBRSinV0A(0), fHistRBRCosV0C(0), fHistRBRSinV0C(0), fHistRBRCosTPC(0), fHistRBRSinTPC(0)
{ // Default constructor
  for (Int_t i = 0; i < 12; ++i) fGeomMatrix[i] = 0;
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0V2ShSh::AliAnalysisTaskEMCALPi0V2ShSh(const char *name) : 
  AliAnalysisTaskSE(name), 
  fESD(0), fAOD(0),
  fCentrality(-999.), fEventPlane(0), fFlatContainer(0),
  fRunNumber(0), fInternalRunNum(0), 
  fTPCFlat(0x0), fV0AFlat(0x0), fV0CFlat(0x0),
  fEMCALClustersListName(""), fESDClusters(0), fAODClusters(0),
  fESDCells(0), fAODCells(0),
  fGeom(0), fGeoName("EMCAL_COMPLETEV1"), fOADBContainer(0),
  fOutputList(0),
  fEPTPC(-999.), fEPTPCResolution(0.),
  fEPV0(-999.), fEPV0A(-999.), fEPV0C(-999.),
  fEPV0AR(-999.), fEPV0CR(-999.), fEPV0R(-999.),
  fEPV0AR4(-999.), fEPV0AR5(-999.), fEPV0AR6(-999.), fEPV0AR7(-999.),
  fEPV0CR0(-999.), fEPV0CR1(-999.), fEPV0CR2(-999.), fEPV0CR3(-999.),
  fHistEPV0MB(0), fHistEPV0RMB(0), fHistEPV0AMB(0), fHistEPV0CMB(0), fHistEPTPCMB(0),
  fHistEPTPC(0), fHistEPTPCResolution(0),
  fHistEPV0(0), fHistEPV0A(0), fHistEPV0C(0),
  fHistEPV0AR(0), fHistEPV0CR(0), fHistEPV0R(0),
  fHistEPV0AR4(0), fHistEPV0AR7(0), fHistEPV0CR0(0), fHistEPV0CR3(0),
  fHistDiffV0A_V0CR0(0), fHistDiffV0A_V0CR3(0), fHistDiffV0CR0_V0CR3(0),
  fHistDiffV0C_V0AR4(0), fHistDiffV0C_V0AR7(0), fHistDiffV0AR4_V0AR7(0), fHistDiffV0AR_V0CR(0),
  fHistClusterEta(0), fHistClusterPhi(0), fHistClusterPhiEta(0),
  fHistClusterE(0), fHistClusterEt(0),
  fHistClusterN(0), fHistClusterM02(0),
  fHistClusterEN(0), fHistClusterEtN(0), 
  fHistClusterEM02(0), fHistClusterEtM02(0),
  fHistClusterdphiV0(0),
  fHistTrackPt(0), fHistTrackEta(0), fHistTrackPhi(0), fHistTrackPhiEta(0),
  fClusterV0(0), fClusterV0A(0), fClusterV0C(0), fClusterTPC(0),
  fHistEvtStat(0), fHistClusterStat(0), fHistRunNumStat(0), fHistCentrality(0), fHistCentralityCorrected(0),
  fHistEPTPCFlatten(0), fHistEPV0AFlatten(0), fHistEPV0CFlatten(0), 
  fHistRBRCosV0A(0), fHistRBRSinV0A(0), fHistRBRCosV0C(0), fHistRBRSinV0C(0), fHistRBRCosTPC(0), fHistRBRSinTPC(0)
{ // Constructor
  for (Int_t i = 0; i < 12; ++i) fGeomMatrix[i] = 0;
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::UserCreateOutputObjects()
{ // Create histograms, called once
    
  fESDClusters = new TObjArray();
  fAODClusters = new TObjArray();
  fGeom = AliEMCALGeometry::GetInstance(fGeoName.Data());
  fOADBContainer = new AliOADBContainer("AliEMCALgeo");
  fOADBContainer->InitFromFile(Form("$ALICE_PHYSICS/OADB/EMCAL/EMCALlocal2master.root"), "AliEMCALgeo");
  fFlatContainer = new AliOADBContainer("phosFlat");
  fFlatContainer->InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root", "phosFlat");

  fOutputList = new TList();
  fOutputList->SetOwner();// Container cleans up all histos (avoids leaks in merging)   

  fHistEvtStat = new TH1I("fHistEvtStat","",10,0,10);
  fHistEvtStat->GetXaxis()->SetBinLabel(1, "All");
  fHistEvtStat->GetXaxis()->SetBinLabel(2, "ESD");
  fHistEvtStat->GetXaxis()->SetBinLabel(3, "AOD");
  fHistEvtStat->GetXaxis()->SetBinLabel(4, "!(ESD/AOD)");
  fHistEvtStat->GetXaxis()->SetBinLabel(5, "Cent");
  fHistEvtStat->GetXaxis()->SetBinLabel(6, "VZERO");
  fHistEvtStat->GetXaxis()->SetBinLabel(7, "Clust");
  fOutputList->Add(fHistEvtStat);

  fHistClusterStat = new TH1I("fHistClusterStat","",10,0,10);
  fHistClusterStat->GetXaxis()->SetBinLabel(1, "All Clust.");
  fHistClusterStat->GetXaxis()->SetBinLabel(2, "EMCal Clust.");
  fHistClusterStat->GetXaxis()->SetBinLabel(3, "Good V1 Clust.");
  fOutputList->Add(fHistClusterStat);

  fHistRunNumStat = new TH1I("fHistRunNumStat", "", 45, 0, 45);
  fHistRunNumStat->GetXaxis()->SetBinLabel(1, "167987");
  fOutputList->Add(fHistRunNumStat);

  fHistCentrality = new TH1D("fHistCentrality", "", 200, 0, 100);
  fOutputList->Add(fHistCentrality);

  fHistCentralityCorrected = new TH1D("fHistCentralityCorrected", "", 200, 0, 100);
  fOutputList->Add(fHistCentralityCorrected);
  
  fHistEPV0MB = new TH1F("fHistEPV0MB","",100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0MB);

  fHistEPV0RMB = new TH1F("fHistEPV0RMB","",100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0RMB);

  fHistEPV0AMB = new TH1F("fHistEPV0AMB","",100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0AMB);

  fHistEPV0CMB = new TH1F("fHistEPV0CMB","",100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0CMB);

  fHistEPTPCMB = new TH1F("fHistEPTPCMB","",100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPTPCMB);

  fHistEPTPC = new TH2F("fHistEPTPC","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPTPC);

  fHistEPTPCResolution = new TH2F("fHistEPTPCResolution","",100,0,100,100,0.0,1.0);
  fOutputList->Add(fHistEPTPCResolution); 

  fHistEPV0 = new TH2F("fHistEPV0","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0);

  fHistEPV0A = new TH2F("fHistEPV0A","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0A);

  fHistEPV0C = new TH2F("fHistEPV0C","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0C);

  fHistEPV0AR = new TH2F("fHistEPV0AR","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0AR);

  fHistEPV0CR = new TH2F("fHistEPV0CR","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0CR);

  fHistEPV0R = new TH2F("fHistEPV0R","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0R);

  fHistEPV0AR4 = new TH2F("fHistEPV0AR4","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0AR4);

  fHistEPV0AR7 = new TH2F("fHistEPV0AR7","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0AR7);

  fHistEPV0CR0 = new TH2F("fHistEPV0CR0","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0CR0);

  fHistEPV0CR3 = new TH2F("fHistEPV0CR3","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0CR3);

  fHistEPV0AFlatten = new TH2F("fHistEPV0AFlatten","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0AFlatten);

  fHistEPV0CFlatten = new TH2F("fHistEPV0CFlatten","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPV0CFlatten);

  fHistEPTPCFlatten = new TH2F("fHistEPTPFlatten","",100,0,100,100,0.0,TMath::Pi());
  fOutputList->Add(fHistEPTPCFlatten);

  fHistDiffV0A_V0CR0 = new TH2F("fHistDiffV0A_V0CR0","",100,0,100,100,-1.0,1.0);
  fOutputList->Add(fHistDiffV0A_V0CR0);

  fHistDiffV0A_V0CR3 = new TH2F("fHistDiffV0A_V0CR3","",100,0,100,100,-1.0,1.0);
  fOutputList->Add(fHistDiffV0A_V0CR3);

  fHistDiffV0CR0_V0CR3 = new TH2F("fHistDiffV0CR0_V0CR3","",100,0,100,100,-1.0,1.0);
  fOutputList->Add(fHistDiffV0CR0_V0CR3);

  fHistDiffV0C_V0AR4 = new TH2F("fHistDiffV0C_V0AR4","",100,0,100,100,-1.0,1.0);
  fOutputList->Add(fHistDiffV0C_V0AR4);

  fHistDiffV0C_V0AR7 = new TH2F("fHistDiffV0C_V0AR7","",100,0,100,100,-1.0,1.0);
  fOutputList->Add(fHistDiffV0C_V0AR7);

  fHistDiffV0AR4_V0AR7 = new TH2F("fHistDiffV0AR4_V0AR7","",100,0,100,100,-1.0,1.0);
  fOutputList->Add(fHistDiffV0AR4_V0AR7);

  fHistDiffV0AR_V0CR = new TH2F("fHistDiffV0AR_V0CR","",100,0,100,100,-1.0,1.0);
  fOutputList->Add(fHistDiffV0AR_V0CR);

  fHistRBRCosV0A = new TProfile("fHistRBRCosV0A","",45,0,45,-1.,1.);
  fOutputList->Add(fHistRBRCosV0A);

  fHistRBRSinV0A = new TProfile("fHistRBRSinV0A","",45,0,45,-1.,1.);
  fOutputList->Add(fHistRBRSinV0A);

  fHistRBRCosV0C = new TProfile("fHistRBRCosV0C","",45,0,45,-1.,1.);
  fOutputList->Add(fHistRBRCosV0C);

  fHistRBRSinV0C = new TProfile("fHistRBRSinV0C","",45,0,45,-1.,1.);
  fOutputList->Add(fHistRBRSinV0C);

  fHistRBRCosTPC = new TProfile("fHistRBRCosTPC","",45,0,45,-1.,1.);
  fOutputList->Add(fHistRBRCosTPC);

  fHistRBRSinTPC = new TProfile("fHistRBRSinTPC","",45,0,45,-1.,1.);
  fOutputList->Add(fHistRBRSinTPC);
  
  fHistClusterEta = new TH1F("fHistClusterEta","Cluster Pseudorapidity Distribution",100,-1.0,1.0);
  fHistClusterEta->GetYaxis()->SetTitle("Entries"); fHistClusterEta->GetXaxis()->SetTitle("#eta");
  fOutputList->Add(fHistClusterEta);

  fHistClusterPhi = new TH1F("fHistClusterPhi","Cluster #phi Distribution",100,0.0,6.29);
  fHistClusterPhi->GetYaxis()->SetTitle("Entries"); fHistClusterPhi->GetXaxis()->SetTitle("#phi [rad]");
  fOutputList->Add(fHistClusterPhi);

  fHistClusterE = new TH1F("fHistClusterE","Cluster Energy Distribution",100,0.0,50.0);
  fHistClusterE->GetYaxis()->SetTitle("Entries"); fHistClusterE->GetXaxis()->SetTitle("Energy [GeV]");
  fOutputList->Add(fHistClusterE);
  
  fHistClusterEt = new TH1F("fHistClusterEt","Cluster Transverse Energy Distribution",100,0.0,50.0);
  fHistClusterEt->GetYaxis()->SetTitle("Entries"); fHistClusterEt->GetXaxis()->SetTitle("Transverse Energy [GeV]");
  fOutputList->Add(fHistClusterEt);
    
  fHistClusterN = new TH1F("fHistClusterN","Cluster N Distribution",30,0.0,30.0);
  fHistClusterN->GetYaxis()->SetTitle("Entries"); fHistClusterN->GetXaxis()->SetTitle("N");
  fOutputList->Add(fHistClusterN);
  
  fHistClusterM02 = new TH1F("fHistClusterM02","Cluster M02 Distribution",100,0.0,10.0);
  fHistClusterM02->GetYaxis()->SetTitle("Entries"); fHistClusterM02->GetXaxis()->SetTitle("M02");
  fOutputList->Add(fHistClusterM02);

  fHistClusterEN = new TH2F("fHistClusterEN","N vs Cluster Energy",100,0.0,50.0,30,0.0,30.0);
  fHistClusterEN->GetYaxis()->SetTitle("N"); fHistClusterEN->GetXaxis()->SetTitle("Energy [GeV]");
  fOutputList->Add(fHistClusterEN);

  fHistClusterEM02 = new TH2F("fHistClusterEM02","Cluster Energy vs M02",100,0.0,50.0,100,0.0,10.0);
  fHistClusterEM02->GetYaxis()->SetTitle("M02"); fHistClusterEM02->GetXaxis()->SetTitle("Energy [GeV]");
  fOutputList->Add(fHistClusterEM02);

  fHistClusterPhiEta = new TH2F("fHistClusterPhiEta","Cluster Pseudorapidity vs #phi",100,-1.0,1.0,100,0.0,6.29);
  fHistClusterPhiEta->GetXaxis()->SetTitle("#eta"); fHistClusterPhiEta->GetYaxis()->SetTitle("#phi [rad]");
  fOutputList->Add(fHistClusterPhiEta);

  fHistClusterEtN = new TH2F("fHistClusterEtN","N vs Cluster Transverse Energy",100,0.0,50.0,30,0.0,30.0);
  fHistClusterEtN->GetYaxis()->SetTitle("N"); fHistClusterEtN->GetXaxis()->SetTitle("Transverse Energy [GeV]");
  fOutputList->Add(fHistClusterEtN);

  fHistClusterEtM02 = new TH2F("fHistClusterEtM02","Cluster Transverse Energy vs M02",100,0.0,50.0,100,0.0,10.0);
  fHistClusterEtM02->GetYaxis()->SetTitle("M02"); fHistClusterEtM02->GetXaxis()->SetTitle("Transverse Energy [GeV]");
  fOutputList->Add(fHistClusterEtM02);

  fHistClusterdphiV0 = new TH1D("fHistClusterdphiV0","Cluster dphiV0 Distribution",100,0.0,TMath::Pi());
  fHistClusterdphiV0->GetYaxis()->SetTitle("Entries"); fHistClusterdphiV0->GetXaxis()->SetTitle("dphiV0 [rad]");
  fOutputList->Add(fHistClusterdphiV0);

  fHistTrackPt = new TH1F("fHistTrackPt","Track Transverse Momentum Distribution",100,0.0,30.0);
  fHistTrackPt->GetYaxis()->SetTitle("Entries"); fHistTrackPt->GetXaxis()->SetTitle("P_{t} [GeV/c]");
  fOutputList->Add(fHistTrackPt);

  fHistTrackEta = new TH1F("fHistTrackEta","Track Pseudorapidity Distribution",100,-1.0,1.0);
  fHistTrackEta->GetYaxis()->SetTitle("Entries"); fHistTrackEta->GetXaxis()->SetTitle("#eta");
  fOutputList->Add(fHistTrackEta);

  fHistTrackPhi = new TH1F("fHistTrackPhi","Track #phi Distribution",100,0.0,6.29);
  fHistTrackPhi->GetYaxis()->SetTitle("Entries"); fHistTrackPhi->GetXaxis()->SetTitle("#phi [rad]");
  fOutputList->Add(fHistTrackPhi);

  fHistTrackPhiEta = new TH2F("fHistTrackPhiEta","Track Pseudorapidity vs #phi",100,-1.0,1.0,100,0.0,6.29);
  fHistTrackPhiEta->GetXaxis()->SetTitle("#eta"); fHistTrackPhiEta->GetYaxis()->SetTitle("#phi [rad]");
  fOutputList->Add(fHistTrackPhiEta);

                    //  Et   M02  cent DeltaPhi    Cos[2*DeltaPhi]
  Int_t    bins[5] = {  500,  500,  100,     100,          100}; // binning
  Double_t min[5]  = {  0.0,  0.0,    0,     0.0,         -1.0}; // min x
  Double_t max[5]  = { 50.0, 10.0,  100,  TMath::Pi(),     1.0}; // max x
   
  fClusterV0 = new THnSparseF("fClusterV0","",5,bins,min,max);
  fClusterV0->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); fClusterV0->GetAxis(1)->SetTitle("M02"); fClusterV0->GetAxis(2)->SetTitle(" Centrality"); 
  fClusterV0->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); fClusterV0->GetAxis(4)->SetTitle("Cos[2*Delta(#phi)]");  
  fOutputList->Add(fClusterV0);
   
  fClusterV0A = new THnSparseF("fClusterV0A","",5,bins,min,max);
  fClusterV0A->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); fClusterV0A->GetAxis(1)->SetTitle("M02"); fClusterV0A->GetAxis(2)->SetTitle(" Centrality"); 
  fClusterV0A->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); fClusterV0A->GetAxis(4)->SetTitle("Cos[2*Delta(#phi)]");  
  fOutputList->Add(fClusterV0A);
   
  fClusterV0C = new THnSparseF("fClusterV0C","",5,bins,min,max);
  fClusterV0C->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); fClusterV0C->GetAxis(1)->SetTitle("M02"); fClusterV0C->GetAxis(2)->SetTitle(" Centrality"); 
  fClusterV0C->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); fClusterV0C->GetAxis(4)->SetTitle("Cos[2*Delta(#phi)]");  
  fOutputList->Add(fClusterV0C);
   
  fClusterTPC = new THnSparseF("fClusterTPC","",5,bins,min,max);
  fClusterTPC->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); fClusterTPC->GetAxis(1)->SetTitle("M02"); fClusterTPC->GetAxis(2)->SetTitle(" Centrality"); 
  fClusterTPC->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); fClusterTPC->GetAxis(4)->SetTitle("Cos[2*Delta(#phi)]");  
  fOutputList->Add(fClusterTPC);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::UserExec(Option_t *) 
{ // Main loop, called for each event
  fESDClusters = 0;
  fESDCells = 0;
  fAODClusters = 0;
  fAODCells = 0;

  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event\n"); return; }
  fHistEvtStat->Fill(0);

  fESD = dynamic_cast<AliESDEvent*>(event);
  if (!fESD) {
    fAOD = dynamic_cast<AliAODEvent*>(event);
    if(!fAOD){
      printf("ERROR: Could not retrieve the event\n");
      return;
    }
  }
  if (fESD) fHistEvtStat->Fill(1);
  if (fAOD) fHistEvtStat->Fill(2);
  if (!fESD && !fAOD) fHistEvtStat->Fill(3);
  
  fRunNumber = event->GetRunNumber();
  fInternalRunNum = GetInternalRunNum(fRunNumber);
  fHistRunNumStat->Fill(fInternalRunNum);
  
  // Centrality
  if (event->GetCentrality()) {
    fCentrality = event->GetCentrality()->GetCentralityPercentile("V0M"); // CL1 ?
    fHistCentrality->Fill(fCentrality);
    if (!IsCentralityDistrCorrected()) return;
    fHistCentralityCorrected->Fill(fCentrality);
    fHistEvtStat->Fill(4);
  } else {
    printf("ERROR: Could not retrieve the Centrality\n");
    return;
  }
  
  // Event Plane
  TObjArray *maps = (TObjArray*)fFlatContainer->GetObject(fRunNumber,"phosFlat");
  fTPCFlat = (AliEPFlattener*)maps->At(0);
  fV0AFlat = (AliEPFlattener*)maps->At(1);
  fV0CFlat = (AliEPFlattener*)maps->At(2);
  Bool_t flattenEP = kTRUE;
  fEventPlane = event->GetEventplane();
  if (fEventPlane) {
    VZEROEventPlane(flattenEP);
    fHistEvtStat->Fill(5);
  } else {
    printf("ERROR: Could not retrieve the Event Plane\n");
    return;
  }

  // Geom
  TObjArray *matEMCAL = (TObjArray*)fOADBContainer->GetObject(fRunNumber,"EmcalMatrices");
  for (Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++) {
    if (fGeoName=="EMCAL_FIRSTYEARV1" && mod>3) break;
    /*if(fESD)
      fGeom->SetMisalMatrix(fESD->GetEMCALMatrix(mod), mod);
      else*/
    // if(fVEvent->GetEMCALMatrix(mod))
    fGeomMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod);
    fGeom->SetMisalMatrix(fGeomMatrix[mod], mod);
  }

  // Cluster Readout
  if (fESD) {
    // TList *l = fESD->GetList();
    // if (1) {//fDebug){
    //   for(int nk=0;nk<l->GetEntries();nk++){
    //     TObject *obj = (TObject*)l->At(nk);
    //     TString oname = obj->GetName();
    // 	  //if(oname.Contains("lus"))
    //     printf("Object %d has a clus array named %s +++++++++\n",nk,oname.Data());
    //   }
    // }
    // fESDClusters = dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
    // // fESDCells = fESD->GetEMCALCells();
    // if(fDebug)
    //   printf("ESD cluster mult= %d\n",fESDClusters->GetEntriesFast());
    // if(fESDClusters->GetEntriesFast()<1){
    //   printf("ERROR: There are no EMCAL clusters in this event\n");
    //   return;
    // }
  } else if (fAOD) {
    if (event->FindListObject(fEMCALClustersListName))
      fAODClusters = dynamic_cast<TClonesArray*>(event->FindListObject(fEMCALClustersListName));
    if (!fAODClusters) {
      AliWarning(Form("Wrong name of list with clusters? <%s>", fEMCALClustersListName.Data()));
      return;
    }

    // fAODClusters = dynamic_cast<TClonesArray*>(fAOD->GetCaloClusters());
    // TObjArray* fAODClusters1 = dynamic_cast<TClonesArray*>(event->FindListObject("caloClusters"));
    // fAODCells = fAOD->GetEMCALCells();

    if (fDebug)
      printf("AOD cluster mult = %d\n", fAODClusters->GetEntriesFast());
    if (!fAODClusters || fAODClusters->GetEntriesFast()<1) {
      printf("ERROR: There are no EMCAL clusters in this event\n");
      return;
    }
    fHistEvtStat->Fill(6);
  }

  FillHistsCluster();  
  
  FillHistsTrack();

  PostData(1, fOutputList);
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALPi0V2ShSh::GetInternalRunNum(Int_t runnumber)
{ // convert runnumber into internal runnumber
  Int_t stdRunNum[45] = {167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 
                         168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168984, 
                         168988, 168992, 169035, 169091, 169094, 169138, 169143, 169144, 169145, 169148, 
                         169156, 169160, 169167, 169238, 169411, 169415, 169417, 169835, 169837, 169846, 
                         169855, 169858, 169859, 169923, 169965};

  for (int i = 0; i < 45; ++i) {
    if (runnumber==stdRunNum[i]) return i;
  }

  return -1;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0V2ShSh::IsCentralityDistrCorrected()
{
  Bool_t bIsNot = kFALSE;
  if (fCentrality<=10) {   // 0-10%
    TRandom3 *rndm = new TRandom3(0);
    Double_t Nrndm = rndm->Uniform(0.,1.);
    if (fCentrality<=1) {
      if(Nrndm > 0.77308) bIsNot = kTRUE;
    } else if (1<fCentrality && fCentrality<=2) {
      if (Nrndm>0.75863) bIsNot = kTRUE;
    } else if (2<fCentrality && fCentrality<=3){
      if (Nrndm>0.76365) bIsNot = kTRUE;
    } else if (3<fCentrality && fCentrality<=4){
      if (Nrndm>0.76763) bIsNot = kTRUE;
    } else if (4<fCentrality && fCentrality<=5){
      if (Nrndm>0.76251) bIsNot = kTRUE;
    } else if (5<fCentrality && fCentrality<=6){
      if (Nrndm>0.79069) bIsNot = kTRUE;
    } else if (6<fCentrality && fCentrality<=7){
      if (Nrndm>0.77669) bIsNot = kTRUE;
    } else if (7<fCentrality && fCentrality<=8){
      if (Nrndm>0.78537) bIsNot = kTRUE;
    } else if (8<fCentrality && fCentrality<=9){
      if (Nrndm>0.82727) bIsNot = kTRUE;
    } else if (9<fCentrality && fCentrality<=10){
      if (Nrndm>1) bIsNot = kTRUE;
    }
    delete rndm; rndm = NULL;
    return bIsNot;
  } else if (10<fCentrality && fCentrality<=50) {  // 10-50%
    TString centfired;
    if (fESD) {
      centfired = fESD->GetFiredTriggerClasses();
    } else {
      centfired = fAOD->GetFiredTriggerClasses();
    }
    if (!centfired.Contains("CVLN_B2-B-NOPF-ALLNOTRD") && 
        !centfired.Contains("CVLN_R1-B-NOPF-ALLNOTRD") && 
        !centfired.Contains("CSEMI_R1-B-NOPF-ALLNOTRD")) return kFALSE;
    else return kTRUE;
  } else return kTRUE; // other%

  return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::VZEROEventPlane(Bool_t flattenEP)
{ // Calculate the V0 Event Plane
  if (fEventPlane->GetQVector()) { 
    // fEPTPC = TVector2::Phi_0_2pi(fEventPlane->GetQVector()->Phi())/2.0; //if(fEPTPC>TMath::Pi()) {fEPTPC-=TMath::Pi();} 
    fEPTPC = fEventPlane->GetEventplane("Q");
  } else fEPTPC = -999;

  if (fEventPlane->GetQsub1() && fEventPlane->GetQsub2()) {
    // fEPTPCResolution = TMath::Cos(2.0*(fEventPlane->GetQsub1()->Phi()/2.0-fEventPlane->GetQsub2()->Phi()/2.0));
    fEPTPCResolution = TMath::Cos(2.*(fEventPlane->GetQsubRes()));
  } else fEPTPCResolution = -1;

  fEPV0  = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0",  fAOD)); if (fEPV0>TMath::Pi()) fEPV0-=TMath::Pi();
  fEPV0A = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0A", fAOD)); if (fEPV0A>TMath::Pi()) fEPV0A-=TMath::Pi();
  fEPV0C = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0C", fAOD)); if (fEPV0C>TMath::Pi()) fEPV0C-=TMath::Pi();

  if (fDebug)
    cout << "fEPV0/fEPV0A/fEPV0C = " << fEPV0 << "/" << fEPV0A << "/" << fEPV0C << endl;

  Double_t qx=0, qy=0, qxr=0, qyr=0;
  fEPV0AR = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 4, 5, 2, qxr, qyr)); if (fEPV0AR>TMath::Pi()) fEPV0AR-=TMath::Pi();
  fEPV0CR = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 2, 3, 2, qx,  qy)); if (fEPV0CR>TMath::Pi()) fEPV0CR-=TMath::Pi();
  qxr += qx; qyr += qy;
  fEPV0R   = TVector2::Phi_0_2pi(TMath::ATan2(qyr,qxr))/2.0; // equals to ring 2-5
  fEPV0AR4 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 4, 2, qx, qy)); if (fEPV0AR4>TMath::Pi()) fEPV0AR4-=TMath::Pi();
  fEPV0AR5 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 5, 2, qx, qy)); if (fEPV0AR5>TMath::Pi()) fEPV0AR5-=TMath::Pi();
  fEPV0AR6 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 6, 2, qx, qy)); if (fEPV0AR6>TMath::Pi()) fEPV0AR6-=TMath::Pi();
  fEPV0AR7 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 7, 2, qx, qy)); if (fEPV0AR7>TMath::Pi()) fEPV0AR7-=TMath::Pi();
  fEPV0CR0 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 0, 2, qx, qy)); if (fEPV0CR0>TMath::Pi()) fEPV0CR0-=TMath::Pi();
  fEPV0CR1 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 1, 2, qx, qy)); if (fEPV0CR1>TMath::Pi()) fEPV0CR1-=TMath::Pi();
  fEPV0CR2 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 2, 2, qx, qy)); if (fEPV0CR2>TMath::Pi()) fEPV0CR2-=TMath::Pi();
  fEPV0CR3 = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fAOD, 3, 2, qx, qy)); if (fEPV0CR3>TMath::Pi()) fEPV0CR3-=TMath::Pi();

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

  fHistEPV0MB->Fill(fEPV0);
  fHistEPV0RMB->Fill(fEPV0R);
  fHistEPV0AMB->Fill(fEPV0A);
  fHistEPV0CMB->Fill(fEPV0C);  
  fHistEPTPCMB->Fill(fEPTPC);

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

  fHistDiffV0A_V0CR0->Fill(fCentrality, TMath::Cos(2.0*(fEPV0A - fEPV0CR0)));
  fHistDiffV0A_V0CR3->Fill(fCentrality, TMath::Cos(2.0*(fEPV0A - fEPV0CR3)));
  fHistDiffV0CR0_V0CR3->Fill(fCentrality, TMath::Cos(2.0*(fEPV0CR0 - fEPV0CR3)));
  fHistDiffV0C_V0AR4->Fill(fCentrality, TMath::Cos(2.0*(fEPV0C - fEPV0AR4)));
  fHistDiffV0C_V0AR7->Fill(fCentrality, TMath::Cos(2.0*(fEPV0C - fEPV0AR7)));
  fHistDiffV0AR4_V0AR7->Fill(fCentrality, TMath::Cos(2.0*(fEPV0AR4 - fEPV0AR7)));   
  fHistDiffV0AR_V0CR->Fill(fCentrality, TMath::Cos(2.0*(fEPV0AR - fEPV0CR)));

  // run-by-run QA
  fHistRBRCosV0A->Fill(fInternalRunNum, TMath::Cos(2*fEPV0A));
  fHistRBRSinV0A->Fill(fInternalRunNum, TMath::Sin(2*fEPV0A));
  fHistRBRCosV0C->Fill(fInternalRunNum, TMath::Cos(2*fEPV0C));
  fHistRBRSinV0C->Fill(fInternalRunNum, TMath::Sin(2*fEPV0C));
  fHistRBRCosTPC->Fill(fInternalRunNum, TMath::Cos(2*fEPTPC));
  fHistRBRSinTPC->Fill(fInternalRunNum, TMath::Sin(2*fEPTPC));
}

//________________________________________________________________________
Double_t  AliAnalysisTaskEMCALPi0V2ShSh::ApplyFlatteningV0A(Double_t phi, Double_t c)
{
  if (fV0AFlat) return fV0AFlat->MakeFlat(phi,c);
  return phi;
}

Double_t  AliAnalysisTaskEMCALPi0V2ShSh::ApplyFlatteningV0C(Double_t phi, Double_t c)
{ 
  if (fV0CFlat) return fV0CFlat->MakeFlat(phi,c);
  return phi;
}

Double_t  AliAnalysisTaskEMCALPi0V2ShSh::ApplyFlatteningTPC(Double_t phi, Double_t c)
{
  if (fTPCFlat) return fTPCFlat->MakeFlat(phi,c);
  return phi;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::FillHistsCluster()
{ // Fill cluster histograms
  TObjArray *clustersArray = fESDClusters;
  if (!clustersArray) clustersArray = fAODClusters;
  if (!clustersArray) return;

  const Int_t nclus = clustersArray->GetEntries();
  if (nclus==0) return;

  Float_t pos[3]; 
  for (Int_t iclus = 0; iclus < nclus; iclus++) {
    fHistClusterStat->Fill(0);
    AliVCluster *clus = (AliVCluster*)clustersArray->At(iclus);
    if (!clus) continue;
    if (!clus->IsEMCAL()) continue;
    fHistClusterStat->Fill(1);

    if (!IsGoodClusterV1(clus)) continue;
    if (!IsGoodPi0(clus)) continue;
    fHistClusterStat->Fill(2);
    
    clus->GetPosition(pos);
    TVector3 vpos(pos[0],pos[1],pos[2]);
    Double_t N = clus->GetNCells();
    Double_t E = clus->E();
    Double_t Et = (E / (TMath::CosH(vpos.Eta())));
    Double_t eta = vpos.Eta();
    Double_t phi = vpos.Phi();
    Double_t M02 = clus->GetM02();

    Double_t dphiV0 = TVector2::Phi_0_2pi(phi-fEPV0R); if(dphiV0>TMath::Pi()) {dphiV0-=TMath::Pi();} // ??
    Double_t dphiV0A = TVector2::Phi_0_2pi(phi-fEPV0A); if(dphiV0A>TMath::Pi()) {dphiV0A-=TMath::Pi();}
    Double_t dphiV0C = TVector2::Phi_0_2pi(phi-fEPV0C); if(dphiV0C>TMath::Pi()) {dphiV0C-=TMath::Pi();}
    Double_t dphiTPC = TVector2::Phi_0_2pi(phi-fEPTPC); if(dphiTPC>TMath::Pi()) {dphiTPC-=TMath::Pi();}

    Double_t DataV0[5];
    DataV0[0] = Et; 	  
    DataV0[1] = M02; 
    DataV0[2] = fCentrality;
    DataV0[3] = dphiV0;
    DataV0[4] = TMath::Cos(2.0*(dphiV0));
    fClusterV0->Fill(DataV0);

    Double_t DataV0A[5];
    DataV0A[0] = Et; 	  
    DataV0A[1] = M02; 
    DataV0A[2] = fCentrality;
    DataV0A[3] = dphiV0A;
    DataV0A[4] = TMath::Cos(2.0*(dphiV0A));
    fClusterV0A->Fill(DataV0A);

    Double_t DataV0C[5];
    DataV0C[0] = Et; 	  
    DataV0C[1] = M02; 
    DataV0C[2] = fCentrality;
    DataV0C[3] = dphiV0C;
    DataV0C[4] = TMath::Cos(2.0*(dphiV0C));
    fClusterV0C->Fill(DataV0C);

    Double_t DataTPC[5];
    DataTPC[0] = Et; 	  
    DataTPC[1] = M02; 
    DataTPC[2] = fCentrality;
    DataTPC[3] = dphiTPC;
    DataTPC[4] = TMath::Cos(2.0*(dphiTPC));
    fClusterTPC->Fill(DataTPC); 

    fHistClusterE->Fill(E);
    fHistClusterEt->Fill(Et);
    fHistClusterM02->Fill(M02);
    fHistClusterN->Fill(N);
    fHistClusterPhi->Fill(phi);
    fHistClusterEta->Fill(eta);
    fHistClusterPhiEta->Fill(eta,phi);
    fHistClusterEN->Fill(E,N);
    fHistClusterEM02->Fill(E,M02);
    fHistClusterEtN->Fill(Et,N);
    fHistClusterEtM02->Fill(Et,M02);
    fHistClusterdphiV0->Fill(dphiV0);
  }

  if (fDebug)
    cout << "Good Pi0 cluster mult = " << fHistClusterStat->GetBinContent(3) << endl;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0V2ShSh::IsGoodClusterV1(const AliVCluster *c) const
{
  if (!c) return kFALSE;

  if (c->GetNCells()<2.) return kFALSE;
  fHistClusterStat->Fill(3);

  if (c->E()<1.) return kFALSE;
  fHistClusterStat->Fill(4);

  Short_t id = -1;
  Double_t maxE = GetMaxCellEnergy(c, id);
  if ((1. - double(GetCrossEnergy(c,id))/maxE) > 0.97) return kFALSE;
  fHistClusterStat->Fill(5);

  Float_t pos1[3];
  c->GetPosition(pos1);
  TVector3 clsPos(pos1);
  Double_t eta = clsPos.Eta();

  if (TMath::Abs(eta) > 0.65) return kFALSE;
  fHistClusterStat->Fill(6);

  if (!IsWithinFiducialVolume(id)) return kFALSE;
  fHistClusterStat->Fill(7);

  if (c->GetM02() <0.3) return kFALSE;
  fHistClusterStat->Fill(8);

  Double_t dr = TMath::Sqrt(c->GetTrackDx() * c->GetTrackDx() + 
                            c->GetTrackDz() * c->GetTrackDz());
  if (dr<0.025) return kFALSE;
  fHistClusterStat->Fill(9);

  return kTRUE;
}

Double_t AliAnalysisTaskEMCALPi0V2ShSh::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{ // Get maximum energy of attached cell
  id = -1;

  AliVCaloCells *cells = 0;
  if (fESD) {
    cells = fESD->GetEMCALCells();
  } else {
    cells = fAOD->GetEMCALCells();
  }
  if (!cells) return 0;

  Double_t maxe = 0;
  const Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Double_t e = cells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}

Bool_t AliAnalysisTaskEMCALPi0V2ShSh::IsWithinFiducialVolume(Short_t id) const
{ // Check if cell is within given fiducial volume
  Double_t fNFiducial = 1;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;

  Bool_t okrow = kFALSE;
  Bool_t okcol = kFALSE;

  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom) return kFALSE;

  Int_t cellAbsId = id;
  geom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
  geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);

  // Check rows/phi
  if (iSupMod < 10) {
    if (iphi >= fNFiducial && iphi < 24-fNFiducial) okrow = kTRUE;
  } else {
    if (iphi >= fNFiducial && iphi < 12-fNFiducial) okrow = kTRUE;
  }

  // Check columns/eta
  Bool_t noEMCALBorderAtEta0 = kTRUE;
  if (!noEMCALBorderAtEta0) {
    if (ieta > fNFiducial && ieta < 48-fNFiducial) okcol = kTRUE;
  } else {
    if (iSupMod%2==0) {
      if (ieta >= fNFiducial) okcol = kTRUE;
    } else {
      if (ieta < 48-fNFiducial) okcol = kTRUE;
    }
  }
  if (okrow && okcol) return kTRUE;

  return kFALSE;
}

Double_t AliAnalysisTaskEMCALPi0V2ShSh::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax) const
{ // Calculate the energy of cross cells around the leading cell
  AliVCaloCells *cells;
  if (fESD) {
    cells = fESD->GetEMCALCells();
  } else {
    cells = fAOD->GetEMCALCells();
  }
  if (!cells) return 0;

  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom) return 0;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  Double_t crossEnergy = 0.;

  geom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = cluster->GetCellAbsId(i);
    geom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t aphidiff = TMath::Abs(iphi-iphis);
    if (aphidiff>1) continue;
    Int_t aetadiff = TMath::Abs(ieta-ietas);
    if (aetadiff>1) continue;
    if ((aphidiff==1 && aetadiff==0) ||
        (aphidiff==0 && aetadiff==1)) {
      crossEnergy += cells->GetCellAmplitude(cellAbsId);
    }
  }

  return crossEnergy;
}

//_____________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0V2ShSh::IsGoodPi0(const AliVCluster *c)
{
  Double_t E = c->E();
  // Double_t Et = (E / (TMath::CosH(vpos.Eta())));
  Double_t M02 = c->GetM02();

  Double_t M02Min = exp(2.135-0.245*E);
  Double_t M02Max = exp(0.0662-0.0201*E) - 0.0955 + 0.00186*E + 9.91/E; // NLM = 1
  // Double_t M02Max = exp(0.353-0.0264*E) - 0.524 + 0.00559*E + 21.9/E; // NLM = 2

  if (M02>M02Max || M02<M02Min) return kFALSE;
  if (fDebug)
    cout << "E/M02 = " << E << "/" << M02 << " (" << M02Min << ", " << M02Max << ")" << endl;

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::FillHistsTrack()
{ // Fill track histograms
  if (fESD) {
    for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* track = (AliESDtrack*)fESD->GetTrack(iTracks);
      if (!track){
        printf("ERROR: Could not retreive esdtrack\n");
        continue;
      }
      fHistTrackPt->Fill(track->Pt());
      fHistTrackEta->Fill(track->Eta());
      fHistTrackPhi->Fill(track->Phi());
      fHistTrackPhiEta->Fill(track->Eta(),track->Phi());
    }
  } else if (fAOD) {
    for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
      AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTracks);
      if (!track){
        printf("ERROR: Could not retreive esdtrack\n");
        continue;
      }
      fHistTrackPt->Fill(track->Pt());
      fHistTrackEta->Fill(track->Eta());
      fHistTrackPhi->Fill(track->Phi());
      fHistTrackPhiEta->Fill(track->Eta(),track->Phi());
    }    
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::Terminate(Option_t *) 
{
}
