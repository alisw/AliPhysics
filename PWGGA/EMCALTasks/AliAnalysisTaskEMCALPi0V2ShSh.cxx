// $Id: AliAnalysisTaskEMCALPi0V2ShSh.cxx$

#include "AliAnalysisTaskEMCALPi0V2ShSh.h"

//Root include files 
//#include <Riostream.h>
#include <TParticle.h>
#include <TRefArray.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TMath.h>
#include <TVirtualFFT.h>

//AliRoot include files 
#include "AliAnalysisTaskSE.h"
#include "AliRunLoader.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliStack.h"
#include "AliEMCALGeometry.h"
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
  AliAnalysisTaskSE(), 
  fEventPlane(0),
  fCentralityV0M(99.),
  fESDClusters(0),
  fAODClusters(0),
  fESDCells(0),
  fAODCells(0),
  fGeom(0),
  fGeoName("EMCAL_COMPLETEV1"),
  fOADBContainer(0),
  fESD(0),
  fAOD(0),
  fOutputList(0),
  fEPTPC(-999.),
  fEPTPCResolution(0.),
  fEPV0(-999.),
  fEPV0A(-999.),
  fEPV0C(-999.),
  fEPV0Ar(-999.),
  fEPV0Cr(-999.),
  fEPV0r(-999.),
  fEPV0A4r(-999.),
  fEPV0A5r(-999.),
  fEPV0A6r(-999.),
  fEPV0A7r(-999.),
  fEPV0C0r(-999.),
  fEPV0C1r(-999.),
  fEPV0C2r(-999.),
  fEPV0C3r(-999.),
  fHistAllcentV0(0),
  fHistAllcentV0r(0),
  fHistAllcentV0A(0),
  fHistAllcentV0C(0),
  fHistAllcentTPC(0),
  fHistEPTPC(0),
  fHistEPTPCResolution(0),
  fHistEPV0(0),
  fHistEPV0A(0),
  fHistEPV0C(0),
  fHistEPV0Ar(0),
  fHistEPV0Cr(0),
  fHistEPV0r(0),
  fHistEPV0A4r(0),
  fHistEPV0A7r(0),
  fHistEPV0C0r(0),
  fHistEPV0C3r(0),
  fHistdifV0A_V0C0r(0),
  fHistdifV0A_V0C3r(0),
  fHistdifV0C0r_V0C3r(0),
  fHistdifV0C_V0A4r(0),
  fHistdifV0C_V0A7r(0),
  fHistdifV0A4r_V0A7r(0),
  fHistdifV0Ar_V0Cr(0),
  fHistClusterEta(0),
  fHistClusterPhi(0),
  fHistClusterE(0),
  fHistClusterEt(0),
  fHistClusterN(0),
  fHistClusterM02(0),
  fHistClusterEN(0),
  fHistClusterEM02(0),
  fHistClusterPhiEta(0),
  fHistClusterEtN(0),
  fHistClusterEtM02(0),
  fHistClusterdphiV0(0),
  fHistTrackPt(0),
  fHistTrackEta(0),
  fHistTrackPhi(0),
  fHistTrackPhiEta(0),
  fClusterPbV0(0),
  fClusterPbV0A(0),
  fClusterPbV0C(0),
  fClusterPbTPC(0)
{
  // Default constructor.
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0V2ShSh::AliAnalysisTaskEMCALPi0V2ShSh(const char *name) : 
  AliAnalysisTaskSE(name), 
  fEventPlane(0), 
  fCentralityV0M(99.), 
  fESDClusters(0),
  fAODClusters(0),
  fESDCells(0),
  fAODCells(0),
  fGeom(0),
  fGeoName("EMCAL_COMPLETEV1"),
  fOADBContainer(0),
  fESD(0), 
  fAOD(0),
  fOutputList(0),
  fEPTPC(-999.), 
  fEPTPCResolution(0.), 
  fEPV0(-999.), 
  fEPV0A(-999.), 
  fEPV0C(-999.),
  fEPV0Ar(-999.),
  fEPV0Cr(-999.),
  fEPV0r(-999.),
  fEPV0A4r(-999.),
  fEPV0A5r(-999.),
  fEPV0A6r(-999.), 
  fEPV0A7r(-999.),
  fEPV0C0r(-999.), 
  fEPV0C1r(-999.), 
  fEPV0C2r(-999.), 
  fEPV0C3r(-999.),
  fHistAllcentV0(0), 
  fHistAllcentV0r(0), 
  fHistAllcentV0A(0), 
  fHistAllcentV0C(0), 
  fHistAllcentTPC(0),
  fHistEPTPC(0), 
  fHistEPTPCResolution(0),
  fHistEPV0(0), 
  fHistEPV0A(0), 
  fHistEPV0C(0), 
  fHistEPV0Ar(0), 
  fHistEPV0Cr(0), 
  fHistEPV0r(0), 
  fHistEPV0A4r(0),
  fHistEPV0A7r(0),
  fHistEPV0C0r(0),
  fHistEPV0C3r(0),
  fHistdifV0A_V0C0r(0),
  fHistdifV0A_V0C3r(0),
  fHistdifV0C0r_V0C3r(0),
  fHistdifV0C_V0A4r(0), 
  fHistdifV0C_V0A7r(0), 
  fHistdifV0A4r_V0A7r(0), 
  fHistdifV0Ar_V0Cr(0),
  fHistClusterEta(0), 
  fHistClusterPhi(0), 
  fHistClusterE(0), 
  fHistClusterEt(0), 
  fHistClusterN(0), 
  fHistClusterM02(0),
  fHistClusterEN(0), 
  fHistClusterEM02(0), 
  fHistClusterPhiEta(0), 
  fHistClusterEtN(0), 
  fHistClusterEtM02(0), 
  fHistClusterdphiV0(0),
  fHistTrackPt(0), 
  fHistTrackEta(0), 
  fHistTrackPhi(0), 
  fHistTrackPhiEta(0),
  fClusterPbV0(0), 
  fClusterPbV0A(0), 
  fClusterPbV0C(0), 
  fClusterPbTPC(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::UserCreateOutputObjects()
{
  // Create histograms, called once.
    
  fESDClusters = new TObjArray();
  fAODClusters = new TObjArray();
  fGeom = AliEMCALGeometry::GetInstance(fGeoName.Data());
  fOADBContainer = new AliOADBContainer("AliEMCALgeo");
  fOADBContainer->InitFromFile(Form("$ALICE_PHYSICS/OADB/EMCAL/EMCALlocal2master.root"),"AliEMCALgeo");

  fOutputList = new TList();
  fOutputList->SetOwner();// Container cleans up all histos (avoids leaks in merging)   
  
  fHistTrackPt = new TH1F("fHistTrackPt","Track Transverse Momentum Distribution of Pb+Pb",100,0.0,30.0);
  fHistTrackPt->GetYaxis()->SetTitle("Entries"); fHistTrackPt->GetXaxis()->SetTitle("P_{t} [GeV/c]");
  fOutputList->Add(fHistTrackPt);

  fHistTrackEta = new TH1F("fHistTrackEta","Track Pseudorapidity Distribution of Pb+Pb",100,-1.0,1.0);
  fHistTrackEta->GetYaxis()->SetTitle("Entries"); fHistTrackEta->GetXaxis()->SetTitle("#eta");
  fOutputList->Add(fHistTrackEta);

  fHistTrackPhi = new TH1F("fHistTrackPhi","Track #phi Distribution of Pb+Pb",100,0.0,6.29);
  fHistTrackPhi->GetYaxis()->SetTitle("Entries"); fHistTrackPhi->GetXaxis()->SetTitle("#phi [rad]");
  fOutputList->Add(fHistTrackPhi);

  fHistTrackPhiEta = new TH2F("fHistTrackPhiEta","Track Pseudorapidity vs #phi of Pb+Pb",100,-1.0,1.0,100,0.0,6.29);
  fHistTrackPhiEta->GetXaxis()->SetTitle("#eta"); fHistTrackPhiEta->GetYaxis()->SetTitle("#phi [rad]");
  fOutputList->Add(fHistTrackPhiEta);
  
  fHistClusterEta = new TH1F("fHistClusterEta","Cluster Pseudorapidity Distribution of Pb+Pb",100,-1.0,1.0);
  fHistClusterEta->GetYaxis()->SetTitle("Entries"); fHistClusterEta->GetXaxis()->SetTitle("#eta");
  fOutputList->Add(fHistClusterEta);

  fHistClusterPhi = new TH1F("fHistClusterPhi","Cluster #phi Distribution of Pb+Pb",100,0.0,6.29);
  fHistClusterPhi->GetYaxis()->SetTitle("Entries"); fHistClusterPhi->GetXaxis()->SetTitle("#phi [rad]");
  fOutputList->Add(fHistClusterPhi);

  fHistClusterPhiEta = new TH2F("fHistClusterPhiEta","Cluster Pseudorapidity vs #phi of Pb+Pb",100,-1.0,1.0,100,0.0,6.29);
  fHistClusterPhiEta->GetXaxis()->SetTitle("#eta"); fHistClusterPhiEta->GetYaxis()->SetTitle("#phi [rad]");
  fOutputList->Add(fHistClusterPhiEta);

  fHistClusterM02 = new TH1F("fHistClusterM02","Cluster M02 Distribution of Pb+Pb",100,0.0,3.0);
  fHistClusterM02->GetYaxis()->SetTitle("Entries"); fHistClusterM02->GetXaxis()->SetTitle("M02");
  fOutputList->Add(fHistClusterM02);

  fHistClusterE = new TH1F("fHistClusterE","Cluster Energy Distribution of Pb+Pb",100,0.0,20.0);
  fHistClusterE->GetYaxis()->SetTitle("Entries"); fHistClusterE->GetXaxis()->SetTitle("Energy [GeV]");
  fOutputList->Add(fHistClusterE);
  
  fHistClusterEt = new TH1F("fHistClusterEt","Cluster Transverse Energy Distribution of Pb+Pb",100,0.0,20.0);
  fHistClusterEt->GetYaxis()->SetTitle("Entries"); fHistClusterEt->GetXaxis()->SetTitle("Transverse Energy [GeV]");
  fOutputList->Add(fHistClusterEt);
  
  fHistClusterEM02 = new TH2F("fHistClusterEM02","Cluster Energy vs M02 of Pb+Pb",100,0.0,20.0,100,0.0,3.0);
  fHistClusterEM02->GetYaxis()->SetTitle("M02"); fHistClusterEM02->GetXaxis()->SetTitle("Energy [GeV]");
  fOutputList->Add(fHistClusterEM02);

  fHistClusterEtM02 = new TH2F("fHistClusterEtM02","Cluster Transverse Energy vs M02 of Pb+Pb",100,0.0,20.0,100,0.0,3.0);
  fHistClusterEtM02->GetYaxis()->SetTitle("M02"); fHistClusterEtM02->GetXaxis()->SetTitle("Transverse Energy [GeV]");
  fOutputList->Add(fHistClusterEtM02);
  
  fHistClusterN = new TH1F("fHistClusterN","Cluster N Distribution of Pb+Pb",30,0.0,30.0);
  fHistClusterN->GetYaxis()->SetTitle("Entries"); fHistClusterN->GetXaxis()->SetTitle("N");
  fOutputList->Add(fHistClusterN);
  
  fHistClusterEN = new TH2F("fHistClusterEN","N vs Cluster Energy of Pb+Pb",100,0.0,20.0,30,0.0,30.0);
  fHistClusterEN->GetYaxis()->SetTitle("N"); fHistClusterEN->GetXaxis()->SetTitle("Energy [GeV]");
  fOutputList->Add(fHistClusterEN);

  fHistClusterEtN = new TH2F("fHistClusterEtN","N vs Cluster Transverse Energy of Pb+Pb",100,0.0,20.0,30,0.0,30.0);
  fHistClusterEtN->GetYaxis()->SetTitle("N"); fHistClusterEtN->GetXaxis()->SetTitle("Transverse Energy [GeV]");
  fOutputList->Add(fHistClusterEtN);

  fHistClusterdphiV0 = new TH1D("fHistClusterdphiV0","Cluster dphiV0 Distribution of Pb+Pb",100,0.0,TMath::Pi());
  fHistClusterdphiV0->GetYaxis()->SetTitle("Entries"); fHistClusterdphiV0->GetXaxis()->SetTitle("dphiV0 [rad]");
  fOutputList->Add(fHistClusterdphiV0);
  
  fHistEPV0 = new TH2F("fHistEPV0","V0 Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0->GetYaxis()->SetTitle("V0 Event Plane Angle [rad]"); fHistEPV0->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0);

  fHistEPV0A = new TH2F("fHistEPV0A","V0A Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0A->GetYaxis()->SetTitle("V0A Event Plane Angle [rad]"); fHistEPV0A->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0A);

  fHistEPV0C = new TH2F("fHistEPV0","V0C Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0C->GetYaxis()->SetTitle("V0C Event Plane Angle [rad]"); fHistEPV0C->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0C);

  fHistEPV0r = new TH2F("fHistEPV0r","V0r Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0r->GetYaxis()->SetTitle("V0r Event Plane Angle [rad]"); fHistEPV0r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0r);

  fHistEPV0Ar = new TH2F("fHistEPV0","V0Ar Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0Ar->GetYaxis()->SetTitle("V0Ar Event Plane Angle [rad]"); fHistEPV0Ar->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0Ar);

  fHistEPV0Cr = new TH2F("fHistEPV0Cr","V0Cr Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0Cr->GetYaxis()->SetTitle("V0Cr Event Plane Angle [rad]"); fHistEPV0Cr->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0Cr);

  fHistEPV0A4r = new TH2F("fHistEPV0Cr","V0A4r Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0A4r->GetYaxis()->SetTitle("V0A4r Event Plane Angle [rad]"); fHistEPV0A4r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0A4r);

  fHistEPV0A7r = new TH2F("fHistEPV0Cr","V0A7r Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0A7r->GetYaxis()->SetTitle("V0A7r Event Plane Angle [rad]"); fHistEPV0A7r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0A7r);

  fHistEPV0C0r = new TH2F("fHistEPV0Cr","V0C0r Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0C0r->GetYaxis()->SetTitle("V0C0r Event Plane Angle [rad]"); fHistEPV0C0r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0C0r);

  fHistEPV0C3r = new TH2F("fHistEPV0Cr","V0C3r Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPV0C3r->GetYaxis()->SetTitle("V0C3r Event Plane Angle [rad]"); fHistEPV0C3r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPV0C3r);

  fHistdifV0A_V0C0r = new TH2F("fHistdifV0A_V0C0r","(V0A - V0C0r) vs Centrality of Pb+Pb",100,0,100,100,-1.0,1.0);
  fHistdifV0A_V0C0r->GetYaxis()->SetTitle("Cos[2*(V0A - V0C0r)]"); fHistdifV0A_V0C0r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistdifV0A_V0C0r);

  fHistdifV0A_V0C3r = new TH2F("fHistdifV0A_V0C3r","(V0A - V0C3r) vs Centrality of Pb+Pb",100,0,100,100,-1.0,1.0);
  fHistdifV0A_V0C3r->GetYaxis()->SetTitle("Cos[2*(V0A - V0C3r)]"); fHistdifV0A_V0C3r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistdifV0A_V0C3r);

  fHistdifV0C0r_V0C3r = new TH2F("fHistdifV0C0r_V0C3r","(V0C0r - V0C3r) vs Centrality of Pb+Pb",100,0,100,100,-1.0,1.0);
  fHistdifV0C0r_V0C3r->GetYaxis()->SetTitle("Cos[2*(V0C0r - V0C3r)]"); fHistdifV0C0r_V0C3r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistdifV0C0r_V0C3r);

  fHistdifV0C_V0A4r = new TH2F("fHistdifV0C_V0A4r","(V0C - V0A4r) vs Centrality of Pb+Pb",100,0,100,100,-1.0,1.0);
  fHistdifV0C_V0A4r->GetYaxis()->SetTitle("Cos[2*(V0C - V0C4r)]"); fHistdifV0C_V0A4r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistdifV0C_V0A4r);

  fHistdifV0C_V0A7r = new TH2F("fHistdifV0C_V0A3r","(V0C - V0A7r) vs Centrality of Pb+Pb",100,0,100,100,-1.0,1.0);
  fHistdifV0C_V0A7r->GetYaxis()->SetTitle("Cos[2*(V0C - V0A7r)]"); fHistdifV0C_V0A7r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistdifV0C_V0A7r);

  fHistdifV0A4r_V0A7r = new TH2F("fHistdifV0A4r_V0A7r","(V0A4r - V0A7r) vs Centrality of Pb+Pb",100,0,100,100,-1.0,1.0);
  fHistdifV0A4r_V0A7r->GetYaxis()->SetTitle("Cos[2*(V0A4r - V0A7r)]"); fHistdifV0A4r_V0A7r->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistdifV0A4r_V0A7r);

  fHistdifV0Ar_V0Cr = new TH2F("fHistdifV0Ar_V0Cr","(V0Ar - V0Cr) vs Centrality of Pb+Pb",100,0,100,100,-1.0,1.0);
  fHistdifV0Ar_V0Cr->GetYaxis()->SetTitle("Cos[2*(V0Ar - V0Cr)]"); fHistdifV0Ar_V0Cr->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistdifV0Ar_V0Cr);

  fHistAllcentV0 = new TH1F("fHistAllcentV0","V0 Event Plane Angle of Pb+Pb",100,0.0,TMath::Pi());
  fHistAllcentV0->GetXaxis()->SetTitle("V0 Event Plane Angle [rad]"); fHistAllcentV0->GetYaxis()->SetTitle("");
  fOutputList->Add(fHistAllcentV0);

  fHistAllcentV0r = new TH1F("fHistAllcentV0r","V0r Event Plane Angle of Pb+Pb",100,0.0,TMath::Pi());
  fHistAllcentV0r->GetXaxis()->SetTitle("V0r Event Plane Angle [rad]"); fHistAllcentV0r->GetYaxis()->SetTitle("");
  fOutputList->Add(fHistAllcentV0r);

  fHistAllcentV0A = new TH1F("fHistAllcentV0A","V0A Event Plane Angle of Pb+Pb",100,0.0,TMath::Pi());
  fHistAllcentV0A->GetXaxis()->SetTitle("V0A Event Plane Angle [rad]"); fHistAllcentV0A->GetYaxis()->SetTitle("");
  fOutputList->Add(fHistAllcentV0A);

  fHistAllcentV0C = new TH1F("fHistAllcentV0C","V0C Event Plane Angle of Pb+Pb",100,0.0,TMath::Pi());
  fHistAllcentV0C->GetXaxis()->SetTitle("V0C Event Plane Angle [rad]"); fHistAllcentV0C->GetYaxis()->SetTitle("");
  fOutputList->Add(fHistAllcentV0C);

  fHistAllcentTPC = new TH1F("fHistAllcentTPC","TPC Event Plane Angle of Pb+Pb",100,0.0,TMath::Pi());
  fHistAllcentTPC->GetXaxis()->SetTitle("TPC Event Plane Angle [rad]"); fHistAllcentTPC->GetYaxis()->SetTitle("");
  fOutputList->Add(fHistAllcentTPC);

  fHistEPTPC = new TH2F("fHistEPTPC","TPC Event Plane Angle vs Centrality of Pb+Pb",100,0,100,100,0.0,TMath::Pi());
  fHistEPTPC->GetYaxis()->SetTitle("Event Plane Angle [rad]"); fHistEPTPC->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPTPC);

  fHistEPTPCResolution = new TH2F("fHistEPTPCResolution","TPC Resolution vs Centrality of Pb+Pb",100,0,100,100,0.0,1.0);
  fHistEPTPCResolution->GetYaxis()->SetTitle("TPC Resolution"); fHistEPTPCResolution->GetXaxis()->SetTitle("V0M Centrality");
  fOutputList->Add(fHistEPTPCResolution);
 

		    //  Et   M02  V0Mcent DeltaPhi    Cos[2*DeltaPhi]
  Int_t    bins[5] = {  500, 250,  100,     100,          100  }; // binning
  Double_t min[5]  = {  0.0, 0.0,    0,     0.0,         -1.0}; // min x
  Double_t max[5]  = { 50.0, 2.5,  100,  TMath::Pi(),     1.0}; // max x
	 
  fClusterPbV0 = new THnSparseF("fClusterPbV0","",5,bins,min,max);
  fClusterPbV0->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); fClusterPbV0->GetAxis(1)->SetTitle("M02"); fClusterPbV0->GetAxis(2)->SetTitle("V0M Centrality"); 
  fClusterPbV0->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); fClusterPbV0->GetAxis(4)->SetTitle("Cos[2*Delta(#phi)]");  
  fOutputList->Add(fClusterPbV0);
 	 
  fClusterPbV0A = new THnSparseF("fClusterPbV0A","",5,bins,min,max);
  fClusterPbV0A->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); fClusterPbV0A->GetAxis(1)->SetTitle("M02"); fClusterPbV0A->GetAxis(2)->SetTitle("V0M Centrality"); 
  fClusterPbV0A->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); fClusterPbV0A->GetAxis(4)->SetTitle("Cos[2*Delta(#phi)]");  
  fOutputList->Add(fClusterPbV0A);
	 
  fClusterPbV0C = new THnSparseF("fClusterPbV0C","",5,bins,min,max);
  fClusterPbV0C->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); fClusterPbV0C->GetAxis(1)->SetTitle("M02"); fClusterPbV0C->GetAxis(2)->SetTitle("V0M Centrality"); 
  fClusterPbV0C->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); fClusterPbV0C->GetAxis(4)->SetTitle("Cos[2*Delta(#phi)]");  
  fOutputList->Add(fClusterPbV0C);
	 
  fClusterPbTPC = new THnSparseF("fClusterPbTPC","",5,bins,min,max);
  fClusterPbTPC->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); fClusterPbTPC->GetAxis(1)->SetTitle("M02"); fClusterPbTPC->GetAxis(2)->SetTitle("V0M Centrality"); 
  fClusterPbTPC->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); fClusterPbTPC->GetAxis(4)->SetTitle("Cos[2*Delta(#phi)]");  
  fOutputList->Add(fClusterPbTPC);


  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  fESDClusters = 0;
  fESDCells = 0;
  fAODClusters = 0;
  fAODCells = 0;

  // Create pointer to reconstructed event
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event\n"); return; }

  fESD = dynamic_cast<AliESDEvent*>(event);
  if (!fESD) {
    fAOD =  dynamic_cast<AliAODEvent*>(event);
    if(!fAOD){
      printf("ERROR: Could not retrieve the event\n");
      return;
    }
  }

  //Bool_t isSelected =0;      
  //isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kSemiCentral));
  //if(!isSelected) { return; }
  
  if(event->GetCentrality()){
    fCentralityV0M = event->GetCentrality()->GetCentralityPercentile("V0M");
  } else {
    printf("ERROR: Could not retrieve the Centrality\n");
    return;
  }	
  
  fEventPlane = event->GetEventplane(); 
  if (fEventPlane) {
     VZEROEventPlane();
  } else {
     printf("ERROR: Could not retrieve the Centrality\n");
     return;
  }
  Int_t   runnumber = InputEvent()->GetRunNumber() ;

  TObjArray *matEMCAL=(TObjArray*)fOADBContainer->GetObject(runnumber,"EmcalMatrices");
  for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
    if(fGeoName=="EMCAL_FIRSTYEARV1" && mod>3)
      break;
    /*if(fESD)
      fGeom->SetMisalMatrix(fESD->GetEMCALMatrix(mod), mod);
      else*/
    // if(fVEvent->GetEMCALMatrix(mod))
    fGeomMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod);
    fGeom->SetMisalMatrix(fGeomMatrix[mod] , mod);
  }

  if(fESD){
    TList *l = fESD->GetList();
    if(1){//fDebug){
      for(int nk=0;nk<l->GetEntries();nk++){
	  TObject *obj = (TObject*)l->At(nk);
	  TString oname = obj->GetName();
	  //if(oname.Contains("lus"))
	  printf("Object %d has a clus array named %s +++++++++\n",nk,oname.Data());
	}
    }
    fESDClusters =  dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
    fESDCells = fESD->GetEMCALCells();
    if(fDebug)
      printf("ESD cluster mult= %d\n",fESDClusters->GetEntriesFast());
    if(fESDClusters->GetEntriesFast()<1){
      printf("ERROR: There are no EMCAL clusters in this event\n");
      return;
    }
  }
  else if(fAOD){
    fAODClusters = dynamic_cast<TClonesArray*>(fAOD->GetCaloClusters());
    fAODCells = fAOD->GetEMCALCells();
    if(fDebug)
      printf("AOD cluster mult= %d\n",fAODClusters->GetEntriesFast());
    if(fAODClusters->GetEntriesFast()<1){
      printf("ERROR: There are no EMCAL clusters in this event\n");
      return;
    }
  }
  FillClusterHists();
  if(fESD)
    FillTrackHists();

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::VZEROEventPlane()
{// Calculate the V0 Event Plane

      if (fEventPlane->GetQVector()) { 
        fEPTPC = TVector2::Phi_0_2pi(fEventPlane->GetQVector()->Phi())/2.0; //if(fEPTPC>TMath::Pi()) {fEPTPC-=TMath::Pi();} 
      } else { fEPTPC = -999.; }

      if (fEventPlane->GetQsub1()&&fEventPlane->GetQsub2()) {
        fEPTPCResolution = TMath::Cos(2.0*(fEventPlane->GetQsub1()->Phi()/2.0-fEventPlane->GetQsub2()->Phi()/2.0)); }
      else { fEPTPCResolution = -1; }

      fEPV0  = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0",  fESD)); if(fEPV0>TMath::Pi()) {fEPV0-=TMath::Pi();}
      fEPV0A = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0A", fESD)); if(fEPV0A>TMath::Pi()) {fEPV0A-=TMath::Pi();}
      fEPV0C = TVector2::Phi_0_2pi(fEventPlane->GetEventplane("V0C", fESD)); if(fEPV0C>TMath::Pi()) {fEPV0C-=TMath::Pi();}

      Double_t qx=0, qy=0, qxr=0, qyr=0;
      fEPV0Ar = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 4, 5, 2, qxr, qyr)); if(fEPV0Ar>TMath::Pi()) {fEPV0Ar-=TMath::Pi();}
      fEPV0Cr = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 2, 3, 2, qx,  qy)); if(fEPV0Cr>TMath::Pi()) {fEPV0Cr-=TMath::Pi();}
      qxr += qx; qyr += qy;
      fEPV0r   = TVector2::Phi_0_2pi(TMath::ATan2(qyr,qxr))/2.0; 

      fEPV0A4r = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 4, 2, qx, qy)); if(fEPV0A4r>TMath::Pi()) {fEPV0A4r-=TMath::Pi();}
      fEPV0A5r = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 5, 2, qx, qy)); if(fEPV0A5r>TMath::Pi()) {fEPV0A5r-=TMath::Pi();}
      fEPV0A6r = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 6, 2, qx, qy)); if(fEPV0A6r>TMath::Pi()) {fEPV0A6r-=TMath::Pi();}
      fEPV0A7r = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 7, 2, qx, qy)); if(fEPV0A7r>TMath::Pi()) {fEPV0A7r-=TMath::Pi();}
      fEPV0C0r = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 0, 2, qx, qy)); if(fEPV0C0r>TMath::Pi()) {fEPV0C0r-=TMath::Pi();}
      fEPV0C1r = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 1, 2, qx, qy)); if(fEPV0C1r>TMath::Pi()) {fEPV0C1r-=TMath::Pi();}
      fEPV0C2r = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 2, 2, qx, qy)); if(fEPV0C2r>TMath::Pi()) {fEPV0C2r-=TMath::Pi();}
      fEPV0C3r = TVector2::Phi_0_2pi(fEventPlane->CalculateVZEROEventPlane(fESD, 3, 2, qx, qy)); if(fEPV0C3r>TMath::Pi()) {fEPV0C3r-=TMath::Pi();}

      fHistEPTPC->Fill(fCentralityV0M,  fEPTPC); 
      if(fEPTPCResolution!=-1) { fHistEPTPCResolution->Fill(fCentralityV0M, fEPTPCResolution); }
      fHistEPV0->Fill(fCentralityV0M, fEPV0);
      fHistEPV0A->Fill(fCentralityV0M, fEPV0A);
      fHistEPV0C->Fill(fCentralityV0M, fEPV0C);
      fHistEPV0Ar->Fill(fCentralityV0M, fEPV0Ar);
      fHistEPV0Cr->Fill(fCentralityV0M, fEPV0Cr);
      fHistEPV0r->Fill(fCentralityV0M, fEPV0r);
      fHistEPV0A4r->Fill(fCentralityV0M, fEPV0A4r);
      fHistEPV0A7r->Fill(fCentralityV0M, fEPV0A7r);
      fHistEPV0C0r->Fill(fCentralityV0M, fEPV0C0r);
      fHistEPV0C3r->Fill(fCentralityV0M, fEPV0C3r);

      fHistAllcentV0->Fill(fEPV0);
      fHistAllcentV0r->Fill(fEPV0r);
      fHistAllcentV0A->Fill(fEPV0A);
      fHistAllcentV0C->Fill(fEPV0C);  
      fHistAllcentTPC->Fill(fEPTPC);

      fHistdifV0A_V0C0r->Fill(fCentralityV0M, TMath::Cos(2.0*(fEPV0A - fEPV0C0r)));
      fHistdifV0A_V0C3r->Fill(fCentralityV0M, TMath::Cos(2.0*(fEPV0A - fEPV0C3r)));
      fHistdifV0C0r_V0C3r->Fill(fCentralityV0M, TMath::Cos(2.0*(fEPV0C0r - fEPV0C3r)));
      fHistdifV0C_V0A4r->Fill(fCentralityV0M, TMath::Cos(2.0*(fEPV0C - fEPV0A4r)));
      fHistdifV0C_V0A7r->Fill(fCentralityV0M, TMath::Cos(2.0*(fEPV0C - fEPV0A7r)));
      fHistdifV0A4r_V0A7r->Fill(fCentralityV0M, TMath::Cos(2.0*(fEPV0A4r - fEPV0A7r)));   
      fHistdifV0Ar_V0Cr->Fill(fCentralityV0M, TMath::Cos(2.0*(fEPV0Ar - fEPV0Cr)));

}  

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::FillClusterHists()
{// Fill cluster histograms.

  Float_t pos[3] ; 

  TObjArray *clusters = fESDClusters;
  if (!clusters)
    clusters = fAODClusters;
  if (!clusters)
    return;

  const Int_t nclus = clusters->GetEntries();

  if(nclus==0)
    return;

  for(Int_t iclus=0;iclus<nclus;iclus++){
    AliVCluster *clus = (AliVCluster *) clusters->At(iclus); // retrieve cluster from esd
    if(!clus)
      continue;
    if(!clus->IsEMCAL()){ 
      printf("ERROR: Current Cluster is not an EMCAL Cluster\n");
      continue;
    }
    clus->GetPosition(pos);
    TVector3 vpos(pos[0],pos[1],pos[2]);
    Double_t Transverse_Energy = ((clus->E())/ (TMath::CosH(vpos.Eta())));    
    Double_t dphiV0 = TVector2::Phi_0_2pi(vpos.Phi()-fEPV0r); if(dphiV0>TMath::Pi()) {dphiV0-=TMath::Pi();}
    Double_t dphiV0A = TVector2::Phi_0_2pi(vpos.Phi()-fEPV0A); if(dphiV0A>TMath::Pi()) {dphiV0A-=TMath::Pi();}
    Double_t dphiV0C = TVector2::Phi_0_2pi(vpos.Phi()-fEPV0C); if(dphiV0C>TMath::Pi()) {dphiV0C-=TMath::Pi();}
    Double_t dphiTPC = TVector2::Phi_0_2pi(vpos.Phi()-fEPTPC); if(dphiTPC>TMath::Pi()) {dphiTPC-=TMath::Pi();}

    Double_t DataV0[5];
    DataV0[0] = Transverse_Energy; 	  
    DataV0[1] = clus->GetM02(); 
    DataV0[2] = fCentralityV0M;
    DataV0[3] = dphiV0;
    DataV0[4] = TMath::Cos(2.0*(dphiV0));
    fClusterPbV0->Fill(DataV0);

    Double_t DataV0A[5];
    DataV0A[0] = Transverse_Energy; 	  
    DataV0A[1] = clus->GetM02(); 
    DataV0A[2] = fCentralityV0M;
    DataV0A[3] = dphiV0A;
    DataV0A[4] = TMath::Cos(2.0*(dphiV0A));
    fClusterPbV0A->Fill(DataV0A);

    Double_t DataV0C[5];
    DataV0C[0] = Transverse_Energy; 	  
    DataV0C[1] = clus->GetM02(); 
    DataV0C[2] = fCentralityV0M;
    DataV0C[3] = dphiV0C;
    DataV0C[4] = TMath::Cos(2.0*(dphiV0C));
    fClusterPbV0C->Fill(DataV0C);

    Double_t DataTPC[5];
    DataTPC[0] = Transverse_Energy; 	  
    DataTPC[1] = clus->GetM02(); 
    DataTPC[2] = fCentralityV0M;
    DataTPC[3] = dphiTPC;
    DataTPC[4] = TMath::Cos(2.0*(dphiTPC));
    fClusterPbTPC->Fill(DataTPC); 

    fHistClusterE->Fill(clus->E());
    fHistClusterEt->Fill(Transverse_Energy);
    fHistClusterM02->Fill(clus->GetM02());
    fHistClusterN->Fill(clus->GetNCells());
    fHistClusterPhi->Fill(vpos.Phi());
    fHistClusterEta->Fill(vpos.Eta());  
    fHistClusterPhiEta->Fill(vpos.Eta(),vpos.Phi()); 
    fHistClusterEN->Fill(clus->E(),clus->GetNCells()); 
    fHistClusterEM02->Fill(clus->E(),clus->GetM02()); 
    fHistClusterEtN->Fill(Transverse_Energy,clus->GetNCells()); 
    fHistClusterEtM02->Fill(Transverse_Energy,clus->GetM02()); 
    fHistClusterdphiV0->Fill(dphiV0);

  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::FillTrackHists()
{// Fill track histograms.
  if(!fESD)
    return;
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
}
  
//________________________________________________________________________
void AliAnalysisTaskEMCALPi0V2ShSh::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations
  // Called once at the end of the query

 /* fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutputList) { Printf("ERROR: could not retrieve TList fOutputList"); return; }
      
  fHistTrackPt = dynamic_cast<TH1F*> (fOutputList->FindObject("fHistTrackPt"));
  if (!fHistTrackPt) { Printf("ERROR: could not retrieve fHistTrackPt"); return;}
  fHistTrackEta = dynamic_cast<TH1F*> (fOutputList->FindObject("fHistTrackEta"));
  if (!fHistTrackEta) { Printf("ERROR: could not retrieve fHistTrackEta"); return;}
  fHistTrackPhi = dynamic_cast<TH1F*> (fOutputList->FindObject("fHistTrackPhi"));
  if (!fHistTrackPhi) { Printf("ERROR: could not retrieve fHistTrackPhi"); return;}
  fHistClusterEta = dynamic_cast<TH1F*> (fOutputList->FindObject("fHistClusterEta"));
  if (!fHistClusterEta) { Printf("ERROR: could not retrieve fHistClusterEta"); return;}
  fHistClusterPhi = dynamic_cast<TH1F*> (fOutputList->FindObject("fHistClusterPhi"));
  if (!fHistClusterPhi) { Printf("ERROR: could not retrieve fHistClusterPhi"); return;}
  fHistClusterE = dynamic_cast<TH1F*> (fOutputList->FindObject("fHistClusterE"));
  if (!fHistClusterE) { Printf("ERROR: could not retrieve fHistClusterE"); return;}
  fHistClusterN = dynamic_cast<TH1F*> (fOutputList->FindObject("fHistClusterN"));
  if (!fHistClusterN) { Printf("ERROR: could not retrieve fHistClusterN"); return;}
  fHistClusterM02 = dynamic_cast<TH1F*> (fOutputList->FindObject("fHistClusterM02"));
  if (!fHistClusterM02) { Printf("ERROR: could not retrieve fHistClusterM02"); return;}
      
    // Get the physics selection histograms with the selection statistics
    //AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    //AliESDInputHandler *inputH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
    //TH2F *histStat = (TH2F*)inputH->GetStatistics();
   
   
  TCanvas *canvas1 = new TCanvas("canvas1","Track P_{T} & #eta & #phi",10,10,1020,510);
  canvas1->Divide(3,1);
  canvas1->cd(1)->SetLogy();
  fHistTrackPt->DrawCopy("");
  canvas1->cd(2);
  fHistTrackEta->DrawCopy("");
  canvas1->cd(3);
  fHistTrackPhi->DrawCopy("");
 
  TCanvas *canvas2 = new TCanvas("canvas2","Cluster #eta & #phi",10,10,1020,510);
  canvas2->Divide(2,1);
  canvas2->cd(1);
  fHistClusterEta->DrawCopy("");
  canvas2->cd(2);
  fHistClusterPhi->DrawCopy("");

  TCanvas *canvas3 = new TCanvas("canvas3","Cluster E & N & M02",10,10,1020,510);
  canvas3->Divide(3,1);
  canvas3->cd(1)->SetLogy();
  fHistClusterE->DrawCopy("");
  canvas3->cd(2)->SetLogy();
  fHistClusterN->DrawCopy("");
  canvas3->cd(3)->SetLogy();
  fHistClusterM02->DrawCopy("");

  TCanvas *canvas4 = new TCanvas("canvas4","Track #phi vs #eta &  Cluster #phi vs #eta",10,10,1020,510);
  canvas4->Divide(2,1);
  canvas4->cd(1)->SetLogz();
  fHistTrackPhiEta->DrawCopy("colorz");
  canvas4->cd(2)->SetLogz();
  fHistClusterPhiEta->DrawCopy("colorz");

  TCanvas *canvas5 = new TCanvas("canvas5","Cluster E vs N &  E vs M02",10,10,1020,510);
  canvas5->Divide(2,1);
  canvas5->cd(1)->SetLogz();
  fHistClusterEN->DrawCopy("colorz");
  canvas5->cd(2)->SetLogz();
  fHistClusterEM02->DrawCopy("colorz");

  TCanvas *canvas6 = new TCanvas("canvas6","Cluster Et & Et vs N & Et vs M02",10,10,1020,510);
  canvas6->cd(1)->SetLogy();
  fHistClusterEt->DrawCopy("");
  canvas6->cd(2)->SetLogz();
  fHistClusterEtN->DrawCopy("colorz");
  canvas6->cd(3)->SetLogz();
  fHistClusterEtM02->DrawCopy("colorz");

  canvas1->SaveAs("lhc11h_2_Track_PT_Eta_Phi.jpg");
  canvas2->SaveAs("lhc11h_2_Cluster_Eta_Phi.jpg");
  canvas3->SaveAs("lhc11h_2_Cluster_E_N_M02.jpg");
  canvas4->SaveAs("lhc11h_2_Track_PhivsEta_Cluster_PhivsEta.jpg");
  canvas5->SaveAs("lhc11h_2_Cluster_EvsN_EvsM02.jpg");
  canvas6->SaveAs("lhc11h_2_Cluster_Et_EtvsN_EtvsM02.jpg"); */
}
