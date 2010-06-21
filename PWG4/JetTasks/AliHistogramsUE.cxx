/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: A.Abrahantes, E.Lopez, S.Vallero                               *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id:$ */
////////////////////////////////////////////////
//--------------------------------------------- 
// Class  to handle histograms for UE analysis
//---------------------------------------------
////////////////////////////////////////////////
#include <TROOT.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliHistogramsUE.h"
#include "AliAnalyseUE.h"
#include "AliLog.h"

ClassImp( AliHistogramsUE)

//____________________________________________________________________
AliHistogramsUE:: AliHistogramsUE():
TObject(),
fBinsPtInHist(0),
fMinJetPtInHist(0.),
fMaxJetPtInHist(0.),
fTrackEtaCut(0.),
fListOfHistos(0x0),  
fhNJets(0x0),
fhEleadingPt(0x0),
fhMinRegPtDist(0x0),
fhRegionMultMin(0x0),
fhMinRegAvePt(0x0), 
fhMinRegSumPt(0x0),            
fhMinRegMaxPtPart(0x0),
fhMinRegSumPtvsMult(0x0),
fhdNdEtaPhiDist(0x0),        
fhFullRegPartPtDistVsEt(0x0), 
fhTransRegPartPtDistVsEt(0x0),
fhRegionSumPtMaxVsEt(0x0),
fhRegionMultMax(0x0),         
fhRegionMultMaxVsEt(0x0),     
fhRegionSumPtMinVsEt(0x0), //fhRegionMultMin(0x0),         
fhRegionMultMinVsEt(0x0),
fhRegionAveSumPtVsEt(0x0),
fhRegionDiffSumPtVsEt(0x0),
fhRegionAvePartPtMaxVsEt(0x0),
fhRegionAvePartPtMinVsEt(0x0),
fhRegionMaxPartPtMaxVsEt(0x0),
fhRegForwardMult(0x0),
fhRegForwardSumPtvsMult(0x0),
fhRegBackwardMult(0x0),
fhRegBackwardSumPtvsMult(0x0),
fhRegForwardPartPtDistVsEt(0x0),
fhRegBackwardPartPtDistVsEt(0x0),
fhRegTransMult(0x0),
fhRegTransSumPtVsMult(0x0),
fhMinRegSumPtJetPtBin(0x0),
fhMaxRegSumPtJetPtBin(0x0),
fhVertexMult(0x0),
fh1Xsec(0x0),
fh1Trials(0x0)
{
  // Default constructor

}

//____________________________________________________________________
AliHistogramsUE:: AliHistogramsUE(TList *list):
TObject(),
fBinsPtInHist(0),
fMinJetPtInHist(0.),
fMaxJetPtInHist(0.),
fTrackEtaCut(0.),
fListOfHistos(0x0),  
fhNJets(0x0),
fhEleadingPt(0x0),
fhMinRegPtDist(0x0),
fhRegionMultMin(0x0),
fhMinRegAvePt(0x0), 
fhMinRegSumPt(0x0),            
fhMinRegMaxPtPart(0x0),
fhMinRegSumPtvsMult(0x0),
fhdNdEtaPhiDist(0x0),        
fhFullRegPartPtDistVsEt(0x0), 
fhTransRegPartPtDistVsEt(0x0),
fhRegionSumPtMaxVsEt(0x0),
fhRegionMultMax(0x0),         
fhRegionMultMaxVsEt(0x0),     
fhRegionSumPtMinVsEt(0x0), //fhRegionMultMin(0x0),         
fhRegionMultMinVsEt(0x0),
fhRegionAveSumPtVsEt(0x0),
fhRegionDiffSumPtVsEt(0x0),
fhRegionAvePartPtMaxVsEt(0x0),
fhRegionAvePartPtMinVsEt(0x0),
fhRegionMaxPartPtMaxVsEt(0x0),
fhRegForwardMult(0x0),
fhRegForwardSumPtvsMult(0x0),
fhRegBackwardMult(0x0),
fhRegBackwardSumPtvsMult(0x0),
fhRegForwardPartPtDistVsEt(0x0),
fhRegBackwardPartPtDistVsEt(0x0),
fhRegTransMult(0x0),
fhRegTransSumPtVsMult(0x0),
fhMinRegSumPtJetPtBin(0x0),
fhMaxRegSumPtJetPtBin(0x0),
fhVertexMult(0x0),
fh1Xsec(0x0),
fh1Trials(0x0)
{
  // Constructor, initialize members from list
	fhNJets = (TH1F*)list->FindObject("hNJets"); 
	fhEleadingPt = (TH1F*)list->FindObject("hEleadingPt");
	fhMinRegPtDist = (TH1F*)list->FindObject("hMinRegPtDist");
	fhRegionMultMin = (TH1F*)list->FindObject("hRegionMultMin");
	fhMinRegAvePt = (TH1F*)list->FindObject("hMinRegAvePt");
	fhMinRegSumPt = (TH1F*)list->FindObject("hMinRegSumPt");
	fhMinRegMaxPtPart = (TH1F*)list->FindObject("hMinRegMaxPtPart");
	fhMinRegSumPtvsMult = (TH1F*)list->FindObject("hMinRegSumPtvsMult");
	fhdNdEtaPhiDist = (TH2F*)list->FindObject("hdNdEtaPhiDist");  
	fhFullRegPartPtDistVsEt = (TH2F*)list->FindObject("hFullRegPartPtDistVsEt");
	fhTransRegPartPtDistVsEt = (TH2F*)list->FindObject("hTransRegPartPtDistVsEt");
	fhRegionSumPtMaxVsEt = (TH1F*)list->FindObject("hRegionSumPtMaxVsEt"); 
	fhRegionMultMax = (TH1I*)list->FindObject("hRegionMultMax");
	fhRegionMultMaxVsEt = (TH1F*)list->FindObject("hRegionMultMaxVsEt");
	fhRegionSumPtMinVsEt = (TH1F*)list->FindObject("hRegionSumPtMinVsEt");
	fhRegionMultMinVsEt = (TH1F*)list->FindObject("hRegionMultMinVsEt");
	fhRegionAveSumPtVsEt = (TH1F*)list->FindObject("hRegionAveSumPtVsEt");
	fhRegionDiffSumPtVsEt = (TH1F*)list->FindObject("hRegionDiffSumPtVsEt");
	fhRegionAvePartPtMaxVsEt = (TH1F*)list->FindObject("hRegionAvePartPtMaxVsEt");
	fhRegionAvePartPtMinVsEt = (TH1F*)list->FindObject("hRegionAvePartPtMinVsEt");
	fhRegionMaxPartPtMaxVsEt = (TH1F*)list->FindObject("hRegionMaxPartPtMaxVsEt");
	fhRegForwardMult = (TH2F*)list->FindObject("hRegForwardMult"); 
	fhRegForwardSumPtvsMult = (TH2F*)list->FindObject("hRegForwardSumPtvsMult");
	fhRegBackwardMult = (TH2F*)list->FindObject("hRegBackwardMult");
	fhRegBackwardSumPtvsMult = (TH2F*)list->FindObject("hRegBackwardSumPtvsMult");
	fhRegForwardPartPtDistVsEt = (TH2F*)list->FindObject("hRegForwardPartPtDistVsEt");
	fhRegBackwardPartPtDistVsEt = (TH2F*)list->FindObject("hRegBackwardPartPtDistVsEt");
	fhRegTransMult = (TH2F*)list->FindObject("hRegTransMult");
	fhRegTransSumPtVsMult = (TH2F*)list->FindObject("hRegTransSumPtVsMult");
	fhMinRegSumPtJetPtBin = (TH2F*)list->FindObject("hMinRegSumPtJetPtBin");
	fhMaxRegSumPtJetPtBin = (TH2F*)list->FindObject("hMaxRegSumPtJetPtBin");
	fhVertexMult = (TH1F*)list->FindObject("hVertexMult");
	fh1Xsec = (TProfile*)list->FindObject("h1Xsec");
	fh1Trials = (TH1F*)list->FindObject("h1Trials"); 

}
//____________________________________________________________________
AliHistogramsUE:: AliHistogramsUE(const AliHistogramsUE & original):
TObject(),
fBinsPtInHist(original.fBinsPtInHist),
fMinJetPtInHist(original.fMinJetPtInHist),
fMaxJetPtInHist(original.fMaxJetPtInHist),
fTrackEtaCut(original.fTrackEtaCut),
fListOfHistos(original.fListOfHistos),  
fhNJets(original.fhNJets),
fhEleadingPt(original.fhEleadingPt),
fhMinRegPtDist(original.fhMinRegPtDist),
fhRegionMultMin(original.fhRegionMultMin),
fhMinRegAvePt(original.fhMinRegAvePt), 
fhMinRegSumPt(original.fhMinRegSumPt),            
fhMinRegMaxPtPart(original.fhMinRegMaxPtPart),
fhMinRegSumPtvsMult(original.fhMinRegSumPtvsMult),
fhdNdEtaPhiDist(original.fhdNdEtaPhiDist),        
fhFullRegPartPtDistVsEt(original.fhFullRegPartPtDistVsEt), 
fhTransRegPartPtDistVsEt(original.fhTransRegPartPtDistVsEt),
fhRegionSumPtMaxVsEt(original.fhRegionSumPtMaxVsEt),
fhRegionMultMax(original.fhRegionMultMax),         
fhRegionMultMaxVsEt(original.fhRegionMultMaxVsEt),     
fhRegionSumPtMinVsEt(original.fhRegionSumPtMinVsEt),         
fhRegionMultMinVsEt(original.fhRegionMultMinVsEt),
fhRegionAveSumPtVsEt(original.fhRegionAveSumPtVsEt),
fhRegionDiffSumPtVsEt(original.fhRegionDiffSumPtVsEt),
fhRegionAvePartPtMaxVsEt(original.fhRegionAvePartPtMaxVsEt),
fhRegionAvePartPtMinVsEt(original.fhRegionAvePartPtMinVsEt),
fhRegionMaxPartPtMaxVsEt(original.fhRegionMaxPartPtMaxVsEt),
fhRegForwardMult(original.fhRegForwardMult),
fhRegForwardSumPtvsMult(original.fhRegForwardSumPtvsMult),
fhRegBackwardMult(original.fhRegBackwardMult),
fhRegBackwardSumPtvsMult(original.fhRegBackwardSumPtvsMult),
fhRegForwardPartPtDistVsEt(original.fhRegForwardPartPtDistVsEt),
fhRegBackwardPartPtDistVsEt(original.fhRegBackwardPartPtDistVsEt),
fhRegTransMult(original.fhRegTransMult),
fhRegTransSumPtVsMult(original.fhRegTransSumPtVsMult),
fhMinRegSumPtJetPtBin(original.fhMinRegSumPtJetPtBin),
fhMaxRegSumPtJetPtBin(original.fhMaxRegSumPtJetPtBin),
fhVertexMult(original.fhVertexMult),
fh1Xsec(original.fh1Xsec),
fh1Trials(original.fh1Trials)
{

  // Copy constructor

}


//______________________________________________________________
AliHistogramsUE & AliHistogramsUE::operator = (const AliHistogramsUE & /*source*/)
{

  // assignment operator
  return *this;

}

//______________________________________________________________
TObjArray* AliHistogramsUE::CreateCanvas(const Int_t ncanv){
	
	// Create canvas for plotting 
        printf("Creating %d canvas ... \n",ncanv);
        TObjArray *arr=new TObjArray;
	TString name;
        for(Int_t i=0;i<ncanv;i++ ){
        name=Form("Canvas %d",i);
        TCanvas *c=new TCanvas(name,name);
        c->SetFillColor(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	arr->Add(c);
        }
return arr;
}

//____________________________________________________________________
TList*  AliHistogramsUE::CreateHistos(Int_t bins, Double_t min, Double_t max, Double_t etacut)
{

  // Create all histograms necessary for UE analysis 
  fBinsPtInHist = bins;
  fMinJetPtInHist = min;
  fMaxJetPtInHist = max;
  fTrackEtaCut= etacut;

  fListOfHistos = new TList();

  //Number of reconstructed clusters  
  fhNJets = new TH1F("hNJets", "Number of clusters",  20, 0, 20);
  fhNJets->SetXTitle("Number of reconstructed clusters");
  fhNJets->SetYTitle("#");
  fhNJets->Sumw2();
  fListOfHistos->Add( fhNJets );                 // At(0) 
  
  //pT distribution of leading clusters
  fhEleadingPt  = new TH1F("hEleadingPt",   "Leading cluster p_{T}",  bins, min,   max);
  fhEleadingPt->SetXTitle("p_{T} of  cluster (GeV/c)");
  fhEleadingPt->SetYTitle("1/N_{ev} dN/dp_{T} (|#eta|<0.5)");
  fhEleadingPt->Sumw2();
  fListOfHistos->Add( fhEleadingPt );            // At(1)
  
  //Track pT distribution in MIN zone
  fhMinRegPtDist = new TH1F("hMinRegPtDist",   "p_{T} distribution in MIN zone",  50,0.,20.);
  fhMinRegPtDist->SetXTitle("Track p_{T} (GeV/c)");
  fhMinRegPtDist->SetYTitle("dN/dp_{T}");
  fhMinRegPtDist->Sumw2();
  fListOfHistos->Add( fhMinRegPtDist );          // At(2)
  
  //Multiplicity in MIN zone
  fhRegionMultMin = new TH1F("hRegionMultMin",      "N_{ch}^{90, min}",  21, -0.5,   20.5);
  fhRegionMultMin->SetXTitle("N_{ch tracks}");
  fhRegionMultMin->Sumw2();
  fListOfHistos->Add( fhRegionMultMin );         // At(3)            
  
  //Average pT in MIN region
  fhMinRegAvePt = new TH1F("hMinRegAvePt", "#LTp_{T}#GT",  50, 0.,   20.);
  fhMinRegAvePt->SetXTitle("p_{T} (GeV/c)");
  fhMinRegAvePt->Sumw2();
  fListOfHistos->Add( fhMinRegAvePt );           // At(4)
  
  //Sum pT in MIN region
  fhMinRegSumPt = new TH1F("hMinRegSumPt", "#Sigma p_{T} ",  50, 0.,   20.);
  fhMinRegSumPt->SetYTitle("Ed^{3}N_{tracks}/dp^{3} (c^{3}/GeV^{2})");  
  fhMinRegSumPt->SetXTitle("#Sigma p_{T} (GeV/c)");
  fhMinRegSumPt->Sumw2();
  fListOfHistos->Add( fhMinRegSumPt );           // At(5)
  
  //Track with maximum pT in MIN region
  fhMinRegMaxPtPart = new TH1F("hMinRegMaxPtPart", "max(p_{T})|_{event} ",  50, 0.,   20.);
  fhMinRegMaxPtPart->SetYTitle("Ed^{3}N_{tracks}/dp^{3} (c^{3}/GeV^{2})");  
  fhMinRegMaxPtPart->SetXTitle("p_{T} (GeV/c)");
  fhMinRegMaxPtPart->Sumw2();
  fListOfHistos->Add( fhMinRegMaxPtPart );       // At(6)
  
  //Sum pT vs. multiplicity in MIN region
  fhMinRegSumPtvsMult = new TH1F("hMinRegSumPtvsMult", "#Sigma p_{T} vs. Multiplicity ",  21, -0.5,   20.5);
  fhMinRegSumPtvsMult->SetYTitle("#Sigma p_{T} (GeV/c)");  
  fhMinRegSumPtvsMult->SetXTitle("N_{charge}");
  fhMinRegSumPtvsMult->Sumw2();
  fListOfHistos->Add( fhMinRegSumPtvsMult );     // At(7);
  
  //Phi correlation track-cluster vs. leading cluster pT 
  fhdNdEtaPhiDist  = new TH2F("hdNdEtaPhiDist", Form("Charge particle density |#eta|<%3.1f vs #Delta#phi", fTrackEtaCut),62, 0.,   2.*TMath::Pi(), bins, min, max);
  fhdNdEtaPhiDist->SetXTitle("#Delta#phi");
  fhdNdEtaPhiDist->SetYTitle("Leading cluster p_{T}");
  fhdNdEtaPhiDist->Sumw2();
  fListOfHistos->Add( fhdNdEtaPhiDist );        // At(8)
  
  //Can be used to get track pT distribution for different cluster pT bins (full region)
  fhFullRegPartPtDistVsEt = new TH2F("hFullRegPartPtDistVsEt", Form( "dN/dp_{T} |#eta|<%3.1f vs Leading cluster p_{T}", fTrackEtaCut),100,0.,50., bins, min, max);
  fhFullRegPartPtDistVsEt->SetYTitle("Leading cluster p_{T}");
  fhFullRegPartPtDistVsEt->SetXTitle("p_{T}");
  fhFullRegPartPtDistVsEt->Sumw2();
  fListOfHistos->Add( fhFullRegPartPtDistVsEt );  // At(9) 
  
  //Can be used to get part pT distribution for different cluster pT bins (transverse region)
  fhTransRegPartPtDistVsEt = new TH2F("hTransRegPartPtDistVsEt", Form( "dN/dp_{T} in tranvese regions |#eta|<%3.1f vs Leading cluster p_{T}", fTrackEtaCut),100,0.,50., bins, min,   max);
  fhTransRegPartPtDistVsEt->SetYTitle("Leading cluster p_{T}");
  fhTransRegPartPtDistVsEt->SetXTitle("p_{T}");
  fhTransRegPartPtDistVsEt->Sumw2();
  fListOfHistos->Add( fhTransRegPartPtDistVsEt );  // At(10)
  
  //Sum pT in MAX region vs. leading-cluster pT
  fhRegionSumPtMaxVsEt = new TH1F("hRegionSumPtMaxVsEt",  "P_{T}^{90, max} vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionSumPtMaxVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegionSumPtMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionSumPtMaxVsEt );      // At(11)
  

  //Sum pT in MIN region vs. leading-cluster pT
  fhRegionSumPtMinVsEt = new TH1F("hRegionSumPtMinVsEt",   "P_{T}^{90, min} vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionSumPtMinVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegionSumPtMinVsEt->Sumw2();
  fListOfHistos->Add( fhRegionSumPtMinVsEt );      // At(12)
 
  //Multiplicity in MAX region
  fhRegionMultMax = new TH1I("hRegionMultMax",      "N_{ch}^{90, max}",  21, -0.5,   20.5);
  fhRegionMultMax->SetXTitle("N_{ch tracks}");
  fhRegionMultMax->Sumw2();
  fListOfHistos->Add( fhRegionMultMax );           // At(13)
  
  //Multiplicity in MAX region vs. leading-cluster pT
  fhRegionMultMaxVsEt = new TH1F("hRegionMultMaxVsEt",  "N_{ch}^{90, max} vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionMultMaxVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegionMultMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionMultMaxVsEt );       // At(14)
  
  //Multiplicity in MIN region vs. leading-cluster pT
  fhRegionMultMinVsEt = new TH1F("hRegionMultMinVsEt",  "N_{ch}^{90, min} vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionMultMinVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegionMultMinVsEt->Sumw2();
  fListOfHistos->Add( fhRegionMultMinVsEt );      // At(15)
 
  //Average sum pT in TRANSVERSE(MIN+MAX) region vs. leading-cluster pT 
  fhRegionAveSumPtVsEt = new TH1F("hRegionAveSumPtVsEt", "(P_{T}^{90, max} + P_{T}^{90, min})/2 vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionAveSumPtVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegionAveSumPtVsEt->Sumw2();
  fListOfHistos->Add( fhRegionAveSumPtVsEt );     // At(16)
  
  //Difference sum pT (MAX-MIN) vs.  leading-cluster pT
  fhRegionDiffSumPtVsEt= new TH1F("hRegionDiffSumPtVsEt", "(P_{T}^{90, max} - P_{T}^{90, min}) vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionDiffSumPtVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegionDiffSumPtVsEt->Sumw2();
  fListOfHistos->Add( fhRegionDiffSumPtVsEt );    // At(17)
  
  //Average track pT in MAX region vs. leading-cluster pT
  fhRegionAvePartPtMaxVsEt = new TH1F("hRegionAvePartPtMaxVsEt", "#LTp_{T}#GT^{90, max} vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionAvePartPtMaxVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegionAvePartPtMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionAvePartPtMaxVsEt );  // At(18)
  
  //Average track pT in MIN region vs. leading-cluster pT
  fhRegionAvePartPtMinVsEt = new TH1F("hRegionAvePartPtMinVsEt", "#LTp_{T}#GT^{90, min} vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionAvePartPtMinVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegionAvePartPtMinVsEt->Sumw2();
  fListOfHistos->Add( fhRegionAvePartPtMinVsEt );   // At(19)
  
  //Maximum track pT in MAX region vs. leading-cluster pT
  fhRegionMaxPartPtMaxVsEt = new TH1F("hRegionMaxPartPtMaxVsEt", "max(p_{T})^{90} vs Leading cluster p_{T}",  bins, min,   max);
  fhRegionMaxPartPtMaxVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegionMaxPartPtMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionMaxPartPtMaxVsEt );    // At(20)
  
  //Multiplicity in FORWARD region
  fhRegForwardMult = new TH2F("hRegForwardMult", "N_{ch}^{forward}", bins, min, max, 21, -0.5,   20.5);
  fhRegForwardMult->SetXTitle("N_{ch tracks}");
  fhRegForwardMult->Sumw2();
  fListOfHistos->Add( fhRegForwardMult );           // At(25)
  
  //Sum pT in FORWARD region vs. multiplicity
  fhRegForwardSumPtvsMult = new TH2F("hRegForwardSumPtvsMult", "Forward #Sigma p_{T} vs. Multiplicity ", bins, min, max, 21, -0.5,   20.5);
  fhRegForwardSumPtvsMult->SetYTitle("#Sigma p_{T} (GeV/c)");  
  fhRegForwardSumPtvsMult->SetXTitle("N_{charge}");
  fhRegForwardSumPtvsMult->Sumw2();
  fListOfHistos->Add( fhRegForwardSumPtvsMult );     // At(26);
  
  //Multiplicity in BACKWARD region
  fhRegBackwardMult = new TH2F("hRegBackwardMult", "N_{ch}^{backward}", bins, min, max, 21, -0.5,   20.5);
  fhRegBackwardMult->SetXTitle("N_{ch tracks}");
  fhRegBackwardMult->Sumw2();
  fListOfHistos->Add( fhRegBackwardMult );           // At(27)
 
  //Sum pT in BACKWARD region vs. multiplicity
  fhRegBackwardSumPtvsMult = new TH2F("hRegBackwardSumPtvsMult", "Backward #Sigma p_{T} vs. Multiplicity ", bins, min, max, 21, -0.5,   20.5);
  fhRegBackwardSumPtvsMult->SetYTitle("#Sigma p_{T} (GeV/c)");  
  fhRegBackwardSumPtvsMult->SetXTitle("N_{charge}");
  fhRegBackwardSumPtvsMult->Sumw2();
  fListOfHistos->Add( fhRegBackwardSumPtvsMult );     // At(28);
  
  //Track pT distribution in FORWARD region vs. leading-cluster pT 
  fhRegForwardPartPtDistVsEt = new TH2F("hRegForwardPartPtDistVsEt", Form( "dN/dP_{T} |#eta|<%3.1f vs Leading cluster p_{T}", fTrackEtaCut), 100,0.,50., bins, min, max);
  fhRegForwardPartPtDistVsEt->SetYTitle("Leading cluster p_{T}");
  fhRegForwardPartPtDistVsEt->SetXTitle("p_{T} (GeV/c)");
  fhRegForwardPartPtDistVsEt->Sumw2();
  fListOfHistos->Add( fhRegForwardPartPtDistVsEt );  // At(29) 
  
  //Track pT distribution in BACKWARD region vs. leading-cluster pT 
  fhRegBackwardPartPtDistVsEt = new TH2F("hRegBackwardPartPtDistVsEt", Form( "dN/dP_{T} |#eta|<%3.1f vs Leading cluster p_{T}", fTrackEtaCut), 100,0.,50., bins, min, max);
  fhRegBackwardPartPtDistVsEt->SetYTitle("Leading cluster p_{T}");
  fhRegBackwardPartPtDistVsEt->SetXTitle("p_{T}");
  fhRegBackwardPartPtDistVsEt->Sumw2();
  fListOfHistos->Add( fhRegBackwardPartPtDistVsEt );  // At(30) 
  
  //Multiplicity in TRANSVERSE (MIN+MAX) region
  fhRegTransMult  = new TH2F("hRegTransMult", "N_{ch}^{transv}", bins, min, max, 21, -0.5,   20.5);
  fhRegTransMult->SetXTitle("N_{ch tracks}");
  fhRegTransMult->Sumw2();
  fListOfHistos->Add( fhRegTransMult );           // At(31)
  
  //Sum pT in TRANSVERSE (MIN+MAX) region vs. multiplicity
  fhRegTransSumPtVsMult = new TH2F("hRegTransSumPtVsMult", "Transverse #Sigma p_{T} vs. Multiplicity ",bins, min, max, 21, -0.5,   20.5);
  fhRegTransSumPtVsMult->SetYTitle("#Sigma p_{T} (GeV/c)");  
  fhRegTransSumPtVsMult->SetXTitle("N_{charge}");
  fhRegTransSumPtVsMult->Sumw2();
  fListOfHistos->Add( fhRegTransSumPtVsMult );     // At(32);
  
  //Sum pT in MIN region per cluster pT bin
  fhMinRegSumPtJetPtBin = new TH2F("hMinRegSumPtJetPtBin", "Transverse Min Reg #Sigma p_{T} per cluster pT bin",bins, min, max, 50, 0.,   20.);
  fhMinRegSumPtJetPtBin->SetXTitle("Leading cluster p_{T}");
  fhMinRegSumPtJetPtBin->Sumw2();
  fListOfHistos->Add( fhMinRegSumPtJetPtBin );           // At(33)
  
  //Sum pT in MAX region per cluster pT bin
  fhMaxRegSumPtJetPtBin = new TH2F("hMaxRegSumPtJetPtBin",      "Transverse Max Reg #Sigma p_{T} per cluster pT bin", bins, min, max, 50, 0.,   20.);
  fhMaxRegSumPtJetPtBin->SetXTitle("Leading cluster p_{T}");
  fhMaxRegSumPtJetPtBin->Sumw2();
  fListOfHistos->Add( fhMaxRegSumPtJetPtBin );           // At(34)
  
  //Multiplicity in main vertex
  fhVertexMult = new TH1F("hVertexMult",      "Multiplicity in Main Vertex", 81, -0.5 , 80.5);
  fhVertexMult->SetXTitle("Main Vertex Multiplicity");
  fhVertexMult->Sumw2();
  fListOfHistos->Add( fhVertexMult ); //At(35)
  
  fh1Xsec = new TProfile("h1Xsec","xsec from pyxsec.root",1,0,1); 
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fh1Xsec->Sumw2();
  fListOfHistos->Add( fh1Xsec );            //At(36)
  
  fh1Trials = new TH1F("h1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fh1Trials->Sumw2();
  fListOfHistos->Add( fh1Trials ); //At(37)

  return fListOfHistos;
}


//____________________________________________________________________
void AliHistogramsUE::DrawUE(Int_t debug){

    // To draw histograms at the end of task running 
    // Normalize histos to region area TODO: 
    // Normalization done at Analysis, taking into account 
    // area variations on per-event basis (cone case)
   
    //HIGH WARNING!!!!!: DO NOT SCALE ANY OF THE ORIGINAL HISTOGRAMS
    //MAKE A COPY, DRAW IT, And later sacale that copy. CAF Issue!!!!!
    
    Int_t binsPtInHist = fhEleadingPt->GetNbinsX();
    Double_t minJetPtInHist = fhEleadingPt->GetXaxis()->GetBinLowEdge(1);
    Double_t maxJetPtInHist = fhEleadingPt->GetXaxis()->GetBinUpEdge(binsPtInHist);
     
    //Sum pT
    TCanvas* c1 = new TCanvas("c1",Form("sumPt dist (%s)", GetTitle()),60,60,1100,700);
    c1->Divide(2,2);
    c1->cd(1);
    TH1F *h1r = new TH1F("hRegionEtvsSumPtMax" , "", binsPtInHist,  minJetPtInHist, maxJetPtInHist);
    //TH1F *h1r = new TH1F();
    h1r->Divide(fhRegionSumPtMaxVsEt,fhEleadingPt,1,1);
    //h1r->Scale( areafactor );
    h1r->SetMarkerStyle(20);
    h1r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h1r->SetYTitle("P_{T}^{90, max}");
    h1r->DrawCopy("p");
    
    c1->cd(2);
    TH1F *h2r = new TH1F("hRegionEtvsSumPtMin" , "", binsPtInHist,  minJetPtInHist, maxJetPtInHist);
    h2r->Divide(fhRegionSumPtMinVsEt,fhEleadingPt,1,1);
    //h2r->Scale( areafactor );
    h2r->SetMarkerStyle(20);
    h2r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h2r->SetYTitle("P_{T}^{90, min}");
    h2r->DrawCopy("p");
    
    c1->cd(3);
    TH1F *h4r = new TH1F("hRegionEtvsDiffPt" , "", binsPtInHist,  minJetPtInHist, maxJetPtInHist);
    //TH1F *h41r = new TH1F("hRegForwvsDiffPt" , "", fbinsPtInHist,  fMinJetPtInHist, fMaxJetPtInHist);
    //TH1F *h42r = new TH1F("hRegBackvsDiffPt" , "", fbinsPtInHist,  fMinJetPtInHist, fMaxJetPtInHist);
    //h41r->Divide(fhRegForwardSumPtVsEt,fhEleadingPt,1,1);
    //h42r->Divide(fhRegBackwardSumPtVsEt,fhEleadingPt,1,1);
    h4r->Divide(fhRegionAveSumPtVsEt,fhEleadingPt,1,1);
    //h4r->Scale(2.); // make average
    //h4r->Scale( areafactor );
    h4r->SetYTitle("#DeltaP_{T}^{90}");
    h4r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h4r->SetMarkerStyle(20);
    h4r->DrawCopy("p");

    c1->cd(4);
    TH1F *h5r = new TH1F("hRegionMultMaxVsEtleading",   "",  binsPtInHist, minJetPtInHist,   maxJetPtInHist);
    TH1F *h6r = new TH1F("hRegionMultMinVsEtleading",   "",  binsPtInHist, minJetPtInHist,   maxJetPtInHist);
    h5r->Divide(fhRegionMultMaxVsEt,fhEleadingPt,1,1);
    h6r->Divide(fhRegionMultMinVsEt,fhEleadingPt,1,1);
    //h5r->Scale( areafactor );
    h5r->SetYTitle("N_{Tracks}^{90}");
    h5r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h5r->SetMarkerStyle(20);
    h5r->DrawCopy("p");
    h6r->SetMarkerStyle(21);
    h6r->SetMarkerColor(2);
    h6r->SetYTitle("N_{Tracks}^{90}");
    h6r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h6r->DrawCopy("p same");
    c1->Update();
    
    //Get Normalization
    //Double_t xsec = fh1Xsec->GetBinContent(1);
    Double_t xsec = fh1Xsec->GetBinContent(1);
    Double_t ntrials = fh1Trials->GetBinContent(1);
    Double_t normFactor = xsec/ntrials;
    if(debug > 1)Printf("xSec %f nTrials %f Norm %f \n",xsec,ntrials,normFactor);
    
    
    //Jet pT distribution
    TCanvas* c2 = new TCanvas("c2","Jet Pt dist",160,160,1200,800);
    TH1 * copy = 0;
    c2->Divide(2,2);
    c2->cd(1);
    fhEleadingPt->SetMarkerStyle(20);
    fhEleadingPt->SetMarkerColor(2);
    //if( normFactor > 0.) fhEleadingPt->Scale(normFactor);
    //fhEleadingPt->Draw("p");
    copy = fhEleadingPt->DrawCopy("p");
    if( normFactor > 0.) copy->Scale(normFactor);
    gPad->SetLogy();
    
    c2->cd(2);
    Int_t xbin1 = fhdNdEtaPhiDist->GetYaxis()->FindFixBin(minJetPtInHist);
    Int_t xbin2 = fhdNdEtaPhiDist->GetYaxis()->FindFixBin(maxJetPtInHist);
    TH1D* dNdEtaPhiDistAllJets = fhdNdEtaPhiDist->ProjectionX("dNdEtaPhiDistAllJets",xbin1,xbin2);
    dNdEtaPhiDistAllJets->SetMarkerStyle(20);
    dNdEtaPhiDistAllJets->SetMarkerColor(2);
    dNdEtaPhiDistAllJets->DrawCopy("p");
    gPad->SetLogy();
    
    c2->cd(3);      
    fhNJets->DrawCopy();
    //c2->cd(4);      
    //fhValidRegion->DrawCopy("p");
    //fhTransRegPartPtDist  = (TH1F*)fListOfHistos->At(2);
    //fhRegionMultMin          = (TH1F*)fListOfHistos->At(3);
    //fhMinRegAvePt     = (TH1F*)fListOfHistos->At(4);
    //fhMinRegSumPt     = (TH1F*)fListOfHistos->At(5);   
    //fhMinRegMaxPtPart     = (TH1F*)fListOfHistos->At(6);
    //fhMinRegSumPtvsMult   = (TH1F*)fListOfHistos->At(7);
    c2->Update(); 
    
    //pT distributions
    TCanvas* c3 = new TCanvas("c3"," pT dist",160,160,1200,800);
    c3->Divide(2,2);
    c3->cd(1);
     //fhTransRegPartPtDist->SetMarkerStyle(20);
     //fhTransRegPartPtDist->SetMarkerColor(2); 
     //fhTransRegPartPtDist->Scale(areafactor/fhTransRegPartPtDist->GetEntries());
     //fhTransRegPartPtDist->DrawCopy("p");
     //gPad->SetLogy();
     

    c3->cd(2); 
    fhMinRegSumPt->SetMarkerStyle(20);
    fhMinRegSumPt->SetMarkerColor(2);  
    //fhMinRegSumPt->Scale(areafactor);
    fhMinRegSumPt->DrawCopy("p");
    gPad->SetLogy();
    
    c3->cd(3);
    fhMinRegAvePt->SetMarkerStyle(20);
    fhMinRegAvePt->SetMarkerColor(2);  
    //fhMinRegAvePt->Scale(areafactor);
    fhMinRegAvePt->DrawCopy("p");
    gPad->SetLogy();
    
    c3->cd(4);
    TH1F *h7r = new TH1F("hRegionMultMinVsMult",   "",  21, -0.5,   20.5);
    h7r->Divide(fhMinRegSumPtvsMult,fhRegionMultMin,1,1);
    h7r->SetMarkerStyle(20);
    h7r->SetMarkerColor(2);   
    h7r->DrawCopy("p");
    c3->Update();
    
    

    //Save canvas
    c1->SaveAs("c1.pdf");
    AliInfo("Canvas 1 saved");
    c2->SaveAs("c2.pdf");
    AliInfo("Canvas 2 saved");
    c3->SaveAs("c3.pdf");
    AliInfo("Canvas 3 saved");
  
}

//____________________________________________________________________
void AliHistogramsUE::FillHistogram(const char* name, Double_t fillX){
  
  // Fill 1D histogram with double
  ((TH1F*)fListOfHistos->FindObject(name))->Fill(fillX);

}

//____________________________________________________________________
void AliHistogramsUE::FillHistogram(const char* name, Int_t fillX){
  
  // Fill 1D histogram with integer
  ((TH1F*)fListOfHistos->FindObject(name))->Fill(fillX);

}

//____________________________________________________________________
void AliHistogramsUE::FillHistogram(const char* name, Double_t fillX, Double_t fillY){
  
 // Case of TH1F with weight or TH2F w/o weight
 TObject *obj = fListOfHistos->FindObject(name);
 if (obj->InheritsFrom("TH1F")){
 	((TH1F*)fListOfHistos->FindObject(name))->Fill(fillX, fillY);
 } else {
 	((TH2F*)fListOfHistos->FindObject(name))->Fill(fillX, fillY);
	}


}

//____________________________________________________________________
void AliHistogramsUE::FillHistogram(const char* name, Double_t fillX, Double_t fillY, Double_t weight){

  // Fill 2D histogram with double and weight   
  ((TH2F*)fListOfHistos->FindObject(name))->Fill(fillX, fillY, weight);

}

//____________________________________________________________________
void AliHistogramsUE::FillHistogram(const char* name, Double_t fillX, Int_t fillY, Double_t weight){
  
  // Fill 2D histogram with integer and weight   
  ((TH2F*)fListOfHistos->FindObject(name))->Fill(fillX, fillY, weight);

}

//____________________________________________________________________
TObjArray*  AliHistogramsUE::GetHistosForPlotting(TString data, TString branches){

        // Instance filled histos for plotting purpose
	printf("Creating histograms ... \n");
        	
	printf("Reading file: %s\n",data.Data());
	
        // Read input files ----------------------------------------- 
        TFile *fdata = new TFile(data.Data());
	TDirectoryFile *ddata[20];
	TList *ldata[20];

	TObjArray  *arrb=branches.Tokenize(";");
	TIter next(arrb);
	TObject *o=0;
	Int_t br=0;
	while ( (o=next()) ){
	   	ddata[br] = (TDirectoryFile*)fdata->Get(Form("PWG4_UE%s",o->GetName()));
		if(!ddata[br]) printf("ERROR: No histo dir found! \n");
	        ldata[br] = (TList*)ddata[br]->Get(Form("histosUE%s",o->GetName()));
		printf("Reading branch: %s\n",o->GetName());
		if(!ldata[br]) printf("ERROR: No histo list found! \n");
		br++;
	}
	
	TObjArray *arr=new TObjArray;

	TH1F *hjets[20];		// accepted leading jets
        TH1F *hnjets[20];          // number of accepted jets
	TH2F *hetaphi[20];         // delta-phi particle-jet correlation
	TH2F *hptfull[20];         // particle pT all regions vs. jet pT 
	TH2F *hpttransv[20];       // particle pT transv. regions vs. jet pT 
	TH1F *hmax[20];    	// sum pT in MAX region
        TH1F *hmin[20];    	// sum pT in MIN region
        TH1F *hmultmax[20];	// multiplicity in MAX region
        TH1F *hmultmin[20];       // multiplicity in MIN region

	for (Int_t i =0; i<br; i++){

	        //Number of jets --------------------------------------------
		hnjets[i] = (TH1F*) ldata[i]->FindObject("fhNJets"); 
       		hnjets[i]->GetXaxis()->SetTitle("Number of jets");
        	hnjets[i]->GetYaxis()->SetTitle("1/n_{ev} dN/dN_{jets}");
        	hnjets[i]->SetMarkerStyle(20); 
        	hnjets[i]->SetMarkerColor(1+i);


	        //Leading jets ----------------------------------------------
        	hjets[i] = (TH1F*) ldata[i]->FindObject("hEleadingPt"); 
        	hjets[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
        	hjets[i]->GetYaxis()->SetTitle("1/n_{ev} dN/dp_{T} (|#eta<0.5|)");
        	hjets[i]->SetMarkerStyle(20); 
        	hjets[i]->SetMarkerColor(1+i);
		hjets[i]->SetMinimum(0.1);
		hjets[i]->SetMaximum(1000.);


	        //Transverse Region MAX -------------------------------------
		hmax[i] = (TH1F*) ldata[i]->FindObject("hRegionSumPtMaxVsEt");
		if (!hmax[i])AliInfo("Histo not found!!!");
		hmax[i]->GetXaxis()->SetTitle("Leading jet P_{T} (GeV/c)");
		hmax[i]->GetYaxis()->SetTitle("P_{T}^{90,max} (GeV/c)");
		hmax[i]->SetMarkerStyle(20);
		hmax[i]->SetMarkerColor(1+i);
        	hmax[i]->Divide(hjets[i]); // normalize for jet spectrum
		hmax[i]->SetMaximum(5.);


	        //Transverse Region MIN -------------------------------------
		hmin[i] = (TH1F*) ldata[i]->FindObject("hRegionSumPtMinVsEt");
		hmin[i]->GetXaxis()->SetTitle("Leading jet P_{T} (GeV/c)");
		hmin[i]->GetYaxis()->SetTitle("P_{T}^{90,min} (GeV/c)");
		hmin[i]->SetMarkerStyle(20);
		hmin[i]->SetMarkerColor(1+i);
        	hmin[i]->SetMaximum(3.);
        	hmin[i]->Divide(hjets[i]); // normalize for jet spectrum


	        //Multiplicity MAX ------------------------------------------
		hmultmax[i] = (TH1F*) ldata[i]->FindObject("hRegionMultMaxVsEt");
		hmultmax[i]->GetXaxis()->SetTitle("Leading Jet P_{T} (GeV/c)");
		hmultmax[i]->GetYaxis()->SetTitle("N_{ch}^{90,max}");
		hmultmax[i]->SetMarkerStyle(20);
		hmultmax[i]->SetMarkerColor(1+i);
		hmultmax[i]->SetMaximum(10.);
        	hmultmax[i]->Divide(hjets[i]); // normalize for jet spectrum


	        //Multiplicity MIN ------------------------------------------
		hmultmin[i] = (TH1F*) ldata[i]->FindObject("hRegionMultMinVsEt");
		hmultmin[i]->GetXaxis()->SetTitle("Leading Jet P_{T} (GeV/c)");
		hmultmin[i]->GetYaxis()->SetTitle("N_{ch}^{90,min}");
		hmultmin[i]->SetMarkerStyle(20);
		hmultmin[i]->SetMarkerColor(1+i);
		hmultmin[i]->SetMaximum(3.);
        	hmultmin[i]->Divide(hjets[i]); // normalize for jet spectrum
            

	        // Phi-correlation with leading jet --------------------------
		hetaphi[i] = (TH2F*) ldata[i]->FindObject("hdNdEtaPhiDist");
		hetaphi[i]->GetXaxis()->SetTitle("#Delta #phi (w.r.t. leading jet)");
		hetaphi[i]->SetMarkerStyle(20);
            

	        // pT distribution in full region vs. jet pT --------------------------
		hptfull[i] = (TH2F*) ldata[i]->FindObject("hFullRegPartPtDistVsEt");
		hptfull[i]->GetYaxis()->SetTitle("Leading Jet P_{T} (GeV/c)");
		hptfull[i]->GetXaxis()->SetTitle("Track P_{T} (GeV/c)");
            

	        // pT distribution in transv region vs. jet pT --------------------------
		hpttransv[i] = (TH2F*) ldata[i]->FindObject("hTransRegPartPtDistVsEt");
		hpttransv[i]->GetYaxis()->SetTitle("Leading Jet P_{T} (GeV/c)");
		hpttransv[i]->GetXaxis()->SetTitle("Track P_{T} (GeV/c)");
            

	        // Return Histos
        	arr->Add(hnjets[i]);      //at 0 number of jets
		arr->Add(hjets[i]);	//at 1 leading jets
		arr->Add(hmax[i]);	//at 2 sum pT MAX
		arr->Add(hmin[i]);	//at 3 sumpT MIN
		arr->Add(hmultmax[i]);   //at 4 multiplicity MAX 
		arr->Add(hmultmin[i]);   //at 5 multiplicity MIN
		arr->Add(hetaphi[i]);     //at 6 phi correlation
        	arr->Add(hptfull[i]);     //at 7 pT distr in full region
        	arr->Add(hpttransv[i]);   //at 8 pT distr in transv region

	}
	return arr;
}

//_______________________________________________________________________________________
void AliHistogramsUE::SetStyle(){

  // Set plotting style
  gPad->SetFrameFillColor(0);
  gPad->SetFillColor(0);
  gPad->SetBorderSize(2);
  gPad->SetGridy();
  gPad->SetFrameBorderMode(0);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

}


//____________________________________________________________________
TList* AliHistogramsUE::GetListOfHistos(){

  // Return list of relevant histograms 
  return fListOfHistos;

}

//____________________________________________________________________
void AliHistogramsUE::PlotBranchesUE(TString file, TString branches, Double_t minJetProjection){

  // Function to be called by external macro to plot analysis from different jet branches 
  
  //Count the branches 
  TObjArray  *arr=branches.Tokenize(";");
  TIter next(arr);
  TObject *o=0;
  Int_t br=0;
  while ( (o=next()) ){
	br++;
        }
  

  //Canvas
  TObjArray *arrC=CreateCanvas(9);
 
  //Histograms 
  TObjArray *arrH=GetHistosForPlotting(file.Data(),branches.Data());
  TH1F *hnjets[20];
  TH1F *hjets[20];
  TH1F *hmax[20];
  TH1F *hmin[20];
  TH1F *hmultmax[20];
  TH1F *hmultmin[20];
  TH2F *hetaphi[20];
  TH2F *hptfull[20];
  TH2F *hpttransv[20];

  for (Int_t i= 0; i<br; i++){
	hnjets[i]=((TH1F*)arrH->At(9*i));
	hjets[i]=((TH1F*)arrH->At(1+(9*i)));
	hmax[i]=((TH1F*)arrH->At(2+(9*i)));
	hmin[i]=((TH1F*)arrH->At(3+(9*i)));
	hmultmax[i]=((TH1F*)arrH->At(4+(9*i)));
	hmultmin[i]=((TH1F*)arrH->At(5+(9*i)));
	hetaphi[i]=((TH2F*)arrH->At(6+(9*i))); 
	hptfull[i]=((TH2F*)arrH->At(7+(9*i))); 
	hpttransv[i]=((TH2F*)arrH->At(8+(9*i))); 
	}

  //Define jet-pT range in projections
  Int_t binstartprojection = hetaphi[0]->GetYaxis()->FindFixBin(minJetProjection); //be careful...  
  Int_t binstopprojection = hetaphi[0]->GetYaxis()->GetNbins();

  Double_t normphi[20]; 

  for (Int_t i= 0; i<br; i++){
	normphi[i]=hjets[i]->Integral(binstartprojection,binstopprojection);	
	hnjets[i]->Scale(1./(hnjets[i]->GetBinWidth(1)*hnjets[i]->GetEntries()));
        hjets[i]->Scale(1./(hjets[i]->GetBinWidth(1)*hnjets[i]->GetEntries()));
        }

  //LEGENDS
  //Legend jets
  TLegend *leg=new TLegend(0.5,0.6,0.89,0.89);  // for jet spectrum
  leg->SetFillColor(0);
  leg->SetHeader("Jet Finders:");
  //Legend density
  TLegend *legd = new TLegend(0.1077586,0.6016949,0.4971264,0.8919492,NULL,"brNDC");
  legd->SetFillColor(0);
  legd->SetHeader("Jet Finders:");
  //Legend pT distributions
  TLegend *legp = new TLegend(0.1364943,0.1292373,0.5258621,0.4194915,NULL,"brNDC");
  legp->SetFillColor(0);
  legp->SetHeader("Jet Finders:");

  arr=branches.Tokenize(";");
  TIter next1(arr);
  o=0;
  Int_t brleg=0;
  while ( (o=next1()) ){
	leg->AddEntry(hjets[brleg],Form("UE%s",o->GetName()),"p");
	legd->AddEntry(hjets[brleg],Form("UE%s",o->GetName()),"p");
	legp->AddEntry(hjets[brleg],Form("UE%s",o->GetName()),"p");
        brleg++;
        }

  //1) NUMBER OF CLUSTERS 
  TCanvas *c0=((TCanvas*)arrC->At(0)); 
  c0->cd();
  SetStyle();
  gPad->SetLogy();
  hnjets[0]->SetMaximum(1.4);
  hnjets[0]->Draw("E1");
  for (Int_t i= 1; i<br; i++){
	hnjets[i]->Draw("same");
        }
  leg->Draw("same");

  //2) LEADING CLUSTERS pT 
  TCanvas *c1=((TCanvas*)arrC->At(1)); 
  c1->cd();
  SetStyle();
  gPad->SetLogy();
  hjets[0]->Draw("E1");
  hjets[0]->SetMaximum(1.4);
  for (Int_t i= 1; i<br; i++){
	hjets[i]->Draw("same");
        }
  leg->Draw("same");

  //3) SUM-pT IN MAX REGION
  TCanvas *c2=((TCanvas*)arrC->At(2)); 
  c2->cd();
  SetStyle();
  hmax[0]->Draw("E1");
  for (Int_t i= 1; i<br; i++){
	hmax[i]->Draw("same");
        }
  leg->Draw("same");

  //4) SUM-pT IN MIN REGION
  TCanvas *c3=((TCanvas*)arrC->At(3)); 
  c3->cd();
  SetStyle();
  hmin[0]->GetYaxis()->SetRangeUser(0.,2.5);
  hmin[0]->Draw("E1");
  for (Int_t i= 1; i<br; i++){
	hmin[i]->Draw("same");
        }
  leg->Draw("same");

  //5) MULTIPLICITY IN MAX REGION
  TCanvas *c4=((TCanvas*)arrC->At(4));
  c4->cd();
  SetStyle();
  hmultmax[0]->GetYaxis()->SetRangeUser(0.,5.);
  hmultmax[0]->Draw("E1");
  for (Int_t i= 1; i<br; i++){
	hmultmax[i]->Draw("same");
        }
  leg->Draw("same");

  //6) MULTIPLICITY IN MIN REGION
  TCanvas *c5=((TCanvas*)arrC->At(5));
  c5->cd();
  SetStyle();
  hmultmin[0]->GetYaxis()->SetRangeUser(0.,2.5);
  hmultmin[0]->Draw("E1");
  for (Int_t i= 1; i<br; i++){
	hmultmin[i]->Draw("same");
        }
  leg->Draw("same");

  //7) JET-TRACK CORRELATION
  TCanvas *c6=((TCanvas*)arrC->At(6));
  c6->cd();
  SetStyle();
  gPad->SetLogy();
  TH1D *tmpetaphi[20];
  for (Int_t i= 0; i<br; i++){
	tmpetaphi[i]=new TH1D();
	tmpetaphi[i]=hetaphi[i]->ProjectionX(Form("data%d",i),binstartprojection);
	tmpetaphi[i]->GetYaxis()->SetTitle("1/n_{lj} dN/d#Delta #phi");
	tmpetaphi[i]->Scale(1./(hetaphi[i]->GetBinWidth(1)*normphi[i]));
	tmpetaphi[i]->SetMarkerColor(i+1);
	tmpetaphi[i]->GetXaxis()->SetLimits(-3.*TMath::Pi()/2.,TMath::Pi()/2.);
	tmpetaphi[i]->GetYaxis()->SetRangeUser(0.5,20.);
	if (i==0) tmpetaphi[i]->Draw("E1");
	else tmpetaphi[i]->Draw("same");
	// evaluate mean multiplicity in transverse regions
	Double_t err1=0.;
	Double_t err2=0.;
	Double_t err3=0.;
	Double_t mean=(tmpetaphi[i]->IntegralAndError(1,5,err1)+tmpetaphi[i]->IntegralAndError(27,36,err2)+tmpetaphi[i]->IntegralAndError(58,62,err3))/(20.);
	err1=TMath::Sqrt(err1*err1+err2*err2+err3*err3)/20.;
	Printf("Branch: %d  MeanTransvMult: %f err: %f",i,mean,err1);
        }
  legd->Draw("same");

  //8) TRACK-pT DISTRIBUTION IN FULL REGION
  TCanvas *c7=((TCanvas*)arrC->At(7));
  c7->cd();
  SetStyle();
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTitle("Full region (projection X)");
  TH1D *tmpfull[20];
  for (Int_t i= 0; i<br; i++){
	tmpfull[i]=new TH1D();
	tmpfull[i]=hptfull[i]->ProjectionX(Form("full%d",i),binstartprojection);
	tmpfull[i]->GetYaxis()->SetTitle("1/n_{lj} dN/dp_{T}");
        tmpfull[i]->Scale(1./(tmpfull[i]->GetBinWidth(1)*normphi[i]));
	tmpfull[i]->SetMarkerStyle(20);
	tmpfull[i]->SetMarkerColor(i+1);
	tmpfull[i]->SetMaximum(100.);
	if (i==0)tmpfull[i]->Draw("E1");
	else tmpfull[i]->Draw("same"); 
        }
  legp->Draw("same");

  //9) TRACK-pT DISTRIBUTION IN TRANSVERSE (MIN+MAX) REGION
  TCanvas *c8=((TCanvas*)arrC->At(8));
  c8->cd();
  SetStyle();
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTitle("Transverse regions (projection X)");
  TH1D *tmptransv[20];
  for (Int_t i= 0; i<br; i++){
	tmptransv[i]=new TH1D();
	tmptransv[i]=hpttransv[i]->ProjectionX(Form("transv%d",i),binstartprojection);
	tmptransv[i]->GetYaxis()->SetTitle("1/n_{lj} dN/dp_{T}");
        tmptransv[i]->Scale(1./(tmptransv[i]->GetBinWidth(1)*normphi[i]));
	tmptransv[i]->SetMarkerStyle(20);
	tmptransv[i]->SetMarkerColor(i+1);
	tmptransv[i]->SetMaximum(10.);
	if (i==0)tmptransv[i]->Draw("E1");
	else tmptransv[i]->Draw("same"); 
        }
  legp->Draw("same");
}
