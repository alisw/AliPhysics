/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

 /************************************ 
 * flow analysis with multi-particle *
 *           correlations            * 
 *                                   * 
 * author: Ante Bilandzic            * 
 *        (abilandzic@gmail.com)     *
 ************************************/ 

#define AliFlowAnalysisWithMultiparticleCorrelations_cxx

#include "AliFlowAnalysisWithMultiparticleCorrelations.h"

using std::endl;
using std::cout;
using std::flush;

//================================================================================================================

ClassImp(AliFlowAnalysisWithMultiparticleCorrelations)

AliFlowAnalysisWithMultiparticleCorrelations::AliFlowAnalysisWithMultiparticleCorrelations(): 
 // 0.) Base:
 fHistList(NULL),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL)
 {
  // Constructor.  
  
  // a) Book grandmother of all lists;
  // b) Initialize all arrays.

  // a) Book grandmother of all lists:
  fHistList = new TList();
  fHistList->SetName("cobjMPC");
  fHistList->SetOwner(kTRUE);

  // b) Initialize all arrays:
  this->InitializeArraysForControlHistograms();
    
 } // end of AliFlowAnalysisWithMultiparticleCorrelations::AliFlowAnalysisWithMultiparticleCorrelations()
 
//================================================================================================================  

AliFlowAnalysisWithMultiparticleCorrelations::~AliFlowAnalysisWithMultiparticleCorrelations()
{
 // Destructor.
 
 delete fHistList;

} // end of AliFlowAnalysisWithMultiparticleCorrelations::~AliFlowAnalysisWithMultiparticleCorrelations()

//================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::Init()
{
 // Well, this is method Init().

 // a) Trick to avoid name clashes, part 1; 
 // b) Cross-check the initial settings before starting this adventure;
 // c) Book all objects;
 // d) Set all flags;
 // *) Trick to avoid name clashes, part 2. 

 // a) Trick to avoid name clashes, part 1: 
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);
 
 // b) Cross-check the initial settings before starting this adventure:
 this->CrossCheckSettings();

 // c) Book all objects:
 this->BookAndNestAllLists();
 this->BookEverythingForControlHistograms();

 // d) Set all flags:
 // ... 

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::Init()

//================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::Make(AliFlowEventSimple *anEvent)
{
 // Running over data only in this method.
 
 // a) Fill control histograms;

 // a) Fill control histograms:
 this->FillControlHistograms(anEvent);
 
 // ...
 
} // end of AliFlowAnalysisWithMultiparticleCorrelations::Make(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::Finish()
{
 // Draw the curtains. 
 
 // ...

} // end of AliFlowAnalysisWithMultiparticleCorrelations::Finish()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()
{
 // Cross-check all initial settings in this method. 
 
 // ...

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;

 // ...

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("Control Histograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookAndNestAllLists()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::WriteHistograms(TString outputFileName)
{
 // Store the final results in output file <outputFileName>.root.

 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete output;

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::WriteHistograms(TString outputFileName)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::WriteHistograms(TDirectoryFile *outputFileName)
{
 // Store the final results in output file <outputFileName>.root.

 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(),TObject::kSingleKey);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::WriteHistograms(TDirectoryFile *outputFileName)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForControlHistograms()
{
 // Book all the stuff for control histograms.

 // a) Book the profile holding all the flags for control histograms;
 // b) Book all control histograms;
 //  b0) Book TH1D *fKinematicsHist[2][3]; 
 //  b1) Book TH1D *fMultDistributionsHist[3];
 //  b2) Book TH2D *fMultCorrelationsHist[3].

 // a) Book the profile holding all the flags for control histograms: TBI stil incomplete 
 fControlHistogramsFlagsPro = new TProfile("fControlHistogramsFlagsPro","Flags for control histograms",4,0,4);
 fControlHistogramsFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsFlagsPro->SetMarkerStyle(25);
 fControlHistogramsFlagsPro->SetLabelSize(0.04);
 fControlHistogramsFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsFlagsPro->SetStats(kFALSE);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(1,"Use default values");
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(2,"nBins (p_{T})");
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(3,"p_{T,min}");
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(4,"p_{T,max}");
 // ...
 fControlHistogramsList->Add(fControlHistogramsFlagsPro);

 // b) Book all control histograms: // TBI add setters for all these values
 //  b0) Book TH1D *fKinematicsHist[2][3]:
 Int_t nBins[2][3] = {{360,1000,1000},{360,1000,1000}}; // [RP,POI][phi,pt,eta]
 Double_t min[2][3] = {{0.,0.,-1.},{0.,0.,-1.}}; // [RP,POI][phi,pt,eta]
 Double_t max[2][3] = {{TMath::TwoPi(),10.,1.},{TMath::TwoPi(),10.,1.}}; // [RP,POI][phi,pt,eta]
 TString name[2][3] = {{"RP,phi","RP,pt","RP,eta"},{"POI,phi","POI,pt","POI,eta"}}; // [RP,POI][phi,pt,eta]
 TString title[2] = {"Reference particles (RP)","Particles of interest (POI)"}; // [RP,POI]
 Int_t lineColor[2] = {kBlue,kRed}; // [RP,POI]
 Int_t fillColor[2] = {kBlue-10,kRed-10}; // [RP,POI]
 TString xAxisTitle[3] = {"#phi","p_{T}","#eta"}; // [phi,pt,eta]
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   fKinematicsHist[rp][ppe] = new TH1D(name[rp][ppe].Data(),title[rp].Data(),nBins[rp][ppe],min[rp][ppe],max[rp][ppe]);
   fKinematicsHist[rp][ppe]->GetXaxis()->SetTitle(xAxisTitle[ppe].Data());
   fKinematicsHist[rp][ppe]->SetLineColor(lineColor[rp]);
   fKinematicsHist[rp][ppe]->SetFillColor(fillColor[rp]);
   fControlHistogramsList->Add(fKinematicsHist[rp][ppe]);
  }
 }

 //  b1) Book TH1D *fMultDistributionsHist[3]: // TBI add setters for all these values
 Int_t nBinsMult[3] = {3000,3000,3000}; // [RP,POI,reference multiplicity]
 Double_t minMult[3] = {0.,0.,0.}; // [RP,POI,reference multiplicity]
 Double_t maxMult[3] = {3000.,3000.,3000.}; // [RP,POI,reference multiplicity]
 TString nameMult[3] = {"Multiplicity (RP)","Multiplicity (POI)","Multiplicity (REF)"}; // [RP,POI,reference multiplicity]
 TString titleMult[3] = {"Reference particles (RP)","Particles of interest (POI)",""}; // [RP,POI,reference multiplicity]
 Int_t lineColorMult[3] = {kBlue,kRed,kGreen+2}; // [RP,POI,reference multiplicity]
 Int_t fillColorMult[3] = {kBlue-10,kRed-10,kGreen-10}; // [RP,POI,reference multiplicity]
 TString xAxisTitleMult[3] = {"Multiplicity (RP)","Multiplicity (POI)","Multiplicity (REF)"}; // [phi,pt,eta]
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  fMultDistributionsHist[rprm] = new TH1D(nameMult[rprm].Data(),titleMult[rprm].Data(),nBinsMult[rprm],minMult[rprm],maxMult[rprm]);
  fMultDistributionsHist[rprm]->GetXaxis()->SetTitle(xAxisTitleMult[rprm].Data());
  fMultDistributionsHist[rprm]->SetLineColor(lineColorMult[rprm]);
  fMultDistributionsHist[rprm]->SetFillColor(fillColorMult[rprm]);
  fControlHistogramsList->Add(fMultDistributionsHist[rprm]);
 } // for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]

 //  b2) Book TH2D *fMultCorrelationsHist[3]: TBI too large objects to store in this way
 // ...
 fMultCorrelationsHist[0] = new TH2D("Multiplicity (RP vs. POI)","Multiplicity (RP vs. POI)",nBinsMult[0],minMult[0],maxMult[0],nBinsMult[1],minMult[1],maxMult[1]);
 fMultCorrelationsHist[0]->GetXaxis()->SetTitle(xAxisTitleMult[0].Data());
 fMultCorrelationsHist[0]->GetYaxis()->SetTitle(xAxisTitleMult[1].Data());
 fControlHistogramsList->Add(fMultCorrelationsHist[0]);

 // ...
 fMultCorrelationsHist[1] = new TH2D("Multiplicity (RP vs. REF)","Multiplicity (RP vs. REF)",nBinsMult[0],minMult[0],maxMult[0],nBinsMult[2],minMult[2],maxMult[2]);
 fMultCorrelationsHist[1]->GetXaxis()->SetTitle(xAxisTitleMult[0].Data());
 fMultCorrelationsHist[1]->GetYaxis()->SetTitle(xAxisTitleMult[2].Data());
 fControlHistogramsList->Add(fMultCorrelationsHist[1]);

 // ...
 fMultCorrelationsHist[2] = new TH2D("Multiplicity (POI vs. REF)","Multiplicity (POI vs. REF)",nBinsMult[1],minMult[1],maxMult[1],nBinsMult[2],minMult[2],maxMult[2]);
 fMultCorrelationsHist[2]->GetXaxis()->SetTitle(xAxisTitleMult[1].Data());
 fMultCorrelationsHist[2]->GetYaxis()->SetTitle(xAxisTitleMult[2].Data());
 fControlHistogramsList->Add(fMultCorrelationsHist[2]);

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::FillControlHistograms(AliFlowEventSimple *anEvent)
{
 // Fill control histograms. 

 // a) Fill TH1D *fKinematicsHist[2][3];
 // b) Fill TH1D *fMultDistributionsHist[3]; 
 // c) Fill TH2D *fMultCorrelationsHist[3].  

 // a) Fill TH1D *fKinematicsHist[2][3]:
 Int_t nTracks = anEvent->NumberOfTracks();
 for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 {
  AliFlowTrackSimple *pTrack = anEvent->GetTrack(t);
  if(!pTrack){printf("\n AAAARGH: pTrack is NULL in MPC::FCH() !!!!");continue;}
  if(pTrack)
  {
   Double_t dPhi = pTrack->Phi(); 
   //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
   //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
   Double_t dPt = pTrack->Pt();
   Double_t dEta = pTrack->Eta();
   Double_t dPhiPtEta[3] = {dPhi,dPt,dEta};
   for(Int_t rp=0;rp<2;rp++) // [RP,POI]
   {
    for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
    {
     if((0==rp && pTrack->InRPSelection()) || (1==rp && pTrack->InPOISelection())) // TBI 
     { 
      fKinematicsHist[rp][ppe]->Fill(dPhiPtEta[ppe]);
     }
    } // for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
   } // for(Int_t rp=0;rp<2;rp++) // [RP,POI]
  } // if(pTrack)  
 } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks

 // b) Fill TH1D *fMultDistributionsHist[3]: 
 Double_t dMultRP = anEvent->GetNumberOfRPs(); // TBI shall I promote these 3 variables into data members? 
 Double_t dMultPOI = anEvent->GetNumberOfPOIs();
 Double_t dMultREF = anEvent->GetReferenceMultiplicity();
 Double_t dMult[3] = {dMultRP,dMultPOI,dMultREF};
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  fMultDistributionsHist[rprm]->Fill(dMult[rprm]);      
 } 

 // c) Fill TH2D *fMultCorrelationsHist[3]:  
 fMultCorrelationsHist[0]->Fill((Int_t)dMultRP,(Int_t)dMultPOI); // RP vs. POI
 fMultCorrelationsHist[1]->Fill((Int_t)dMultRP,(Int_t)dMultREF); // RP vs. refMult
 fMultCorrelationsHist[2]->Fill((Int_t)dMultPOI,(Int_t)dMultREF); // POI vs. refMult

} // void AliFlowAnalysisWithMultiparticleCorrelations::FillControlHistograms(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForControlHistograms()
{
 // Initialize all arrays for control histograms.

 // a) Initialize TH1D *fKinematicsHist[2][3];
 // b) Initialize TH1D *fMultDistributionsHist[3]; 
 // c) Initialize TH2D *fMultCorrelationsHist[3].  

 // a) Initialize TH1D *fKinematicsHist[2][3]:
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   fKinematicsHist[rp][ppe] = NULL;
  } 
 } 

 // b) Initialize TH1D *fMultDistributionsHist[3]:
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  fMultDistributionsHist[rprm] = NULL;      
 } 

 // c) Initialize TH2D *fMultCorrelationsHist[3]: 
 for(Int_t r=0;r<3;r++) // [RP vs. POI, RP vs. refMult, POI vs. refMult]  
 {
  fMultCorrelationsHist[r] = NULL; 
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetOutputHistograms(TList *histList)
{
 // Get pointers for everything and everywhere from the base list. 

 if(histList)
 {	
  this->SetHistList(histList);
  if(!fHistList)
  {
   printf("\n AAAARGH (MPC): fHistList is NULL in AFAWMPC::GOH() !!!!\n\n"); exit(0);
  }
 } else 
   {
    printf("\n AAAARGH (MPC): histList is NULL in AFAWMPC::GOH() !!!!\n\n"); exit(0);
   }
   
} // void AliFlowAnalysisWithMultiparticleCorrelations::GetOutputHistograms(TList *histList)

//=======================================================================================================================



