
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

//------------------------------------------------------------------------------
// Implementation of AliPerformancePtCalib class. It fills histograms with ESD or
// TPC track information to study charge/pt spectra.
// To analyse the output of AliPerformancePtCalib use AliPerfAnalyzeInvPt class:
// Projection of charge/pt vs theta and vs phi resp. histoprams will be fitted
// with either polynomial or gaussian fit function to extract minimum position of
// charge/pt.
// Fit options and theta, phi bins can be set by user.
// Attention: use the Set functions of AliPerformancePtCalib  when running
// AliPerformancePtCalib::Analyse()
// The result of the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliPerformancePtCalib.
//
// Author: S.Schuchmann 11/13/2009 
//         sschuchm@ikf.uni-frankfurt.de
// Updated: J. Salzwedel 01/12/2014
//------------------------------------------------------------------------------

/*
 
// after running the performance task, read the file, and get component

TFile f("Output.root");
AliPerformancePtCalib * compObj = (AliPerformancePtCalib*)coutput->FindObject("AliPerformancePtCalib");
 
// set phi and theta bins for fitting and analyse comparison data
compObj->SetProjBinsTheta(thetaBins,nThetaBins,minPhi, maxPhi);
compObj->SetProjBinsPhi(phiBins,nPhiBins,minTheta,maxTheta);
compObj->SetMakeFitOption(kFALSE,exclRange,fitRange);
compObj->SetDoRebin(rebin);
compObj->Analyse();
//details see functions of this class

// the output histograms/graphs will be stored in the folder "folderRes" 
compObj->GetAnalysisFolder()->ls("*");

// user can save whole comparison object (or only folder with anlysed histograms) 
// in the seperate output file (e.g.)
TFile fout("Analysed_InvPt.root","recreate");
compObj->Write(); // compObj->GetAnalysisFolder()->Write();
fout.Close();

*/
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TFolder.h"

#include "AliVEvent.h" 
#include "AliVTrack.h"
#include "AliExternalTrackParam.h"
#include "AliVfriendEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"

#include "AliPerfAnalyzeInvPt.h"
#include "AliPerformancePtCalib.h"

using namespace std;

ClassImp(AliPerformancePtCalib)

//________________________________________________________________________
AliPerformancePtCalib::AliPerformancePtCalib(TRootIOCtor* b):
   AliPerformanceObject(b),
   
   // option parameter for AliPerformancePtCalib::Analyse()
   fNThetaBins(0), 
   fNPhiBins(0),
   fMaxPhi(0),
   fMinPhi(0),
   fMaxTheta(0),
   fMinTheta(0),
   fRange(0),
   fExclRange(0),
   fFitGaus(0) ,
   fDoRebin(0),
   fRebin(0),
   // option parameter for user defined charge/pt shift
   fShift(0),
   fDeltaInvP(0),
   //options for cuts
   fOptTPC(0),
   fESDcuts(0),
   fPions(0),
   fEtaAcceptance(0),
   fList(0),
   // histograms
   fHistInvPtPtThetaPhi(0),
   fHistPtShift0(0),
   fHistPrimaryVertexPosX(0),
   fHistPrimaryVertexPosY(0),
   fHistPrimaryVertexPosZ(0),
   fHistTrackMultiplicity(0),
   fHistTrackMultiplicityCuts(0),

   fHistTPCMomentaPosP(0),
   fHistTPCMomentaNegP(0),
   fHistTPCMomentaPosPt(0),
   fHistTPCMomentaNegPt(0),
   fHistUserPtShift(0),	
   fHistdedxPions(0),
   //esd track cuts													 
   fESDTrackCuts(0),
   //pid
   fESDpid(0),
   // analysis folder 
   fAnalysisFolder(0)
{
   for(Int_t i=0; i<100; i++) {
     fThetaBins[i] = 0.0;
     fPhiBins[i] = 0.0;
   }
}

//________________________________________________________________________
AliPerformancePtCalib::AliPerformancePtCalib(const Char_t* name, const Char_t* title):
   AliPerformanceObject(name,title),
   
   // option parameter for AliPerformancePtCalib::Analyse()
   fNThetaBins(0), 
   fNPhiBins(0),
   fMaxPhi(0),
   fMinPhi(0),
   fMaxTheta(0),
   fMinTheta(0),
   fRange(0),
   fExclRange(0),
   fFitGaus(0) ,
   fDoRebin(0),
   fRebin(0),
   // option parameter for user defined charge/pt shift
   fShift(0),
   fDeltaInvP(0),
   //options for cuts
   fOptTPC(0),
   fESDcuts(0),
   fPions(0),
   fEtaAcceptance(0),
   fList(0),
   // histograms
   fHistInvPtPtThetaPhi(0),
   fHistPtShift0(0),
   fHistPrimaryVertexPosX(0),
   fHistPrimaryVertexPosY(0),
   fHistPrimaryVertexPosZ(0),
   fHistTrackMultiplicity(0),
   fHistTrackMultiplicityCuts(0),

   fHistTPCMomentaPosP(0),
   fHistTPCMomentaNegP(0),
   fHistTPCMomentaPosPt(0),
   fHistTPCMomentaNegPt(0),
   fHistUserPtShift(0),	
   fHistdedxPions(0),
   //esd track cuts													 
   fESDTrackCuts(0),
   //pid
   fESDpid(0),
   // analysis folder 
   fAnalysisFolder(0)
   
  
{
   // Constructor
   fShift = kFALSE;
   fDeltaInvP = 0.00;
   //options for cuts
   fOptTPC =  kTRUE;                      // read TPC tracks yes/no
   fESDcuts = kFALSE;
   fPions = kFALSE;
   
   //esd track cut options
   fEtaAcceptance = 0.8;
   
   // options for function AliPerformancePtCalibMC::Analyse()
   fFitGaus = kFALSE;// use gaussian function for fitting charge/pt yes/no
   fNThetaBins = 0; //number of theta bins
   fNPhiBins = 0; //number of phi bins
   fMaxPhi = 6.5;// max phi for 2D projection on theta and charge/pt axis
   fMinPhi = 0.0;// min phi for 2D projection on theta and charge/pt axis
   fMaxTheta = 3.0;// max theta for 2D projection on phi and charge/pt axis
   fMinTheta = 0.0;// min theta for 2D projection on phi and charge/pt axis
   fRange = 0; //fit range around 0
   fExclRange =0; //range of rejection of points around 0
   fDoRebin = kFALSE;
   fRebin = 0;
  
   for(Int_t i=0; i<100; i++) {
     fThetaBins[i] = 0.0;
     fPhiBins[i] = 0.0;
   }


   Init();
}
//________________________________________________________________________
AliPerformancePtCalib::~AliPerformancePtCalib() { 
   //
   // destructor
   //

   if(fList) delete fList;
   // histograms
   if( fHistInvPtPtThetaPhi)  delete fHistInvPtPtThetaPhi;fHistInvPtPtThetaPhi=0;
   if(fHistPtShift0)          delete fHistPtShift0;fHistPtShift0=0; 
   if(fHistPrimaryVertexPosX) delete fHistPrimaryVertexPosX;fHistPrimaryVertexPosX=0; 
   if(fHistPrimaryVertexPosY) delete fHistPrimaryVertexPosY;fHistPrimaryVertexPosY=0; 
   if(fHistPrimaryVertexPosZ) delete fHistPrimaryVertexPosZ;fHistPrimaryVertexPosZ=0; 
   if(fHistTrackMultiplicity) delete fHistTrackMultiplicity;fHistTrackMultiplicity=0; 
   if(fHistTrackMultiplicityCuts) delete fHistTrackMultiplicityCuts;fHistTrackMultiplicityCuts=0; 

   if(fHistTPCMomentaPosP)  delete fHistTPCMomentaPosP;fHistTPCMomentaPosP=0; 
   if(fHistTPCMomentaNegP)  delete fHistTPCMomentaNegP;fHistTPCMomentaNegP=0; 
   if(fHistTPCMomentaPosPt) delete fHistTPCMomentaPosPt;fHistTPCMomentaPosPt=0; 
   if(fHistTPCMomentaNegPt) delete fHistTPCMomentaNegPt ;fHistTPCMomentaNegPt=0; 
   if(fHistUserPtShift)     delete fHistUserPtShift;fHistUserPtShift=0; 
   if(fHistdedxPions)       delete fHistdedxPions;fHistdedxPions=0;
   //esd track cuts
   if(fESDTrackCuts) delete fESDTrackCuts;fESDTrackCuts=0;
   //pid
   if(fESDpid) delete fESDpid;fESDpid=0;
   //analysis folder 
   if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0; 
}

//________________________________________________________________________
void AliPerformancePtCalib::Init() 
{
   // Create histograms
   // Called once

   fList = new TList();
   // init folder
   fAnalysisFolder = CreateFolder("folderPt_TPC","Analysis Pt Resolution Folder");
   
   // Primary Vertex:
   fHistPrimaryVertexPosX       = new TH1F("fHistPrimaryVertexPosX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
   fList->Add(fHistPrimaryVertexPosX);
   fHistPrimaryVertexPosY       = new TH1F("fHistPrimaryVertexPosY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
   fList->Add(fHistPrimaryVertexPosY);
   fHistPrimaryVertexPosZ       = new TH1F("fHistPrimaryVertexPosZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-2.0,2.0);
   fList->Add(fHistPrimaryVertexPosZ);
   // Multiplicity:
   fHistTrackMultiplicity     = new TH1F("fHistTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 250, 0, 250);
   fList->Add(fHistTrackMultiplicity);
   fHistTrackMultiplicityCuts = new TH1F("fHistTrackMultiplicityCuts", "Multiplicity distribution;Number of tracks after cuts;Events", 250, 0, 250);
   fList->Add(fHistTrackMultiplicityCuts);
 
   // momentum histos
   //pt shift 0 only needed if shift in 1/pt is applied
   fHistPtShift0 = new TH1F("fHistPtShift0","1/pt dN/pt vs. pt of ESD track  ",800,-40.0,40.0);
   fList->Add(fHistPtShift0);
 
   // THnSparse for 1/pt and pt spectra vs angles
   const   Int_t invPtDims = 4;
   fMaxPhi = 6.52;
   fMinPhi = 0.0;
   fMaxTheta = 3.0;
   fMinTheta = 0.0;
   
   Double_t xminInvPt[invPtDims] = {-4.5,-40.0,fMinTheta,fMinPhi};
   Double_t xmaxInvPt[invPtDims] = {4.5,40.0,fMaxTheta,fMaxPhi};
   Int_t  binsInvPt[invPtDims] = {450,400,150,163};

  
   fHistInvPtPtThetaPhi = new THnSparseF("fHistInvPtPtThetaPhi","1/pt vs pt vs #theta vs #phi ",invPtDims,binsInvPt,xminInvPt,xmaxInvPt);
   fList->Add(fHistInvPtPtThetaPhi);
   
   // momentum test histos
   fHistTPCMomentaPosP  =  new TH2F("fHistTPCMomentaPosP","TPC p vs global esd track p pos",300,0.0,15.0,300,0.0,15.0);
   fList->Add(fHistTPCMomentaPosP);
   fHistTPCMomentaNegP  =  new TH2F("fHistTPCMomentaNegP","TPC p vs global esd track p neg",300,0.0,15.0,300,0.0,15.0);
   fList->Add(fHistTPCMomentaNegP);
   fHistTPCMomentaPosPt =  new TH2F("fHistTPCMomentaPosPt","TPC pt vs global esd track pt pos",300,0.0,15.0,300,0.0,15.0);
   fList->Add(fHistTPCMomentaPosPt);
   fHistTPCMomentaNegPt =  new TH2F("fHistTPCMomentaNegPt","TPC pt vs global esd track pt neg",300,0.0,15.0,300,0.0,15.0);
   fList->Add(fHistTPCMomentaNegPt);

   //user pt shift check
   fHistUserPtShift = new TH1F("fHistUserPtShift","user defined shift in 1/pt",100,-0.5,1.5);
   fList->Add(fHistUserPtShift);
   //pid by dedx
   fHistdedxPions = new TH2F ("fHistdedxPions","dEdx of pions ident. via kPID vs signed Pt",300,-15.05,15.05,200,0.0,400.0);
   fList->Add(fHistdedxPions);
   //pid
   fESDpid =  new AliESDpid();

   // esd track cuts  
   fESDTrackCuts =NULL;
}

//________________________________________________________________________
void AliPerformancePtCalib::SetPtShift(const Double_t shiftVal ) {
   //set user defined shift in charge/pt

   if(shiftVal) { fShift=kTRUE; fDeltaInvP = shiftVal; } 
}

//________________________________________________________________________
void AliPerformancePtCalib::Exec(AliMCEvent* const /*mcEvent*/, AliVEvent *const vEvent, AliVfriendEvent * const /*vfriendEvent*/, const Bool_t /*bUseMC*/, const Bool_t /*bUseVfriend*/)
{
   //exec: read esd or tpc

   if(!fESDTrackCuts)  {
     Printf("no esd track cut");
     return;
   }
   
   if (!vEvent) {
      Printf("ERROR: Event not available");
      return;
   }

   if (!(vEvent->GetNumberOfTracks()))  return;

   
   //vertex info for cut
   const AliVVertex *vtx = vEvent->GetPrimaryVertex();
   if (!vtx->GetStatus()) return; 

     
   //histo fo user defined shift in charge/pt 
   if(fShift) fHistUserPtShift->Fill(fDeltaInvP);
   
   //track multiplicity
   fHistTrackMultiplicity->Fill(vEvent->GetNumberOfTracks());


   // read primary vertex info
   Double_t tPrimaryVtxPosition[3];
   vtx->GetXYZ(tPrimaryVtxPosition);
   fHistPrimaryVertexPosX->Fill(tPrimaryVtxPosition[0]);
   fHistPrimaryVertexPosY->Fill(tPrimaryVtxPosition[1]);
   fHistPrimaryVertexPosZ->Fill(tPrimaryVtxPosition[2]);


   //_fill histos for pt spectra and shift of transverse momentum
   Int_t count=0;
 
   for(Int_t j = 0;j<vEvent->GetNumberOfTracks();j++){// track loop
     AliVTrack *vTrack = (AliVTrack*)vEvent->GetTrack(j);
     if(!vTrack) continue;

     if(fESDcuts){
       if(!fESDTrackCuts->AcceptVTrack(vTrack)) continue;
     }
      
      
     // fill histos
     if(fOptTPC){ //TPC tracks

       const AliExternalTrackParam *tpcTrack;

       if(vTrack){
	 tpcTrack = vTrack->GetTPCInnerParam(); 
       }

       if(!tpcTrack) continue;
       if(fabs(tpcTrack->Eta())>= fEtaAcceptance) continue;
      
       Double_t signedPt = tpcTrack->GetSignedPt();

       // pid
       if(fPions){
	  
	 fESDpid->GetTPCResponse().SetBetheBlochParameters(0.0283086,2.63394e+01,5.04114e-11, 2.12543e+00,4.88663e+00);

	 if( TMath::Abs(fESDpid->NumberOfSigmasTPC(vTrack,AliPID::kPion)) >1) continue;
	 fHistdedxPions->Fill(signedPt, vTrack->GetTPCsignal());
       }
	 
       Double_t invPt = 0.0;
       if(signedPt) {
	 invPt = 1.0/signedPt;

	 fHistPtShift0->Fill(signedPt);
	
	 if(fShift){Printf("user shift of momentum SET to non zero value!");
	   invPt += fDeltaInvP; //shift momentum for tests
	   if(invPt) signedPt = 1.0/invPt;
	   else {
	     continue;
	   }
	 }
	 Double_t theta = tpcTrack->Theta();
	 Double_t phi = tpcTrack->Phi();
	    
	 Double_t momAng[4] = {invPt,signedPt,theta,phi};
	 fHistInvPtPtThetaPhi->Fill(momAng);

	 Double_t pTPC = tpcTrack->P();
	 Double_t pESD = vTrack->P();
	 Double_t ptESD  = vTrack->GetSignedPt();
	
	 if(vTrack->GetSign()>0){
	   //compare momenta ESD track and TPC track
	   fHistTPCMomentaPosP->Fill(fabs(pESD),fabs(pTPC));
	   fHistTPCMomentaPosPt->Fill(fabs(ptESD),fabs(signedPt));
	 }
	 else{
	   fHistTPCMomentaNegP->Fill(fabs(pESD),fabs(pTPC));
	   fHistTPCMomentaNegPt->Fill(fabs(ptESD),fabs(signedPt));
	 }
	 count++;
       }
       else {
	 continue;
       }
     }
      
     else{// Global tracks
       if(fabs(vTrack->Eta())> fEtaAcceptance) continue;
       Double_t invPt = 0.0;
       Double_t signedPt = vTrack->GetSignedPt();
       if(signedPt){
	 invPt = 1.0/signedPt; 

	 fHistPtShift0->Fill(signedPt);
	  
	 if(fShift){Printf("user shift of momentum SET to non zero value!");
	   invPt += fDeltaInvP;//shift momentum for tests
	   if(invPt) signedPt = 1.0/invPt;
	   else continue;
	 }
	 Double_t theta = vTrack->Theta();
	 Double_t phi = vTrack->Phi();
	 Double_t momAng[4] = {invPt,signedPt,theta,phi};
	 fHistInvPtPtThetaPhi->Fill(momAng);
	 count++;
       }
     }
   }
    
   fHistTrackMultiplicityCuts->Fill(count);
}    
//______________________________________________________________________________________________________________________

void AliPerformancePtCalib::Analyse()
{
   // analyse charge/pt spectra in bins of theta and phi. Bins can be set by user
   
   THnSparseF *copyTHnSparseTheta = (THnSparseF*)fHistInvPtPtThetaPhi->Clone("copyTHnSparseTheta");
   if(!copyTHnSparseTheta) return;
   copyTHnSparseTheta->GetAxis(3)->SetRangeUser(fMinPhi,fMaxPhi);
   TH2F *histInvPtTheta = (TH2F*)copyTHnSparseTheta->Projection(2,0);
      
   THnSparseF *copyTHnSparsePhi = (THnSparseF*)fHistInvPtPtThetaPhi->Clone("copyTHnSparsePhi");
   if(!copyTHnSparsePhi) return;
   copyTHnSparsePhi->GetAxis(2)->SetRangeUser(fMinTheta,fMaxTheta);
   TH2F *histInvPtPhi   = (TH2F*)copyTHnSparsePhi->Projection(3,0);
   
   AliPerfAnalyzeInvPt *ana = new  AliPerfAnalyzeInvPt("AliPerfAnalyzeInvPt","AliPerfAnalyzeInvPt");
   if(!ana) return;
  
     
   TH1::AddDirectory(kFALSE);
 
   ana->SetProjBinsTheta(fThetaBins,fNThetaBins);
   ana->SetProjBinsPhi(fPhiBins,fNPhiBins);
   ana->SetMakeFitOption(fFitGaus,fExclRange,fRange);
   if(fDoRebin) ana->SetDoRebin(fRebin);		   
   TObjArray *aFolderObj = new TObjArray;
   if(!aFolderObj) return;
   
   ana->StartAnalysis(histInvPtTheta,histInvPtPhi,aFolderObj);
  
   // export objects to analysis folder
   fAnalysisFolder = ExportToFolder(aFolderObj);

   // delete only TObjArray
   if(aFolderObj) delete aFolderObj;
   if(ana) delete ana;
   
}

//______________________________________________________________________________________________________________________
TFolder* AliPerformancePtCalib::ExportToFolder(TObjArray * array) 
{
   // recreate folder every time and export objects to new one
   //
   AliPerformancePtCalib * comp=this;
   TFolder *folder = comp->GetAnalysisFolder();

   TString name, title;
   TFolder *newFolder = 0;
   Int_t i = 0;
   Int_t size = array->GetSize();

   if(folder) { 
      // get name and title from old folder
      name = folder->GetName();  
      title = folder->GetTitle();  

      // delete old one
      delete folder;

      // create new one
      newFolder = CreateFolder(name.Data(),title.Data());
      newFolder->SetOwner();

      // add objects to folder
      while(i < size) {
	 newFolder->Add(array->At(i));
	 i++;
      }
   }

   return newFolder;
}

//______________________________________________________________________________________________________________________
Long64_t AliPerformancePtCalib::Merge(TCollection* const list) 
{
   // Merge list of objects (needed by PROOF)

   if (!list)
      return 0;

   if (list->IsEmpty())
      return 1;

   TIterator* iter = list->MakeIterator();
   TObject* obj = 0;

   // collection of generated histograms
   Int_t count=0;
   while((obj = iter->Next()) != 0) 
      {
	 AliPerformancePtCalib* entry = dynamic_cast<AliPerformancePtCalib*>(obj);
	 if (entry == 0) continue; 
	 fHistInvPtPtThetaPhi->Add(entry->fHistInvPtPtThetaPhi);
  
	 fHistPtShift0->Add(entry->fHistPtShift0);
	 fHistPrimaryVertexPosX->Add(entry->fHistPrimaryVertexPosX);
	 fHistPrimaryVertexPosY->Add(entry->fHistPrimaryVertexPosY);
	 fHistPrimaryVertexPosZ->Add(entry->fHistPrimaryVertexPosZ);
	 fHistTrackMultiplicity->Add(entry->fHistTrackMultiplicity);
	 fHistTrackMultiplicityCuts->Add(entry->fHistTrackMultiplicityCuts);
  
	 fHistTPCMomentaPosP->Add(entry->fHistTPCMomentaPosP);
	 fHistTPCMomentaNegP->Add(entry->fHistTPCMomentaNegP);
	 fHistTPCMomentaPosPt->Add(entry->fHistTPCMomentaPosPt);
	 fHistTPCMomentaNegPt->Add(entry->fHistTPCMomentaNegPt);
	 fHistdedxPions->Add(entry->fHistdedxPions);
	 count++;
      }
  
   return count;
}

//______________________________________________________________________________________________________________________
TFolder* AliPerformancePtCalib::CreateFolder(TString name,TString title) { 
   // create folder for analysed histograms
   //
   TFolder *folder = 0;
   folder = new TFolder(name.Data(),title.Data());

   return folder;
}
//______________________________________________________________________________________________________________________
void AliPerformancePtCalib::SetProjBinsPhi(const Double_t *phiBinArray,const Int_t nphBins, const  Double_t minTheta, const  Double_t maxTheta){
   // set phi bins for Analyse()
   //set phi bins as array and set number of this array which is equal to number of bins analysed
   //the last analysed bin will always be the projection from first to last bin in the array
   if(nphBins){
      fNPhiBins = nphBins;
  
      for(Int_t k = 0;k<fNPhiBins;k++){
	 fPhiBins[k] = phiBinArray[k];
      }
      Printf("AliPerformancePtCalib::SetProjBinsPhi: number of bins in phi set to %i",fNPhiBins);
   }
   else  Printf("Warning AliPerformancePtCalib::SetProjBinsPhi: number of bins in phi NOT set!!! Default values are taken.");

   if(fabs(minTheta-maxTheta)<0.001){
      Printf("AliPerformancePtCalib::SetProjBinsPhi: theta range for projection on phi and charge/pt is too small. whole range of theta selected.");
   }
   else{
      Printf("AliPerformancePtCalib::SetProjBinsPhi: theta range for projection on phi and charge/pt is selected by user: %1.3f - %1.3f rad",minTheta,maxTheta);
      fMinTheta = minTheta;
      fMaxTheta = maxTheta;
   }
}

//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalib::SetProjBinsTheta(const Double_t *thetaBinArray, const Int_t nthBins, const Double_t minPhi, const Double_t maxPhi){
   // set theta bins for Analyse()
   //set theta bins as array and set number of this array which is equal to number of bins analysed
   //the last analysed bin will always be the projection from first to last bin in the array
   if(nthBins){
      fNThetaBins = nthBins;
      for(Int_t k = 0;k<fNThetaBins;k++){
	 fThetaBins[k] = thetaBinArray[k];
      }
      Printf("AliPerformancePtCalib::SetProjBinsTheta: number of bins in theta set to %i",fNThetaBins);
   }
   else  Printf("Warning AliPerformancePtCalib::SetProjBinsTheta: number of bins in theta NOT set!!! Default values are taken.");
   
   if(fabs(minPhi-maxPhi)<0.001){
      Printf("AliPerformancePtCalib::SetProjBinsTheta: phi range for projection on theta and charge/pt is too small. whole range of phi selected.");
   }
   else{
      Printf("AliPerformancePtCalib::SetProjBinsTheta: phi range for projection on theta and charge/pt is selected by user: %1.3f - %1.3f rad",minPhi,maxPhi);
      fMinPhi = minPhi;
      fMaxPhi = maxPhi;
   }
}
//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalib::SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR ){
   //set the fit options:
   //for usage of gaussian function instead of polynomial (default) set setGausFit=kTRUE
   //set the range of rejection of points around 0 via exclusionR
   //set the fit range around 0 with fitR
  
   fRange = fitR;
   fFitGaus = setGausFit;
   fExclRange  = exclusionR;
  
   if(fFitGaus) Printf("AliPerformancePtCalib::SetMakeFitOption: set MakeGausFit with fit range %2.3f and exclusion range in fabs(1/pt): %2.3f",fRange,fExclRange);
   else  Printf("AliPerformancePtCalib::SetMakeFitOption: set standard polynomial fit with fit range %2.3f and exclusion range in fabs(1/pt): %2.3f",fRange,fExclRange);
 
}
