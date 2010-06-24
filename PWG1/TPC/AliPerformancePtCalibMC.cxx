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
// Implementation of AliPerformancePtCalibMC class. It compares ESD, TPC track
// momenta with MC information to study charge/pt spectra.
// The output can be analysed with AliPerfAnalyzeInvPt via AliPerformancePtCalibMC::Analyse():
// Projection of charge/pt vs theta and vs phi resp. histoprams will be fitted with either
// polynomial or gaussian fit function to extract minimum position of charge/pt.
// Fit options and theta, phi bins can be set by user.
// Attention: use the Set functions of AliPerformancePtCalibMC when running
// AliPerformancePtCalibMC::Analyse()
// The result of the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliPerformancePtCalibMC.
//
// Author: S.Schuchmann 11/13/2009 
//------------------------------------------------------------------------------

/*
 
// after running the performance task, read the file, and get component

TFile f("Output.root");
AliPerformancePtCalibMC *compObj = (AliPerformancePtCalibMC*)coutput->FindObject("AliPerformancePtCalibMC");
 
// set phi and theta bins for fitting and analyse comparison data
compObj->SetProjBinsTheta(thetaBins,nThetaBins,minPhi, maxPhi);
compObj->SetProjBinsPhi(phiBins,nPhiBins,minTheta,maxTheta);
compObj->SetMakeFitOption(kFALSE,exclRange,fitRange);
compObj->SetDoRebin(rebin);
//compObj->SetAnaMCOff();
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

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDfriendTrack.h"
#include "AliESDfriend.h"

#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"

#include "AliPerformancePtCalibMC.h"
#include "AliPerfAnalyzeInvPt.h"


using namespace std;

ClassImp(AliPerformancePtCalibMC)

//________________________________________________________________________
   AliPerformancePtCalibMC::AliPerformancePtCalibMC() :
      AliPerformanceObject("AliPerformancePtCalibMC"),
      // option parameter for AliPerformancePtCalibMC::Analyse()
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
      fAnaMC(0),
      // option parameter for user defined charge/pt shift
      fShift(0),
      fDeltaInvP(0),
      //options for cuts
      fOptTPC(0),
      fESDcuts(0),
      fPhysSel(0),
      fEtaAcceptance(0),
      fTrigger(AliTriggerAnalysis::kMB1),
      fPhysicsSelection(0),
      fCutsRC(0),
      fCutsMC(0),
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
      fHistInvPtPtThetaPhiMC(0),
      fHistInvPtMCESD(0),
      fHistInvPtMCTPC(0),
      fHistPtMCESD(0),
      fHistPtMCTPC(0),
      fHistMomresMCESD(0),
      fHistMomresMCTPC(0),
      fHistTPCMomentaPosInvPtMC(0),
      fHistTPCMomentaNegInvPtMC(0),
      fHistTPCMomentaPosPtMC(0),
      fHistTPCMomentaNegPtMC(0),
      fHistESDMomentaPosInvPtMC(0),
      fHistESDMomentaNegInvPtMC(0),
      fHistESDMomentaPosPtMC(0), 
      fHistESDMomentaNegPtMC(0),
      fHistUserPtShift(0),
      
      //esd track cuts
      fESDTrackCuts(0),
      // analysis folder 
      fAnalysisFolder(0)
{
   // Dummy constructor
   
   
   fShift = kFALSE;                       // shift in charge/pt yes/no
   fDeltaInvP = 0.00;                     // shift value
   //options for cuts

   fOptTPC =  kTRUE;                      // read TPC tracks yes/no
   fESDcuts = kFALSE;
   fPhysSel = kFALSE;
   fCutsRC = NULL;
   fCutsMC = NULL;

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
   fAnaMC = kTRUE; // analyse MC tracks yes/no
   fDoRebin = kFALSE;
   fRebin = 0;
   
   Init();
} 

//________________________________________________________________________
AliPerformancePtCalibMC::AliPerformancePtCalibMC(const char *name= "AliPerformancePtCalibMC", const char *title="AliPerformancePtCalibMC"):
   AliPerformanceObject(name,title),
   // option parameter for AliPerformancePtCalibMC::Analyse()
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
      fAnaMC(0),
      // option parameter for user defined charge/pt shift
      fShift(0),
      fDeltaInvP(0),
      //options for cuts
      fOptTPC(0),
      fESDcuts(0),
      fPhysSel(0),
      fEtaAcceptance(0),
      fTrigger(AliTriggerAnalysis::kMB1),
      fPhysicsSelection(0),
      fCutsRC(0),
      fCutsMC(0),
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
      fHistInvPtPtThetaPhiMC(0),
      fHistInvPtMCESD(0),
      fHistInvPtMCTPC(0),
      fHistPtMCESD(0),
      fHistPtMCTPC(0),
      fHistMomresMCESD(0),
      fHistMomresMCTPC(0),
      fHistTPCMomentaPosInvPtMC(0),
      fHistTPCMomentaNegInvPtMC(0),
      fHistTPCMomentaPosPtMC(0),
      fHistTPCMomentaNegPtMC(0),
      fHistESDMomentaPosInvPtMC(0),
      fHistESDMomentaNegInvPtMC(0),
      fHistESDMomentaPosPtMC(0), 
      fHistESDMomentaNegPtMC(0),
      fHistUserPtShift(0),
      
      //esd track cuts
      fESDTrackCuts(0),
      // analysis folder 
      fAnalysisFolder(0)
{
   // Dummy constructor
   
   
   fShift = kFALSE;                       // shift in charge/pt yes/no
   fDeltaInvP = 0.00;                     // shift value
   //options for cuts

   fOptTPC =  kTRUE;                      // read TPC tracks yes/no
   fESDcuts = kFALSE;
   fPhysSel = kFALSE;
   fCutsRC = NULL;
   fCutsMC = NULL;

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
   fAnaMC = kTRUE; // analyse MC tracks yes/no
   fDoRebin = kFALSE;// flag for rebin
   fRebin = 0;// bins for rebin
   
   Init();
}

//________________________________________________________________________
AliPerformancePtCalibMC::~AliPerformancePtCalibMC() { 
   //
   // destructor
   //

   if(fList) delete fList;
   // histograms
   if(fHistInvPtPtThetaPhi)       delete fHistInvPtPtThetaPhi;fHistInvPtPtThetaPhi=0;
   if(fHistPtShift0)              delete fHistPtShift0;fHistPtShift0=0; 
   if(fHistPrimaryVertexPosX)     delete fHistPrimaryVertexPosX;fHistPrimaryVertexPosX=0; 
   if(fHistPrimaryVertexPosY)     delete fHistPrimaryVertexPosY;fHistPrimaryVertexPosY=0; 
   if(fHistPrimaryVertexPosZ)     delete fHistPrimaryVertexPosZ;fHistPrimaryVertexPosZ=0; 
   if(fHistTrackMultiplicity)     delete fHistTrackMultiplicity;fHistTrackMultiplicity=0; 
   if(fHistTrackMultiplicityCuts) delete fHistTrackMultiplicityCuts;fHistTrackMultiplicityCuts=0; 

   if(fHistTPCMomentaPosP)    delete fHistTPCMomentaPosP;fHistTPCMomentaPosP=0; 
   if(fHistTPCMomentaNegP)    delete fHistTPCMomentaNegP;fHistTPCMomentaNegP=0; 
   if(fHistTPCMomentaPosPt)   delete fHistTPCMomentaPosPt;fHistTPCMomentaPosPt=0; 
   if(fHistTPCMomentaNegPt)   delete fHistTPCMomentaNegPt ;fHistTPCMomentaNegPt=0; 
   if(fHistUserPtShift)       delete fHistUserPtShift;fHistUserPtShift=0;
   if(fHistInvPtPtThetaPhiMC) delete fHistInvPtPtThetaPhiMC;fHistInvPtPtThetaPhiMC=0;
   if(fHistInvPtMCESD)  delete fHistInvPtMCESD;fHistInvPtMCESD=0;
   if(fHistInvPtMCTPC)  delete fHistInvPtMCTPC;fHistInvPtMCTPC=0;
   if(fHistPtMCESD)     delete fHistPtMCESD;fHistPtMCESD=0;
   if(fHistPtMCTPC)     delete fHistPtMCTPC;fHistPtMCTPC=0;
   if(fHistMomresMCESD) delete fHistMomresMCESD;fHistMomresMCESD=0;
   if(fHistMomresMCTPC) delete fHistMomresMCTPC;fHistMomresMCTPC=0;
   if(fHistTPCMomentaPosInvPtMC) delete fHistTPCMomentaPosInvPtMC;fHistTPCMomentaPosInvPtMC=0;
   if(fHistTPCMomentaNegInvPtMC) delete fHistTPCMomentaNegInvPtMC;fHistTPCMomentaNegInvPtMC=0;
   if(fHistTPCMomentaPosPtMC)    delete fHistTPCMomentaPosPtMC;fHistTPCMomentaPosPtMC=0;
   if(fHistTPCMomentaNegPtMC)    delete fHistTPCMomentaNegPtMC;fHistTPCMomentaNegPtMC=0;
   if(fHistESDMomentaPosInvPtMC) delete fHistESDMomentaPosInvPtMC;fHistESDMomentaPosInvPtMC=0;
   if(fHistESDMomentaNegInvPtMC) delete fHistESDMomentaNegInvPtMC;fHistESDMomentaNegInvPtMC=0;
   if(fHistESDMomentaPosPtMC)    delete fHistESDMomentaPosPtMC;fHistESDMomentaPosPtMC=0;
   if(fHistESDMomentaNegPtMC)    delete fHistESDMomentaNegPtMC;fHistESDMomentaNegPtMC=0;
   

   
   //esd track cuts
   if(fESDTrackCuts)   delete fESDTrackCuts;
   if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;


   
}

//________________________________________________________________________
void AliPerformancePtCalibMC::Init() 
{
   // Create histograms
   // Called once
   
   fList = new TList();
 
   // init folder
   fAnalysisFolder = CreateFolder("folderPt_TPC","Analysis Pt Resolution Folder");
   fList->Add(fAnalysisFolder);
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
   fHistPtShift0   = new TH1F("fHistPtShift0","1/pt dN/pt vs. pt of ESD track  ",800,-20.0,20.0);
   fList->Add(fHistPtShift0);
   const   Int_t invPtDims = 4;
   fMaxPhi=6.5;
   fMinPhi=0.0;
   fMaxTheta=3.0;
   fMinTheta=0.0;
   Double_t xminInvPt[invPtDims] = {-4.5,-20.0,fMinTheta,fMinPhi};
   Double_t xmaxInvPt[invPtDims] = {4.5,20.0,fMaxTheta,fMaxPhi};
   Int_t  binsInvPt[invPtDims] = {900,800,300,325};
   
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
   
   // momentum test histos MC
   fHistTPCMomentaPosInvPtMC = new TH2F("fHistTPCMomentaPosInvPtMC","TPC-MC of 1/pt vs global ESD-MC of 1/pt pos",500, -10.0, 10.0,500, -10.0,10.0);
   fList->Add(fHistTPCMomentaPosInvPtMC);
   fHistTPCMomentaNegInvPtMC = new TH2F("fHistTPCMomentaNegInvPtMC","TPC-MC of 1/pt vs global ESD-MC 1/pt neg",500, -10.0, 10.0,500, -10.0, 10.0);
   fList->Add(fHistTPCMomentaNegInvPtMC);
   fHistTPCMomentaPosPtMC    = new TH2F("fHistTPCMomentaPosPtMC","TPC-MC of pt vs global ESD-MC of pt pos",600,-4.0,44.0,600,-4.0,44.0);
   fList->Add(fHistTPCMomentaPosPtMC);
   fHistTPCMomentaNegPtMC    = new TH2F("fHistTPCMomentaNegPtMC","TPC-MC of pt vs global ESD-MC of pt neg",600,-4.0,44.0,600,-4.0,44.0);
   fList->Add(fHistTPCMomentaNegPtMC);
   fHistESDMomentaPosInvPtMC = new TH1F("fHistESDMomentaPosInvPtMC","ESD-MC of 1/pt ",500, -10.0, 10.0);
   fList->Add(fHistESDMomentaPosInvPtMC);
   fHistESDMomentaNegInvPtMC = new TH1F("fHistESDMomentaNegInvPtMC","ESD-MC of 1/pt",500, -10.0, 10.0);
   fList->Add(fHistESDMomentaNegInvPtMC);
   fHistESDMomentaPosPtMC    = new TH1F("fHistESDMomentaPosPtMC","ESD-MC of pt ",600,-4.0,44.0);
   fList->Add(fHistESDMomentaPosPtMC);
   fHistESDMomentaNegPtMC    = new TH1F("fHistESDMomentaNegPtMC","ESD-MC of pt ",600,-4.0,44.0);
   fList->Add(fHistESDMomentaNegPtMC);

   // MC only info
   fHistInvPtPtThetaPhiMC = new THnSparseF("fHistInvPtPtThetaPhiMC","MC 1/pt vs pt vs #theta vs #phi ",invPtDims,binsInvPt,xminInvPt,xmaxInvPt);
   fList->Add(fHistInvPtPtThetaPhiMC);

 
   //correlation histos MC ESD or TPC
   fHistInvPtMCESD  = new TH2F("fHistInvPtMCESD","inv pt ESD vs MC",900, 0.0, 9.0,900, 0.0, 9.0);
   fList->Add(fHistInvPtMCESD);
   fHistPtMCESD     = new TH2F("fHistPtMCESD"," pt ESD vs MC",300, 0.0, 15.0,300, 0.0, 15.0);
   fList->Add(fHistPtMCESD);
   fHistInvPtMCTPC  = new TH2F("fHistInvPtMCTPC","inv pt TPC vs MC",900, 0.0, 9.0,900, 0.0, 9.0);
   fList->Add(fHistInvPtMCTPC);
   fHistPtMCTPC     = new TH2F("fHistPtMCTPC"," pt TPC vs MC",300, 0.0, 15.0,300, 0.0, 15.0);
   fList->Add(fHistPtMCTPC);
   fHistMomresMCESD = new TH2F("fHistMomresMCESD"," (pt ESD - pt MC)/ptMC vs pt MC",300, 0.0, 15.0,400, -2.0, 2.0);
   fList->Add(fHistMomresMCESD);
   fHistMomresMCTPC = new TH2F("fHistMomresMCTPC"," (pt TPC - pt MC)/ptMC vs pt MC",300, 0.0, 15.0,400, -2.0, 2.0);
   fList->Add(fHistMomresMCTPC);


   //user pt shift check
   fHistUserPtShift = new TH1F("fHistUserPtShift","user defined shift in 1/pt",100,-0.5,1.5);
   fList->Add(fHistUserPtShift);
   
   // esd track cuts  
   fESDTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
   
  
}

//________________________________________________________________________
void AliPerformancePtCalibMC::SetPtShift(const Double_t shiftVal ) {

   //set user defined shift in charge/pt
   
   if(shiftVal) { fShift=kTRUE; fDeltaInvP = shiftVal; } 
}

//________________________________________________________________________
void AliPerformancePtCalibMC::Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, AliESDfriend *const /*esdFriend*/, const Bool_t /*bUseMC*/, const Bool_t /*bUseESDfriend*/)
{
   //exec: read MC and esd or tpc tracks
   
   AliStack* stack = NULL;
 
   if (!esdEvent) {
      Printf("ERROR: Event not available");
      return;
   }

   if (!(esdEvent->GetNumberOfTracks())) {
      Printf(" PtCalibMC task: There is no track in this event");
      return;
   }
   fHistTrackMultiplicity->Fill(esdEvent->GetNumberOfTracks());

   if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
   }    
   stack = mcEvent->Stack();
   if (!stack) {
      Printf("ERROR: Could not retrieve stack");
      return;
   }

  if (GetTriggerClass()){
      Bool_t isEventTriggered = esdEvent->IsTriggerClassFired(GetTriggerClass());
      if(!isEventTriggered) return;
   }
   
  /*
   // trigger selection
   Bool_t isEventTriggered = kTRUE;
   AliPhysicsSelection *trigSel = NULL;
   AliTriggerAnalysis *trigAna = NULL;
   
   if(fPhysSel){ trigSel = GetPhysicsTriggerSelection();
   if(!trigSel) {
   Printf("cannot get trigSel \n");
   return;
   }
   }
   
//   if(esdEvent) isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);// crashes
   
   if(fTrigger == AliTriggerAnalysis::kV0AND) 
   {Printf("physics selection V0AND");//test
   trigAna = trigSel->GetTriggerAnalysis();
   if(!trigAna) 
   return;
   
   isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent,GetTrigger());
	
   }
   
   if(fPhysSel && !isEventTriggered) return;
    
  
   //vertex info for cut
   const AliESDVertex *vtx = esdEvent->GetPrimaryVertex();
   if (!vtx->GetStatus()) return ;
   */


   
   if(fShift) fHistUserPtShift->Fill(fDeltaInvP);
  
   // read primary vertex info
   Double_t tPrimaryVtxPosition[3];
   // Double_t tPrimaryVtxCov[3];
   const AliESDVertex *primaryVtx = esdEvent->GetPrimaryVertexTPC();
 
   tPrimaryVtxPosition[0] = primaryVtx->GetXv();
   tPrimaryVtxPosition[1] = primaryVtx->GetYv();
   tPrimaryVtxPosition[2] = primaryVtx->GetZv();
  
   fHistPrimaryVertexPosX->Fill(tPrimaryVtxPosition[0]);
   fHistPrimaryVertexPosY->Fill(tPrimaryVtxPosition[1]);
   fHistPrimaryVertexPosZ->Fill(tPrimaryVtxPosition[2]);
 

   //fill histos for pt spectra and shift of transverse momentum
   Int_t count=0;
 
   for(Int_t j = 0;j<esdEvent->GetNumberOfTracks();j++){
      AliESDtrack *esdTrack = esdEvent->GetTrack(j);
      if(!esdTrack) continue;

      //esd track cuts
      if(fESDcuts){
	 if(!fESDTrackCuts->AcceptTrack(esdTrack)) continue;
      }
      
      // get MC info 
      Int_t esdLabel = esdTrack->GetLabel();
      if(esdLabel<0) continue;	
      TParticle *  partMC = stack->Particle(esdLabel);
      if (!partMC) continue;
  
      // fill correlation histos MC ESD
      Double_t pESD  = esdTrack->GetP();
      Double_t ptESD = esdTrack->GetSignedPt();
    
      if(!ptESD || !(partMC->Pt()) ) continue;
      Double_t mcPt = partMC->Pt();
      Double_t invPtMC = 1.0/mcPt;
      Int_t signMC = partMC->GetPdgCode();
      //MC only
      if(signMC>0) signMC = 1; 
      else signMC = -1;

      //fill MC information
      Double_t thetaMC = partMC->Theta();
      Double_t phiMC = partMC->Phi();
      
      Double_t momAngMC[4] = {signMC*(fabs(invPtMC)),signMC*(fabs(mcPt)),thetaMC,phiMC};

      // fill only if MC track is in eta acceptance of TPC in order to be compareable to TPC tracks
      if(fabs( partMC->Eta())<= fEtaAcceptance) {
	 fHistInvPtPtThetaPhiMC->Fill(momAngMC);
	 
	 //correlation histos MC ESD
	 fHistInvPtMCESD->Fill(fabs(invPtMC),fabs(1.0/ptESD));
	 fHistPtMCESD->Fill(fabs(mcPt),fabs(ptESD));
      }

      // fill histos TPC or ESD
      if(fOptTPC){
	 //TPC tracks and MC tracks
	 const AliExternalTrackParam *tpcTrack = esdTrack->GetTPCInnerParam(); 
	 if(!tpcTrack) continue;
	 if(fabs(tpcTrack->Eta())>  fEtaAcceptance) continue;
      
	 Double_t signedPt = tpcTrack->GetSignedPt();
	 Double_t invPt = 0.0;
	 if(signedPt) {
	    invPt = 1.0/signedPt;
	
	    fHistPtShift0->Fill(signedPt);//changed

	    if(fShift){Printf("user shift of momentum SET to non zero value!");
	       invPt += fDeltaInvP; //shift momentum for tests
	       if(invPt) signedPt = 1.0/invPt;
	       else continue;
	    }


            Double_t theta = tpcTrack->Theta();
            Double_t phi = tpcTrack->Phi();

	    Double_t momAng[4] = {invPt,signedPt,theta,phi};
	    fHistInvPtPtThetaPhi->Fill(momAng);

	    //correlation histos MC TPC
	    fHistInvPtMCTPC->Fill(fabs(invPtMC),fabs(invPt));
	    fHistPtMCTPC->Fill(fabs(mcPt),fabs(signedPt));
	
	    //compare to MC info
	    Double_t  ptDiffESD = (fabs(ptESD)-fabs(mcPt))/pow(mcPt,2);
	    Double_t  ptDiffTPC = (fabs(signedPt)-fabs(mcPt))/pow(mcPt,2);
	    Double_t  invPtDiffESD = fabs(1.0/ptESD)-1.0/fabs(mcPt);
	    Double_t  invPtDiffTPC = fabs(invPt)-1.0/fabs(mcPt);
	    Double_t  pTPC  = tpcTrack->GetP();
	
	    if(esdTrack->GetSign()>0){//compare momenta ESD track and TPC track
	       fHistTPCMomentaPosP->Fill(fabs(pESD),fabs(pTPC));
	       fHistTPCMomentaPosPt->Fill(fabs(ptESD),fabs(signedPt));
	       fHistTPCMomentaPosInvPtMC->Fill(invPtDiffESD,invPtDiffTPC);
	       fHistTPCMomentaPosPtMC->Fill(ptDiffESD,ptDiffTPC);
	    }
	    else{
	       fHistTPCMomentaNegP->Fill(fabs(pESD),fabs(pTPC));
	       fHistTPCMomentaNegPt->Fill(fabs(ptESD),fabs(signedPt));
	       fHistTPCMomentaNegInvPtMC->Fill(invPtDiffESD,invPtDiffTPC);
	       fHistTPCMomentaNegPtMC->Fill(ptDiffESD,ptDiffTPC);
	    }
	    fHistMomresMCESD->Fill((fabs(mcPt)-fabs(ptESD))/fabs(mcPt),fabs(mcPt));
	    fHistMomresMCTPC->Fill((fabs(mcPt)-fabs(signedPt))/fabs(mcPt),fabs(mcPt));
	    count++;
	 }
	 else continue;
      }
   
      else{
	 // ESD tracks and MC tracks
	 if(fabs(esdTrack->Eta())> fEtaAcceptance) continue;
	 Double_t invPt = 0.0;
      
	 if(ptESD) {
	    invPt = 1.0/ptESD; 
	    fHistPtShift0->Fill(ptESD);//changed
	
	    if(fShift){Printf("user shift of momentum SET to non zero value!");
	       invPt += fDeltaInvP; //shift momentum for tests
	       if(invPt) ptESD = 1.0/invPt; 
	       else continue;
	    }

	    Double_t theta = esdTrack->Theta();
	    Double_t phi = esdTrack->Phi();

	    Double_t momAng[4] = {invPt,ptESD,theta,phi};
	    fHistInvPtPtThetaPhi->Fill(momAng);

	    //differences MC ESD tracks
	    Double_t ptDiffESD = (fabs(ptESD)-fabs(mcPt))/pow(mcPt,2);
	    Double_t invPtdiffESD = fabs(1.0/ptESD)-1.0/fabs(mcPt);
	    if(esdTrack->GetSign()>0){   
	       fHistESDMomentaPosInvPtMC->Fill(invPtdiffESD);
	       fHistESDMomentaPosPtMC->Fill(ptDiffESD);
	    }
	    else{
	       fHistESDMomentaNegInvPtMC->Fill(invPtdiffESD);
	       fHistESDMomentaNegPtMC->Fill(ptDiffESD);
	    }	
	
	    fHistMomresMCESD->Fill((fabs(mcPt)-fabs(ptESD))/fabs(mcPt),fabs(mcPt));
	    count++;
	 }
      }
   }
    
   fHistTrackMultiplicityCuts->Fill(count);
  
}    

//______________________________________________________________________________________________________________________

void AliPerformancePtCalibMC::Analyse()
{
  
   // analyse charge/pt spectra in bins of theta and phi. Bins can be set by user
   
 
   THnSparseF *copyTHnSparseTheta;
   THnSparseF *copyTHnSparsePhi;
   
   if(fAnaMC){
      Printf("AliPerformancePtCalibMC::Analyse: analysing MC!");
      copyTHnSparseTheta = (THnSparseF*)fHistInvPtPtThetaPhiMC->Clone("copyTHnSparseThetaMC");
      copyTHnSparsePhi   = (THnSparseF*)fHistInvPtPtThetaPhiMC->Clone("copyTHnSparsePhiMC");
   }
   else {
      Printf("AliPerformancePtCalibMC::Analyse: analysing ESD (reco)!");
      copyTHnSparseTheta = (THnSparseF*)fHistInvPtPtThetaPhi->Clone("copyTHnSparseTheta");
      copyTHnSparsePhi   = (THnSparseF*)fHistInvPtPtThetaPhi->Clone("copyTHnSparsePhi");
   }
   
   copyTHnSparseTheta->GetAxis(3)->SetRangeUser(fMinPhi,fMaxPhi);
   copyTHnSparsePhi->GetAxis(2)->SetRangeUser(fMinTheta,fMaxTheta);

   TH2F *histInvPtTheta = (TH2F*)copyTHnSparseTheta->Projection(2,0);
   TH2F *histInvPtPhi   = (TH2F*)copyTHnSparsePhi->Projection(3,0);
   
   AliPerfAnalyzeInvPt *ana = new  AliPerfAnalyzeInvPt("AliPerfAnalyzeInvPt","AliPerfAnalyzeInvPt");
  
   TH1::AddDirectory(kFALSE);
 
   ana->SetProjBinsTheta(fThetaBins,fNThetaBins);
   ana->SetProjBinsPhi(fPhiBins,fNPhiBins);
   ana->SetMakeFitOption(fFitGaus,fExclRange,fRange);
   if(fDoRebin) ana->SetDoRebin(fRebin);
   TObjArray *aFolderObj = new TObjArray;

   ana->StartAnalysis(histInvPtTheta,histInvPtPhi, aFolderObj);

   // export objects to analysis folder
   fAnalysisFolder = ExportToFolder(aFolderObj);

   // delete only TObjArray
   if(aFolderObj) delete aFolderObj;
   if(ana) delete ana;
  
}

//______________________________________________________________________________________________________________________
TFolder* AliPerformancePtCalibMC::ExportToFolder(TObjArray * array) 
{
   // recreate folder avery time and export objects to new one
   //
   AliPerformancePtCalibMC * comp=this;
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
Long64_t AliPerformancePtCalibMC::Merge(TCollection* const list) 
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
	 AliPerformancePtCalibMC* entry = dynamic_cast<AliPerformancePtCalibMC*>(obj);
	 if (!entry) continue; 
	 fHistInvPtPtThetaPhi->Add(entry->fHistInvPtPtThetaPhi);

	 fHistInvPtPtThetaPhiMC->Add(entry->fHistInvPtPtThetaPhiMC);

	 fHistInvPtMCESD->Add(entry->fHistInvPtMCESD);
	 fHistPtMCESD->Add(entry->fHistPtMCESD);
	 fHistInvPtMCTPC->Add(entry->fHistInvPtMCTPC);
	 fHistPtMCTPC->Add(entry->fHistPtMCTPC);
	 fHistMomresMCESD->Add(entry->fHistMomresMCESD);
	 fHistMomresMCTPC->Add(entry->fHistMomresMCTPC);

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
	 fHistTPCMomentaPosInvPtMC->Add(entry->fHistTPCMomentaPosInvPtMC);
	 fHistTPCMomentaNegInvPtMC->Add(entry->fHistTPCMomentaNegInvPtMC);
	 fHistTPCMomentaPosPtMC->Add(entry->fHistTPCMomentaPosPtMC);
	 fHistTPCMomentaNegPtMC->Add(entry->fHistTPCMomentaNegPtMC);
	 fHistESDMomentaPosInvPtMC->Add(entry->fHistESDMomentaPosInvPtMC);
	 fHistESDMomentaNegInvPtMC->Add(entry->fHistESDMomentaNegInvPtMC);
	 fHistESDMomentaPosPtMC->Add(entry->fHistESDMomentaPosPtMC);
	 fHistESDMomentaNegPtMC->Add(entry->fHistESDMomentaNegPtMC);
	 count++;
      }
  
   return count;
}

//______________________________________________________________________________________________________________________
TFolder* AliPerformancePtCalibMC::CreateFolder(TString name,TString title) { 
   // create folder for analysed histograms
   //
   TFolder *folder = 0;
   folder = new TFolder(name.Data(),title.Data());

   return folder;
}


// set variables for Analyse()
//______________________________________________________________________________________________________________________
void AliPerformancePtCalibMC::SetProjBinsPhi(const Double_t *phiBinArray,const Int_t nphBins, const  Double_t minTheta, const  Double_t maxTheta){
   // set phi bins for Analyse()
   // set phi bins as array and set number of this array which is equal to number of bins analysed
   // the last analysed bin will always be the projection from first to last bin in the array
   if(nphBins){
      fNPhiBins = nphBins;
  
      for(Int_t k = 0;k<fNPhiBins;k++){
	 fPhiBins[k] = phiBinArray[k];
      }
      Printf("AliPerformancePtCalibMC::SetProjBinsPhi: number of bins in phi set to %i",fNPhiBins);
   }
   else  Printf("Warning AliPerformancePtCalibMC::SetProjBinsPhi: number of bins in phi NOT set!!! Default values are taken.");

   if(fabs(minTheta-maxTheta)<0.001){
      Printf("AliPerformancePtCalibMC::SetProjBinsPhi: theta range for projection for projection on phi and charge/pt is too small. whole range of theta selected.");
   }
   else{
      Printf("AliPerformancePtCalibMC::SetProjBinsPhi: theta range for projection on phi and charge/pt is selected by user: %1.3f - %1.3f rad",minTheta,maxTheta);
      fMinTheta = minTheta;
      fMaxTheta = maxTheta;
   }
}
//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalibMC::SetProjBinsTheta(const Double_t *thetaBinArray, const Int_t nthBins, const Double_t minPhi, const Double_t maxPhi){
   // set theta bins for Analyse()
   // set theta bins as array and set number of this array which is equal to number of bins analysed
   // the last analysed bin will always be the projection from first to last bin in the array
   if(nthBins){
      fNThetaBins = nthBins;
      for(Int_t k = 0;k<fNThetaBins;k++){
	 fThetaBins[k] = thetaBinArray[k];
      }
      Printf("AliPerformancePtCalibMC::SetProjBinsTheta: number of bins in theta set to %i",fNThetaBins);
   }
   else  Printf("Warning AliPerformancePtCalibMC::SetProjBinsTheta: number of bins in theta NOT set!!! Default values are taken.");
   
   if(fabs(minPhi-maxPhi)<0.001){
      Printf("AliPerformancePtCalibMC::SetProjBinsTheta: phi range for projection for projection on theta and charge/pt is too small. whole range of phi selected.");
   }
   else{
      Printf("AliPerformancePtCalibMC::SetProjBinsTheta: phi range for projection on theta and charge/pt is selected by user: %1.3f - %1.3f rad",minPhi,maxPhi);
      fMinPhi = minPhi;
      fMaxPhi = maxPhi;
   }
}

//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalibMC::SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR ){

   // set the fit options:
   // for usage of gaussian function instead of polynomial (default) set setGausFit=kTRUE
   // set the range of rejection of points around 0 via exclusionR
   // set the fit range around 0 with fitR
   
   fFitGaus = setGausFit;
   fExclRange  = exclusionR;
   fRange = fitR;
  
   if(fFitGaus) Printf("AliPerformancePtCalibMC::SetMakeFitOption: set MakeGausFit with fit range %2.3f and exclusion range in fabs(1/pt): %2.3f",fRange,fExclRange);
   else  Printf("AliPerformancePtCalibMC::SetMakeFitOption: set standard polynomial fit with fit range %2.3f and exclusion range in fabs(1/pt): %2.3f",fRange,fExclRange);
 
}
