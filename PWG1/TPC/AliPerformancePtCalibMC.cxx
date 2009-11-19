//------------------------------------------------------------------------------
// Implementation of AliPerformancePtCalibMC class. It compares ESD, TPC track
// momenta with MC information.
// The output can be analysed with AliPerfAnalyzeInvPt* via AliPerformancePtCalibMC::Analyse():
// Projection of 1/pt vs theta and vs phi resp. histoprams will be fitted with either
// polynomial or gaussian fit function to extract minimum position of 1/pt.
// Fit options and theta, phi bins can be set by user.
// Attention: use the Set* functions of AliPerformancePtCalibMC when running
// AliPerformancePtCalibMC::Analyse()
// The result of the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliPerformancePtCalib*.
//
// Author: S.Schuchmann 11/13/2009 
//------------------------------------------------------------------------------

/*
 
// after running comparison task, read the file, and get component
gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
LoadMyLibs();

TFile f("Output.root");
AliPerformancePtCalib * compObj = (AliPerformancePtCalibMC*)coutput->FindObject("AliPerformancePtCalibMC");
 
// analyse comparison data
compObj->Analyse();

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
      fRange(0),
      fExclRange(0),
      fFitGaus(0) ,
      fAnaMC(0),
      // option parameter for user defined charge/pt shift
      fShift(0),
      fDeltaInvP(0),
      //options for cuts
      fOptTPC(0),
      fESDcuts(0),
      fRefitTPC(0),
      fRefitITS(0),
      fDCAcut(0),
      fEtaAcceptance(0),
      fMinPt(0),
      fMaxPt(0),
      fMinNClustersTPC(0),
      fMaxChi2PerClusterTPC(0),
      fMaxDCAtoVertexXY(0),
      fMaxDCAtoVertexZ(0),
      fAcceptKinkDaughters(0),
      fRequireSigmaToVertex(0),
      fDCAToVertex2D(0),
      
      fCutsRC(0),
      fCutsMC(0),
      
      fList(0),
      // histograms
      fHistInvPtTheta(0),
      fHistInvPtPhi(0),
      fHistPtTheta(0),
      fHistPtPhi(0),
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
      fHistInvPtThetaMC(0),
      fHistInvPtPhiMC(0),
      fHistPtThetaMC(0),
      fHistPtPhiMC(0),
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
      fESDTrackCuts(0),
      // analysis folder 
      fAnalysisFolder(0)
{
   // Dummy constructor
   
   
   fShift = kFALSE;                       // shift in charge/pt yes/no
   fDeltaInvP = 0.00;                     // shift value
   //options for cuts
   fOptTPC =  kTRUE;                      // read TPC tracks yes/no
   fESDcuts = kTRUE;                      // read ESD track cuts
   fRefitTPC = kFALSE;                    // require TPC refit
   fRefitITS = kFALSE;                    // require ITS refit
   fDCAcut = kTRUE;                       // apply DCA cuts
   
   fCutsRC = NULL;
   fCutsMC = NULL;
   fEtaAcceptance = 0.8;
   fMinPt=0.15;  // GeV/c
   fMaxPt=1.e10; // GeV/c 
   fMinNClustersTPC = 50;
   fMaxChi2PerClusterTPC = 4.0;
   fMaxDCAtoVertexXY = 2.4; // cm
   fMaxDCAtoVertexZ  = 3.0; // cm
   fAcceptKinkDaughters = kFALSE;
   fRequireSigmaToVertex = kFALSE;
   fDCAToVertex2D = kTRUE;
   
   // options for function AliPerformancePtCalibMC::Analyse()
   fFitGaus = kFALSE;// use gaussian function for fitting charge/pt yes/no
   fNThetaBins = 0; //number of theta bins
   fNPhiBins = 0; //number of phi bins
   fRange = 0; //fit range around 0
   fExclRange =0; //range of rejection of points around 0
   fAnaMC = kTRUE; // analyse MC tracks yes/no
   
   Init();
} 

//________________________________________________________________________
AliPerformancePtCalibMC::AliPerformancePtCalibMC(const char *name= "AliPerformancePtCalibMC", const char *title="AliPerformancePtCalibMC")://,Int_t analysisMode=0,Bool_t hptGenerator=kFALSE) :
   AliPerformanceObject(name,title),
   // option parameter for AliPerformancePtCalibMC::Analyse()
   fNThetaBins(0), 
   fNPhiBins(0),
   fRange(0),
   fExclRange(0),
   fFitGaus(0) ,
   fAnaMC(0),
   // option parameter for user defined 1/pt shift
   fShift(0),
   fDeltaInvP (0),
   //options for cuts
   fOptTPC(0),
   fESDcuts(0),
   fRefitTPC(0),
   fRefitITS(0),
   fDCAcut(0),
   fEtaAcceptance(0),
   fMinPt(0),
   fMaxPt(0),
   fMinNClustersTPC(0),
   fMaxChi2PerClusterTPC(0),
   fMaxDCAtoVertexXY(0),
   fMaxDCAtoVertexZ(0),
   fAcceptKinkDaughters(0),
   fRequireSigmaToVertex(0),
   fDCAToVertex2D(0),
   
   fCutsRC(0),
   fCutsMC(0),

   fList(0),
   // histograms
   fHistInvPtTheta(0),
   fHistInvPtPhi(0),
   fHistPtTheta(0),
   fHistPtPhi(0),
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
   fHistInvPtThetaMC(0),
   fHistInvPtPhiMC(0),
   fHistPtThetaMC(0),
   fHistPtPhiMC(0),
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
   
   fESDTrackCuts(0),
   // analysis folder 
   fAnalysisFolder(0)
   
{
   // Constructor
    
   fShift = kFALSE;                       // shift in charge/pt yes/no
   fDeltaInvP = 0.00;                     // shift value
   //options for cuts
   fOptTPC =  kTRUE;                      // read TPC tracks yes/no
   fESDcuts = kTRUE;                      // read ESD track cuts
   fRefitTPC = kFALSE;                    // require TPC refit
   fRefitITS = kFALSE;                    // require ITS refit
   fDCAcut = kTRUE;                       // apply DCA cuts
 
   fCutsRC = NULL;
   fCutsMC = NULL;
   fEtaAcceptance = 0.8;
   fMinPt=0.15;  // GeV/c
   fMaxPt=1.e10; // GeV/c 
   fMinNClustersTPC = 50;
   fMaxChi2PerClusterTPC = 4.0;
   fMaxDCAtoVertexXY = 2.4; // cm
   fMaxDCAtoVertexZ  = 3.0; // cm
   fAcceptKinkDaughters = kFALSE;
   fRequireSigmaToVertex = kFALSE;
   fDCAToVertex2D = kTRUE;
   
   // options for function AliPerformancePtCalibMC::Analyse()
   fFitGaus = kFALSE;// use gaussian function for fitting charge/pt yes/no
   fNThetaBins = 0; //number of theta bins
   fNPhiBins = 0; //number of phi bins
   fRange = 0; //fit range around 0
   fExclRange =0; //range of rejection of points around 0
   fAnaMC = kTRUE; // analyse MC tracks yes/no

   Init();
}

//________________________________________________________________________
AliPerformancePtCalibMC::~AliPerformancePtCalibMC() { 
   //
   // destructor
   //
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
   fHistPtShift0   = new TH1F("fHistPtShift0","1/pt dN/pt vs. pt of ESD track  ",600,0.0,6.0);
   fList->Add(fHistPtShift0);
   fHistInvPtTheta = new TH2F("fHistInvPtTheta","#theta vs 1/pt ",900, -4.5, 4.5,300,0.0,3.0);
   fList->Add(fHistInvPtTheta);
   fHistInvPtPhi   = new TH2F("fHistInvPtPhi","#phi vs 1/pt",900, -4.5, 4.5,325,0.0,6.5);
   fList->Add(fHistInvPtPhi);
   fHistPtTheta    = new TH2F("fHistPtTheta"," #theta vs pt ",300, 0.0, 15.0,300,0.0,3.0);
   fList->Add(fHistPtTheta);
   fHistPtPhi      = new TH2F("fHistPtPhi"," #phi vs  pt ",300, 0.0,15.0,325,0.0,6.5);
   fList->Add(fHistPtPhi);
  
   // mom test histos
   fHistTPCMomentaPosP  =  new TH2F("fHistTPCMomentaPosP","TPC p vs global esd track p pos",300,0.0,15.0,300,0.0,15.0);
   fList->Add(fHistTPCMomentaPosP);
   fHistTPCMomentaNegP  =  new TH2F("fHistTPCMomentaNegP","TPC p vs global esd track p neg",300,0.0,15.0,300,0.0,15.0);
   fList->Add(fHistTPCMomentaNegP);
   fHistTPCMomentaPosPt =  new TH2F("fHistTPCMomentaPosPt","TPC pt vs global esd track pt pos",300,0.0,15.0,300,0.0,15.0);
   fList->Add(fHistTPCMomentaPosPt);
   fHistTPCMomentaNegPt =  new TH2F("fHistTPCMomentaNegPt","TPC pt vs global esd track pt neg",300,0.0,15.0,300,0.0,15.0);
   fList->Add(fHistTPCMomentaNegPt);
   
   // mom test histos MC
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

   // MC info
   fHistInvPtThetaMC = new TH2F("fHistInvPtThetaMC","theta vs inv pt MC",900, -4.5, 4.5,300,0.0,3.0);
   fList->Add(fHistInvPtThetaMC);
   fHistInvPtPhiMC   = new TH2F("fHistInvPtPhiMC","phi vs inv pt MC",900, -4.5, 4.5,325,0.0,6.5);
   fList->Add(fHistInvPtPhiMC);
   fHistPtThetaMC    = new TH2F("fHistPtThetaMC","theta vs pt MC",300, 0.0, 15.0,300,0.0,3.0);
   fList->Add(fHistPtThetaMC);
   fHistPtPhiMC      = new TH2F("fHistPtPhiMC"," phi vs pt MC",300, 0.0,15.0,325,0.0,6.5);
   fList->Add(fHistPtPhiMC);
 
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
  
   //   //fESDTrackCuts->DefineHistoqgrams(1);
 
   fESDTrackCuts->SetRequireSigmaToVertex(fRequireSigmaToVertex);
   fESDTrackCuts->SetRequireTPCRefit(fRefitTPC);
   fESDTrackCuts->SetAcceptKinkDaughters(fAcceptKinkDaughters);
   fESDTrackCuts->SetMinNClustersTPC((Int_t)fMinNClustersTPC);
   fESDTrackCuts->SetMaxChi2PerClusterTPC(fMaxChi2PerClusterTPC);
   fESDTrackCuts->SetMaxDCAToVertexXY(fMaxDCAtoVertexXY);
   fESDTrackCuts->SetMaxDCAToVertexZ(fMaxDCAtoVertexZ);
   fESDTrackCuts->SetDCAToVertex2D(fDCAToVertex2D);
   fESDTrackCuts->SetPtRange(fMinPt,fMaxPt); 
 
  
}

//________________________________________________________________________
void AliPerformancePtCalibMC::SetPtShift(const Double_t shiftVal ) {

   //set user defined shift in charge/pt
   
   if(shiftVal) { fShift=kTRUE; fDeltaInvP = shiftVal; } 
}

//________________________________________________________________________
void AliPerformancePtCalibMC::Exec(AliMCEvent* const mcEvent, AliESDEvent* const esdEvent , AliESDfriend*, Bool_t, Bool_t) 
{
   //exec: read MC and esd or tpc tracks
   
   AliStack* stack;
 
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
      if(fESDcuts){Printf("esd cuts aplied");
	 if(!fESDTrackCuts->AcceptTrack(esdTrack)) continue;
      }
    
      //more track cuts
      if(fRefitTPC) if(AddTPCcuts(esdTrack)) continue;
      if(fRefitITS) if(AddITScuts(esdTrack)) continue;
      if(fDCAcut)   if(AddDCAcuts(esdTrack)) continue ;
    
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

      //fill MC histos
      fHistInvPtThetaMC->Fill(signMC*fabs(invPtMC),partMC->Theta());
      fHistInvPtPhiMC->Fill(signMC*fabs(invPtMC),partMC->Phi());
      fHistPtThetaMC->Fill(fabs(mcPt),partMC->Theta(),fabs(invPtMC));
      fHistPtPhiMC->Fill(fabs(mcPt),partMC->Phi(),fabs(invPtMC));

      //correlation histos MC ESD
      fHistInvPtMCESD->Fill(fabs(invPtMC),fabs(1.0/ptESD));
      fHistPtMCESD->Fill(fabs(mcPt),fabs(ptESD));


      // fill histos
      if(fOptTPC){
	 //TPC tracks and MC tracks
	 const AliExternalTrackParam *tpcTrack = esdTrack->GetTPCInnerParam(); 
	 if(!tpcTrack) continue;
	 if(fabs(tpcTrack->Eta())>  fEtaAcceptance) continue;
      
	 Double_t signedPt = tpcTrack->GetSignedPt();
	 Double_t invPt = 0.0;
	 if(signedPt) {
	    invPt = 1.0/signedPt;
	
	    fHistPtShift0->Fill(fabs(signedPt));

	    if(fShift){
	       invPt += fDeltaInvP; //shift momentum for tests
	       if(invPt) signedPt = 1.0/invPt;
	       else continue;
	    }

	    fHistInvPtTheta->Fill(invPt,tpcTrack->Theta());
	    fHistInvPtPhi->Fill(invPt,tpcTrack->Phi());
	    fHistPtTheta->Fill(fabs(signedPt),tpcTrack->Theta());
	    fHistPtPhi->Fill(fabs(signedPt),tpcTrack->Phi());


	    //correlation histos MC TPC
	    fHistInvPtMCTPC->Fill(fabs(invPtMC),fabs(invPt));
	    fHistPtMCTPC->Fill(fabs(mcPt),fabs(signedPt));
	
	    //compare to MC info
	    Double_t  ptDiffESD = (fabs(ptESD)-fabs(mcPt))/pow(mcPt,2);
	    Double_t  ptDiffTPC = (fabs(signedPt)-fabs(mcPt))/pow(mcPt,2);
	    Double_t  invPtDiffESD = fabs(1.0/ptESD)-1.0/fabs(mcPt);
	    Double_t  invPtDiffTPC = fabs(invPt)-1.0/fabs(mcPt);
	    Double_t pTPC  = tpcTrack->GetP();
	
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
	 Double_t invPt = 0.0;
      
	 if(ptESD) {
	    invPt = 1.0/ptESD; 
	    fHistPtShift0->Fill(fabs(ptESD));
	
	    if(fShift){
	       invPt += fDeltaInvP; //shift momentum for tests
	       if(invPt) ptESD = 1.0/invPt; 
	       else continue;
	    }
	    fHistInvPtTheta->Fill(invPt,esdTrack->Theta());
	    fHistInvPtPhi->Fill(invPt,esdTrack->Phi());
	    fHistPtTheta->Fill(ptESD,esdTrack->Theta());
	    fHistPtPhi->Fill(ptESD,esdTrack->Phi());

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
Bool_t AliPerformancePtCalibMC::AddTPCcuts(const AliESDtrack *esdTrack){
   // apply TPC cuts
   
   Bool_t cut = kFALSE;
  
   if (!(esdTrack->GetStatus()&AliESDtrack::kTPCrefit)) cut=kTRUE; // TPC refit
   if (esdTrack->GetTPCNcls()<fMinNClustersTPC) cut=kTRUE; // min. nb. TPC clusters
   if(cut) return kTRUE;
   return kFALSE;
}
//______________________________________________________________________________________________________________________
Bool_t AliPerformancePtCalibMC::AddDCAcuts(const AliESDtrack *esdTrack){
   //apply DCA cuts
   Bool_t cut = kFALSE;
  
   Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z and impact parameters:
   esdTrack->GetImpactParameters(dca,cov);
   if(TMath::Abs(dca[0])>fMaxDCAtoVertexXY || TMath::Abs(dca[1])>fMaxDCAtoVertexZ) cut=kTRUE;
   if(esdTrack->GetKinkIndex(0)>0) cut=kTRUE;
   if(cut) return kTRUE;
   return kFALSE;
}

//______________________________________________________________________________________________________________________
Bool_t AliPerformancePtCalibMC::AddITScuts(const AliESDtrack *esdTrack){
   //apply ITS cuts
   Bool_t cut = kFALSE;
  
   if (!(esdTrack->GetStatus()&AliESDtrack::kITSrefit)) cut=kTRUE; // ITS refit
   Int_t clusterITS[200]; 
   if(esdTrack->GetITSclusters(clusterITS)<2) cut=kTRUE;  // min. nb. ITS clusters //3
  
   if(cut) return kTRUE;
   return kFALSE;
}

//______________________________________________________________________________________________________________________

void AliPerformancePtCalibMC::Analyse()
{
  
   // analyse charge/pt spectra in bins of theta and phi. Bins can be set by user
   
   AliPerfAnalyzeInvPt *ana = new  AliPerfAnalyzeInvPt("AliPerfAnalyzeInvPt","AliPerfAnalyzeInvPt");
  
   TH1::AddDirectory(kFALSE);
 
   ana->SetProjBinsTheta(fThetaBins,fNThetaBins);
   ana->SetProjBinsPhi(fPhiBins,fNPhiBins);
   ana->SetMakeFitOption(fFitGaus,fExclRange,fRange);
  
   TObjArray *aFolderObj = new TObjArray;
   if(fAnaMC){
      Printf("AliPerformancePtCalibMC: analysing MC!");
      ana->StartAnalysis(fHistInvPtThetaMC,fHistInvPtPhiMC, aFolderObj);
   }
   else {
      Printf("AliPerformancePtCalibMC: analysing data!");
      ana->StartAnalysis(fHistInvPtTheta,fHistInvPtPhi, aFolderObj);
   }
  
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
  
	 fHistInvPtTheta->Add(entry->fHistInvPtTheta);
	 fHistInvPtPhi->Add(entry-> fHistInvPtPhi);
	 fHistPtTheta->Add(entry->fHistPtTheta);
	 fHistPtPhi->Add(entry->fHistPtPhi);

	 fHistInvPtThetaMC->Add(entry->fHistInvPtThetaMC);
	 fHistInvPtPhiMC->Add(entry->fHistInvPtPhiMC);
	 fHistPtThetaMC->Add(entry->fHistPtThetaMC);
	 fHistPtPhiMC->Add(entry->fHistPtPhiMC);
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
void AliPerformancePtCalibMC::SetProjBinsPhi(const Double_t *phiBinArray,const Int_t nphBins){

   // set phi bins for Analyse()
   //set phi bins as array and set number of this array which is equal to number of bins analysed
   //the last analysed bin will always be the projection from first to last bin in the array
   if(nphBins){
      fNPhiBins = nphBins;
  
      for(Int_t k = 0;k<fNPhiBins;k++){
	 fPhiBins[k] = phiBinArray[k];
      }
      Printf("AliPerformancePtCalibMC: number of bins in phi set to %i",fNPhiBins);
   }
   else  Printf("Warning AliPerformancePtCalibMC::SetProjBinsPhi:  number of bins in phi NOT set!!! Default values are taken.");
}
//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalibMC::SetProjBinsTheta(const Double_t *thetaBinArray, const Int_t nthBins){
   // set theta bins for Analyse()
   //set theta bins as array and set number of this array which is equal to number of bins analysed
   //the last analysed bin will always be the projection from first to last bin in the array
   if(nthBins){
      fNThetaBins = nthBins;
      for(Int_t k = 0;k<fNThetaBins;k++){
	 fThetaBins[k] = thetaBinArray[k];
      }
      Printf("AliPerformancePtCalibMC: number of bins in theta set to %i",fNThetaBins);
   }
   else  Printf("Warning AliPerformancePtCalibMC::SetProjBinsTheta:  number of bins in theta NOT set!!! Default values are taken.");
}
//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalibMC::SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR ){

   //set the fit options:
   //for usage of gaussian function instead of polynomial (default) set setGausFit=kTRUE
   //set the range of rejection of points around 0 via exclusionR
   //set the fit range around 0 with fitR
   
   fFitGaus = setGausFit;
   fExclRange  = exclusionR;
   fRange = fitR;
  
   if(fFitGaus) Printf("AliPerformancePtCalibMC:set MakeGausFit with fit range %2.3f and exclusion range in 1/pt: %2.3f",fRange,fExclRange);
   else  Printf("AliPerformancePtCalibMC: set standard polynomial fit with fit range %2.3f and exclusion range in 1/pt: %2.3f",fRange,fExclRange);
 
}
//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalibMC::SetESDcutValues(const Double_t * esdCutValues){
   // set ESD cut values as an array of size 6
    
   fMinPt                = esdCutValues[0]; 
   fMaxPt                = esdCutValues[1];
   fMinNClustersTPC      = esdCutValues[2];
   fMaxChi2PerClusterTPC = esdCutValues[3];
   fMaxDCAtoVertexXY     = esdCutValues[4];
   fMaxDCAtoVertexZ      = esdCutValues[5];
}
