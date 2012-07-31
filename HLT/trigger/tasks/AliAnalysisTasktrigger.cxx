//-*- Mode: C++ -*-

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliTracker.h" 
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliKineTrackCuts.h"
#include "AliMCParticle.h"

#include "AliAnalysisTasktrigger.h"

using namespace std;

// Study trigger efficiencies for high-pt trigger
// Author: Jochen Thaeder <jochen@thaeder.de>

ClassImp(AliAnalysisTasktrigger)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisTasktrigger::AliAnalysisTasktrigger(const char *name) :
  AliAnalysisTaskSE(name), 
  fRandom(NULL),
  fMC(NULL), fESD(NULL), fESDHLT(NULL), 
  fESDTrackCuts(NULL), fESDHLTTrackCuts(NULL), fMCTrackCuts(NULL),
  fIsSelected(kFALSE), fIsSelectedHLT(kFALSE), fIsSelectedMC(kFALSE),  
  fIsSelectedTask(kFALSE),
  fOutputContainer(NULL),
  fPtCount(NULL), fMultiplicity(NULL),
  fTrigger(NULL) {
  // Constructor
      
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1, TObjArray::Class());
}

//________________________________________________________________________
const Int_t     AliAnalysisTasktrigger::fgkNSettings       = 7;
//________________________________________________________________________
const Double_t  AliAnalysisTasktrigger::fgkTriggerPt[]     = {1.0, 2.0, 2.5, 3.0, 5.0, 7.0, 10.0};


//________________________________________________________________________
const Int_t     AliAnalysisTasktrigger::fgkNTrigger        = 10;

//________________________________________________________________________
const Char_t   *AliAnalysisTasktrigger::fgkTrigger[]       = { 
  "p_{t}> 1.0",              "p_{t}> 2.0",              "p_{t}> 2.5",
  "p_{t}> 3.0",              "p_{t}> 5.0",              "p_{t}> 7.0",
  "p_{t}> 10.0",
  "S1",                      "S2",                      "S3",
};

//________________________________________________________________________
const Int_t     AliAnalysisTasktrigger::fgkNSelectionCuts  = 5;

//________________________________________________________________________
const Char_t   *AliAnalysisTasktrigger::fgkSelectionCuts[] = { 
  "All Events", 
  "AliPhysSel", 
  "AliPhysSel - PrimVertex",
  "AliPhysSel - PrimVertex - Track (OFF)",
  "AliPhysSel - PrimVertex - Track (HLT)",
};

//________________________________________________________________________
AliAnalysisTasktrigger::~AliAnalysisTasktrigger() {
  // Destructor

  if (fRandom)
    delete fRandom;
  fRandom = NULL;

  if (fESDTrackCuts)
    delete fESDTrackCuts;
  fESDTrackCuts = NULL;

  if (fESDHLTTrackCuts)
    delete fESDHLTTrackCuts;
  fESDHLTTrackCuts = NULL;

  if (fMCTrackCuts)
    delete fMCTrackCuts;
  fMCTrackCuts = NULL;

  if (fPtCount)
    delete[] fPtCount;
  fPtCount = NULL;

  if (fMultiplicity)
    delete[] fMultiplicity;
  fMultiplicity = NULL;

  if (fTrigger)
    delete[] fTrigger;
  fTrigger = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisTasktrigger::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  fOutputContainer = new TObjArray(1);
  fOutputContainer->SetName(GetName()) ;
  fOutputContainer->SetOwner(kTRUE) ;

  SetupTrigHistograms();
  SetupPtHistograms();
  SetupMultHistograms();

  SetupESDTrackCuts();

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  fTrigger       = new Bool_t[kNModes*fgkNTrigger];
  fPtCount       = new Int_t[kNModes*fgkNTrigger];

  fMultiplicity  = new Int_t[fgkNSelectionCuts];
  
  fRandom = new TRandom3(0);
}

//________________________________________________________________________
void AliAnalysisTasktrigger::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  // -- Setup Event
  // ----------------
  if ( !SetupEvent() ) 
    return;
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  for (Int_t mode = 0; mode < kNModes; ++mode) {

    if (mode == kMC && !fMC) continue;

    // -- Fill Cut Studies
    // ---------------------
    FillCutStudies(mode);

    // -- Fill Counter loop
    // ----------------------
    FillCounters(mode);
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  // -- Evaluate Trigger
  // ---------------------
  EvaluateTrigger();

  // -- Fill Trigger Histograms
  // ----------------------------
  FillTriggerHistograms();

  // -- Fill Trigger Studies
  // -------------------------
  FillTriggerStudies();
  FillTriggerStudiesMC();
    
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  
  // -- Post output data.
  PostData(1, fOutputContainer);
}      


/*
 * ---------------------------------------------------------------------------------
 *                            Trigger Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisTasktrigger::EvaluateTrigger() {
  // Evaluate Counters and Trigger

  for (Int_t mode = 0; mode < kNModes; ++mode) {
    Int_t offset = mode*fgkNTrigger;

    // -- Single Particle Triggers
    for ( Int_t idx = 0; idx < fgkNSettings; ++idx ) {
     
      if ( fPtCount[offset+idx] >= 1 )
	fTrigger[offset+idx] = kTRUE;
    }

    // -- Scenarios
    // ------------------
    Bool_t s1 = kFALSE;
    Bool_t s2 = kFALSE;
    Bool_t s3 = kFALSE;

    if (!fRandom)
      printf(" ERROR NO RANDOM !!!! \n");

    // S1 : Reduce to 300 Hz --  [ 150 Hz min bias ]
    // ---------------------------------------------
  
    //  > 1 GeV | 800 Hz * 0.1         ~ 80 Hz
    if ( fPtCount[offset+0] >= 1 )
      if ( fRandom->Rndm() < 0.1 )
	s1 = kTRUE;
    //  > 2 GeV | 230 Hz * 0.22        ~ 50 Hz
    if ( fPtCount[offset+1] >= 1 )
      if ( fRandom->Rndm() < 0.22 )
	s1 = kTRUE;
    //  > 3 GeV | 60 Hz * 0.5          ~ 30 Hz
    if ( fPtCount[offset+3] >= 1 )
      if ( fRandom->Rndm() < 0.5 )
	s1 = kTRUE;
    //  > 5 GeV |                      ~ 8  Hz
    if ( fPtCount[offset+4] >= 1 )
      s1 = kTRUE;


    // S2 : Reduce to 400 Hz --  [ 200 Hz min bias ]
    // ---------------------------------------------
  
    //  > 1 GeV | 800 Hz * 0.14        ~110 Hz
    if ( fPtCount[offset+0] >= 1 )
      if ( fRandom->Rndm() < 0.14 )
	s2 = kTRUE;
    //  > 2 GeV | 230 Hz * 0.39        ~ 90 Hz
    if ( fPtCount[offset+1] >= 1 )
      if ( fRandom->Rndm() < 0.39 )
	s2 = kTRUE;
    //  > 3 GeV | 60 Hz * 0.5          ~ 30 Hz
    if ( fPtCount[offset+3] >= 1 )
      if ( fRandom->Rndm() < 0.5 )
	s2 = kTRUE;
    //  > 5 GeV |                      ~ 8  Hz
    if ( fPtCount[offset+4] >= 1 )
      s2 = kTRUE;


    // S3 : Reduce to 500 Hz --  [ 200 Hz min bias ]
    // ---------------------------------------------
 
    //  > 1 GeV | 800 Hz * 0.2         ~160 Hz
    if ( fPtCount[offset+0] >= 1 )
      if ( fRandom->Rndm() < 0.2 )
	s3 = kTRUE;
    //  > 2 GeV | 230 Hz * 0.61        ~140 Hz
    if ( fPtCount[offset+1] >= 1 )
      if ( fRandom->Rndm() < 0.61 )
	s3 = kTRUE;
    //  > 3 GeV |                      ~ 60 Hz
    if ( fPtCount[offset+3] >= 1 )
      s3 = kTRUE;

    fTrigger[offset+fgkNSettings]   = s1;
    fTrigger[offset+fgkNSettings+1] = s2;
    fTrigger[offset+fgkNSettings+2] = s3;
  } 

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Setup Cuts Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisTasktrigger::SetupESDTrackCuts() {
  // Setup ESD cuts

  // -- Pure Off-line track cuts
  //  fESDTrackCuts    = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(kTRUE);

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  // -- HLT adopted track cuts
  fESDHLTTrackCuts = new AliESDtrackCuts;

  // -- turn off criteria
  fESDHLTTrackCuts->SetDCAToVertex2D(kFALSE);
  fESDHLTTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDHLTTrackCuts->SetRequireTPCRefit(kFALSE);
  fESDHLTTrackCuts->SetRequireITSRefit(kFALSE);

  // -- CleanSample
  fESDHLTTrackCuts->SetMaxDCAToVertexXY(3.0);
  fESDHLTTrackCuts->SetMaxDCAToVertexZ(10.0);
  fESDHLTTrackCuts->SetEtaRange(-1.0,1.0);
  fESDHLTTrackCuts->SetMinNClustersTPC(60);

  // -- Pre Trigger Bias
  fESDHLTTrackCuts->SetPtRange(0.3);         

  fESDTrackCuts = fESDHLTTrackCuts;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  fMCTrackCuts = new AliKineTrackCuts;
  fMCTrackCuts->SetEtaRange(-1.0,1.0);
  fMCTrackCuts->SetPtRange(0.3);

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                              Fill Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisTasktrigger::FillCutStudies( Int_t mode ) {
  // Fill histograms for cut studies

  for ( Int_t idx = 0; idx < fgkNSelectionCuts; idx++ )
    fMultiplicity[idx] = 0;

  Char_t *title = "";
  Bool_t hasVertex = kFALSE;
  AliESDEvent *esd = NULL;

  if (mode == kOFF) {
    title = "OFF";
    hasVertex = fIsSelected;
    esd = fESD;
  }
  else if (mode == kHLT) {
    title = "HLT";
    hasVertex = fIsSelectedHLT;
    esd = fESDHLT;
  }
  else if (mode == kMC) {
    title = "MC";
    hasVertex = fIsSelectedMC;
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
 
  if (mode == kOFF || mode == kHLT) {

    for (Int_t idx = 0; idx < esd->GetNumberOfTracks(); idx++) {
      
      AliESDtrack* track = esd->GetTrack(idx);
      if (!track) continue;
      
      // -- "All Events"
      // --------------------------------------------
      Int_t cut = 0; 
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(track->Pt());
      ++(fMultiplicity[cut]);
      
      // -- "AliPhysSel"
      // -------------------------------------------- 
      if (!fIsSelectedTask) 
	continue;
      
      cut = 1;
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(track->Pt());
      ++(fMultiplicity[cut]);
      
      // -- "AliPhysSel - PrimVertex"
      // --------------------------------------------
      if (!hasVertex) 
	continue;
      
      cut = 2;
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(track->Pt());
      ++(fMultiplicity[cut]);
      
      // -- "AliPhysSel - PrimVertex - Track (OFF)"
      // --------------------------------------------
      if (fESDTrackCuts->AcceptTrack(track)) {
	cut = 3;
	(static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(track->Pt());
	++(fMultiplicity[cut]);
      }

      // -- "AliPhysSel - PrimVertex - Track (HLT)"
      // --------------------------------------------
      if (fESDHLTTrackCuts->AcceptTrack(track)) {
	cut = 4;
	(static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(track->Pt());
	++(fMultiplicity[cut]);
      }

    } // for (Int_t idx = 0; idx < esd->GetNumberOfTracks(); idx++) {
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  else {

    AliStack* stack = fMC->Stack();  
    if (!stack) {
      printf("Error : No Stack. \n");
      return;
    }
    
    for (Int_t iterStack = 0; iterStack <  stack->GetNtrack(); iterStack++) {

      TParticle *particle = GetChargedPhysicalPrimary(stack, iterStack);
      if (!particle)
	continue;

      // -- "All Events"
      Int_t cut = 0; 
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(particle->Pt());
      ++(fMultiplicity[cut]);
      
      // -- "AliPhysSel" 
      // -----------------
      if (!fIsSelectedTask) 
	continue;
      
      cut = 1;
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(particle->Pt());
      ++(fMultiplicity[cut]);
      
      // -- "AliPhysSel - PrimVertex"
      // --------------------------------------------
      if (!hasVertex) 
	continue;

      cut = 2;
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(particle->Pt());
      ++(fMultiplicity[cut]);

      // -- "AliPhysSel - PrimVertex - Track (MC)"
      // --------------------------------------------
      if (!fMCTrackCuts->IsSelected(particle))
	continue;

      cut = 3;
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(particle->Pt());
      ++(fMultiplicity[cut]);
      
      // -- "AliPhysSel - PrimVertex - Track (MC) - Findable"
      // ------------------------------------------------------
      if (!IsFindableMC(iterStack, 60.))
	continue;
      
      cut = 4;
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sPt_%d", title, cut))))->Fill(particle->Pt());
      ++(fMultiplicity[cut]);
    
    } // for (Int_t iterStack = 0; iterStack < stack->GetNtrack(); iterStack++) {
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  
  // -- Fill Multiplicity
  for ( Int_t idx = 0; idx < fgkNSelectionCuts; idx++ )
    (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHist%sMult_%d", title, idx))))->Fill(fMultiplicity[idx]);

  return;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::FillCounters( Int_t mode ) {
  // Fill counters for trigger
  
  Int_t offset   = mode*fgkNTrigger ;
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  if (mode == kOFF || mode == kHLT) {

    AliESDtrackCuts *cuts = fESDTrackCuts;
    AliESDEvent *esd      = fESD;

    if (mode == kHLT) {
      esd   = fESDHLT;
      cuts  = fESDHLTTrackCuts;
    }
    
    for (Int_t idx = 0; idx < esd->GetNumberOfTracks(); idx++) {
      AliESDtrack* track = esd->GetTrack(idx);
      if (!track) continue;
      
      // -- ESD cuts
      if ( !cuts->AcceptTrack(track) ) continue;

      Int_t   nClustersTPC = track->GetTPCclusters(0);

      //      Float_t ptTr        = TMath::Abs(track->GetSigned1Pt()) - 3.5 * TMath::Sqrt( track->GetSigma1Pt2());

      // -- Different Pt settings
      for ( Int_t set=0; set < fgkNSettings ; set++ ) {

	if      ( fgkTriggerPt[set] ==  3. ) { if ( nClustersTPC < 100 ) continue; }
	else if ( fgkTriggerPt[set] ==  5. ) { if ( nClustersTPC < 110 ) continue; }
	else if ( fgkTriggerPt[set] ==  7. ) { if ( nClustersTPC < 120 ) continue; }
	else if ( fgkTriggerPt[set] == 10. ) { if ( nClustersTPC < 140 ) continue; }
	else                                 { if ( nClustersTPC <  80 ) continue; }
	
	//	if ( ptTr < 1/fgkTriggerPt[set] )
	//	  ++(fPtCount[offset+set]);

       	if ( track->Pt() > fgkTriggerPt[set] )
	  ++(fPtCount[offset+set]);
      }
      
    } // for (Int_t idx = 0; idx < esd->GetNumberOfTracks(); idx++) {
  } // if (mode == kOFF || mode == kHLT) {
  else {

    AliStack* stack = fMC->Stack();
    if (!stack) {
      printf("Error : No Stack. \n");
      return;
    }
    
    for (Int_t iterStack = 0; iterStack < stack->GetNtrack(); iterStack++) {

      TParticle *particle = GetChargedPhysicalPrimary(stack, iterStack);
      if (!particle)
	continue;
      
      // -- MC cuts
      if ( !fMCTrackCuts->IsSelected(particle) ) continue;

      // -- Different Pt settings
      for ( Int_t set=0; set < fgkNSettings ; set++ ) {

	if      ( fgkTriggerPt[set] ==  3. ) { if ( !IsFindableMC(iterStack, 100.) ) continue; }
	else if ( fgkTriggerPt[set] ==  5. ) { if ( !IsFindableMC(iterStack, 110.) ) continue; }
	else if ( fgkTriggerPt[set] ==  7. ) { if ( !IsFindableMC(iterStack, 120.) ) continue; }
	else if ( fgkTriggerPt[set] == 10. ) { if ( !IsFindableMC(iterStack, 140.) ) continue; }
	else                                 { if ( !IsFindableMC(iterStack,  80.) ) continue; }

	if ( particle->Pt() > fgkTriggerPt[set] )
	  ++(fPtCount[offset+set]);
      }
    } // for (Int_t iterStack = 0; iterStack < stack->GetNtrack(); iterStack++) {
  }
    
  return;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::FillTriggerHistograms() {
  // Fill histograms for trigger
  
  Int_t nTracks = fESD->GetNumberOfTracks();
  Int_t nTracksH = fESDHLT->GetNumberOfTracks();

  Int_t nTracksM = 0;

  if (fMC)
    nTracksM = fMC->Stack()->GetNtrack();

  static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFTriggered"))->Fill(0);
  static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFTriggered"))->Fill(fgkNTrigger+1, nTracks);
  static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTTriggered"))->Fill(0);
  static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTTriggered"))->Fill(fgkNTrigger+1, nTracksH);
  static_cast<TH1F*>(fOutputContainer->FindObject("fNMCTriggered"))->Fill(0);
  static_cast<TH1F*>(fOutputContainer->FindObject("fNMCTriggered"))->Fill(fgkNTrigger+1, nTracksM);

  for ( Int_t idx = 0; idx < fgkNTrigger; idx++ ) {
    if (fTrigger[(kOFF*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFTriggered"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFTriggered"))->Fill(fgkNTrigger+idx+2,nTracks);
    }
    if (fTrigger[(kHLT*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTTriggered"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTTriggered"))->Fill(fgkNTrigger+idx+2,nTracksH);
    }
    if (fTrigger[(kMC*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNMCTriggered"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNMCTriggered"))->Fill(fgkNTrigger+idx+2,nTracksM);
    }

    if (fTrigger[(kHLT*fgkNTrigger)+idx] && !fTrigger[(kOFF*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTFakeToOFF"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTFakeToOFF"))->Fill(fgkNTrigger+idx+2,nTracks);
    }
    if (fTrigger[(kHLT*fgkNTrigger)+idx] && !fTrigger[(kMC*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTFakeToMC"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTFakeToMC"))->Fill(fgkNTrigger+idx+2,nTracks);
    }
    if (fTrigger[(kOFF*fgkNTrigger)+idx] && !fTrigger[(kMC*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFFakeToMC"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFFakeToMC"))->Fill(fgkNTrigger+idx+2,nTracks);
    }


    if (!fTrigger[(kHLT*fgkNTrigger)+idx] && fTrigger[(kOFF*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTMissToOFF"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTMissToOFF"))->Fill(fgkNTrigger+idx+2,nTracks);
    }
    if (!fTrigger[(kHLT*fgkNTrigger)+idx] && fTrigger[(kMC*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTMissToMC"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTMissToMC"))->Fill(fgkNTrigger+idx+2,nTracks);
    }
    if (!fTrigger[(kOFF*fgkNTrigger)+idx] && fTrigger[(kMC*fgkNTrigger)+idx]) {
      static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFMissToMC"))->Fill(idx+1);
      static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFMissToMC"))->Fill(fgkNTrigger+idx+2,nTracks);
    }

  }

  if ( fIsSelected ) {
    static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFTriggeredSel"))->Fill(0);
    static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFTriggeredSel"))->Fill(fgkNTrigger+1, nTracks);
    static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTTriggeredSel"))->Fill(0);
    static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTTriggeredSel"))->Fill(fgkNTrigger+1, nTracksH);
    static_cast<TH1F*>(fOutputContainer->FindObject("fNMCTriggeredSel"))->Fill(0);
    static_cast<TH1F*>(fOutputContainer->FindObject("fNMCTriggeredSel"))->Fill(fgkNTrigger+1, nTracksM);
    
    for ( Int_t idx = 0; idx < fgkNTrigger; idx++ ) {
      if (fTrigger[(kOFF*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFTriggeredSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFTriggeredSel"))->Fill(fgkNTrigger+idx+2,nTracks);
      }
      if (fTrigger[(kHLT*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTTriggeredSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTTriggeredSel"))->Fill(fgkNTrigger+idx+2,nTracksH);
      }
      if (fTrigger[(kMC*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNMCTriggeredSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNMCTriggeredSel"))->Fill(fgkNTrigger+idx+2,nTracksM);
      }

      if (fTrigger[(kHLT*fgkNTrigger)+idx] && !fTrigger[(kOFF*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTFakeToOFFSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTFakeToOFFSel"))->Fill(fgkNTrigger+idx+2,nTracks);
      }
      if (fTrigger[(kHLT*fgkNTrigger)+idx] && !fTrigger[(kMC*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTFakeToMCSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTFakeToMCSel"))->Fill(fgkNTrigger+idx+2,nTracks);
      }
      if (fTrigger[(kOFF*fgkNTrigger)+idx] && !fTrigger[(kMC*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFFakeToMCSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFFakeToMCSel"))->Fill(fgkNTrigger+idx+2,nTracks);
      }

      if (!fTrigger[(kHLT*fgkNTrigger)+idx] && fTrigger[(kOFF*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTMissToOFFSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTMissToOFFSel"))->Fill(fgkNTrigger+idx+2,nTracks);
      }
      if (!fTrigger[(kHLT*fgkNTrigger)+idx] && fTrigger[(kMC*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTMissToMCSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNHLTMissToMCSel"))->Fill(fgkNTrigger+idx+2,nTracks);
      }
      if (!fTrigger[(kOFF*fgkNTrigger)+idx] && fTrigger[(kMC*fgkNTrigger)+idx]) {
	static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFMissToMCSel"))->Fill(idx+1);
	static_cast<TH1F*>(fOutputContainer->FindObject("fNOFFMissToMCSel"))->Fill(fgkNTrigger+idx+2,nTracks);
      }
    }
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::FillTriggerStudies() {
  // Fill histograms for trigger studies

  // -- "AliPhysSel - PrimVertex"
  // --------------------------------------------
  if (!fIsSelected) 
    return;

  Int_t multiplicity = 0;
  
  for (Int_t idx = 0; idx < fESD->GetNumberOfTracks(); idx++) {
    AliESDtrack* track = fESD->GetTrack(idx);
    if (!track) continue;

    // -- "AliPhysSel - PrimVertex - Track (OFF)"
    // --------------------------------------------
    if (!fESDTrackCuts->AcceptTrack(track))
      continue;

    // -- Fill Triggered Pt Spectra
    for ( Int_t trg= 0; trg < fgkNTrigger; trg++ ) {
      
      // -- ESD
      if ( fTrigger[(kOFF*fgkNTrigger)+trg] )
	(static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistOFF_PtTriggered_OFF_%d",trg))))->Fill(track->Pt());
      
      // -- HLT
      if ( fTrigger[(kHLT*fgkNTrigger)+trg] )
	(static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistOFF_PtTriggered_HLT_%d",trg))))->Fill(track->Pt());

      // -- MC
      if ( fTrigger[(kMC*fgkNTrigger)+trg] )
	(static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistOFF_PtTriggered_MC_%d",trg))))->Fill(track->Pt());
      
    } // for ( Int_t trg= 0; trg < fgkNTrigger; trg++ ) {
    
    ++multiplicity;
    
  } // for (Int_t idx = 0; idx < esd->GetNumberOfTracks(); idx++) {

  // -- Fill Triggered Mult Spectra
  for ( Int_t trg= 0; trg < fgkNTrigger; trg++ ) {
    
    // -- ESD
    if ( fTrigger[(kOFF*fgkNTrigger)+trg] )
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistOFF_MultTriggered_OFF_%d",trg))))->Fill(multiplicity);
    
    // -- HLT
    if ( fTrigger[(kHLT*fgkNTrigger)+trg] )
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistOFF_MultTriggered_HLT_%d",trg))))->Fill(multiplicity);

    // -- MC
    if ( fTrigger[(kMC*fgkNTrigger)+trg] )
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistOFF_MultTriggered_MC_%d",trg))))->Fill(multiplicity);

  } // for ( Int_t trg= 0; trg < fgkNTrigger; trg++ ) {

  return;
}


//________________________________________________________________________
void AliAnalysisTasktrigger::FillTriggerStudiesMC() {
  // Fill histograms for trigger studies

  if (!fMC)
    return;

  // -- "AliPhysSel - PrimVertex"
  // --------------------------------------------
  if (!fIsSelectedMC) 
    return;

  Int_t multiplicity = 0;

  AliStack* stack = fMC->Stack();
  if (!stack) {
    printf("Error : No Stack. \n");
    return;
  }
    
  for (Int_t iterStack = 0; iterStack < stack->GetNtrack(); iterStack++) {
    
    TParticle *particle = GetChargedPhysicalPrimary(stack, iterStack);
    if (!particle)
      continue;

    if ( ! fMCTrackCuts->IsSelected(particle) )
      continue;

    if ( ! IsFindableMC(iterStack, 60.) )
      continue;

    for ( Int_t trg= 0; trg < fgkNTrigger; trg++ ) {
      
      // -- ESD
      if ( fTrigger[(kOFF*fgkNTrigger)+trg] )
	(static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistMC_PtTriggered_OFF_%d",trg))))->Fill(particle->Pt());
      
      // -- HLT
      if ( fTrigger[(kHLT*fgkNTrigger)+trg] )
	(static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistMC_PtTriggered_HLT_%d",trg))))->Fill(particle->Pt());
      
      // -- MC
      if ( fTrigger[(kMC*fgkNTrigger)+trg] )
	(static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistMC_PtTriggered_MC_%d",trg))))->Fill(particle->Pt());
      
    } // for ( Int_t trg= 0; trg < fgkNTrigger; trg++ ) {
    
    ++multiplicity;
    
  } // for (Int_t iterStack = 0; iterStack < stack->GetNtrack(); iterStack++) {

  // -- Fill Triggered Mult Spectra
  for ( Int_t trg= 0; trg < fgkNTrigger; trg++ ) {
    
    // -- ESD
    if ( fTrigger[(kOFF*fgkNTrigger)+trg] )
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistMC_MultTriggered_OFF_%d",trg))))->Fill(multiplicity);
    
    // -- HLT
    if ( fTrigger[(kHLT*fgkNTrigger)+trg] )
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistMC_MultTriggered_HLT_%d",trg))))->Fill(multiplicity);

    // -- MC
    if ( fTrigger[(kMC*fgkNTrigger)+trg] )
      (static_cast<TH1F*>(fOutputContainer->FindObject(Form("fHistMC_MultTriggered_MC_%d",trg))))->Fill(multiplicity);

  } // for ( Int_t trg= 0; trg < fgkNTrigger; trg++ ) {

  return;
}


/*
 * ---------------------------------------------------------------------------------
 *                            Setup Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Bool_t AliAnalysisTasktrigger::SetupEvent() {
  // Setup Reading of event

  // -- Clear Counters and Triggers
  // --------------------------------
  for ( Int_t idx = 0; idx < kNModes*fgkNTrigger; idx++ ) {
    fTrigger[idx] = kFALSE;
    fPtCount[idx] = 0;
  }  
  
  fIsSelected = fIsSelectedHLT = fIsSelectedMC = fIsSelectedTask = kFALSE;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  // -- ESD Event Handler
  // ----------------------

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    printf("ERROR: Could not get ESDInputHandler");
    return kFALSE;
  } 

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
 
  // -- MC Event Handler  - MC Event
  // ---------------------------------
 
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  // -- Events
  // ------------
  fESD = esdH->GetEvent();
  if (!fESD) {
    printf("ERROR: fESD not available \n");
    return kFALSE;
  }

  fESDHLT = esdH->GetHLTEvent();
  if (!fESDHLT) {
    printf("ERROR: fESDHLT not available \n");
    return kFALSE;
  }

  fMC = MCEvent();
  if ( mcH && !fMC) {
    printf("ERROR: fMC not available \n");
    return kFALSE;
  }

  if ( !fESD->GetPrimaryVertexTracks() ){
    printf("ERROR: No Vertex \n");
    return kFALSE;
  }
  
  if ( !fESDHLT->GetPrimaryVertexTracks() ){
    printf("ERROR: No HLT Vertex \n");
    return kFALSE;
  }

  if ( fMC && !fMC->GetPrimaryVertex() ){
    printf("ERROR: No MC Vertex \n");
    return kFALSE;
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  // -- Physics Selection
  // ----------------------

  fIsSelectedTask = esdH->IsEventSelected() & AliVEvent::kMB;  
  fIsSelected     = fIsSelectedTask && (fESD->GetPrimaryVertexTracks())->GetStatus();
  fIsSelectedHLT  = fIsSelectedTask && (fESDHLT->GetPrimaryVertexTracks())->GetStatus();
  if (fMC)
    fIsSelectedMC   = fIsSelectedTask && (fMC->GetPrimaryVertex()); //->GetStatus();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::SetupTrigHistograms() {

  Int_t    n = 2*(fgkNTrigger+1);
  Double_t s = 0.;
  Double_t e = static_cast<Double_t>(2*(fgkNTrigger+1));

  AddTriggerHist(new TH1F("fNOFFTriggered",    "N events OFF triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNOFFTriggeredSel", "N events OFF triggered - OFF selected;;N Events", n, s, e));

  AddTriggerHist(new TH1F("fNHLTTriggered",    "N events HLT triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNHLTTriggeredSel", "N events HLT triggered - OFF selected;;N Events", n, s, e));

  AddTriggerHist(new TH1F("fNMCTriggered",     "N events MC triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNMCTriggeredSel",  "N events MC triggered - OFF selected;;N Events", n, s, e));

  // -- -- -- -- 
  
  AddTriggerHist(new TH1F("fNHLTFakeToOFF",    "N events HLT fake (to OFF) triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNHLTFakeToOFFSel", "N events HLT fake (to OFF) triggered - OFF selected;;N Events", n, s, e));

  AddTriggerHist(new TH1F("fNHLTFakeToMC",     "N events HLT fake (to MC) triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNHLTFakeToMCSel",  "N events HLT fake (to MC) triggered - OFF selected;;N Events", n, s, e));

  AddTriggerHist(new TH1F("fNOFFFakeToMC",     "N events OFF fake (to MC) triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNOFFFakeToMCSel",  "N events OFF fake (to MC) triggered - OFF selected;;N Events", n, s, e));

  // -- -- -- -- 

  AddTriggerHist(new TH1F("fNHLTMissToOFF",    "N events HLT miss (to OFF) triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNHLTMissToOFFSel", "N events HLT miss (to OFF) triggered - MC selected;;N Events", n, s, e));

  AddTriggerHist(new TH1F("fNHLTMissToMC",     "N events HLT miss (to MC) triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNHLTMissToMCSel",  "N events HLT miss (to MC) triggered - MC selected;;N Events", n, s, e));

  AddTriggerHist(new TH1F("fNOFFMissToMC",     "N events OFF miss (to MC) triggered;;N Events", n, s, e));
  AddTriggerHist(new TH1F("fNOFFMissToMCSel",  "N events OFF miss (to MC) triggered - MC selected;;N Events", n, s, e));

  return;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::SetupPtHistograms() {

  Int_t    n = 50;
  Double_t s = 0.1;
  Double_t e = 50.;
  
  // -- Create LogPtBinning
  // ------------------------
  Double_t logMin = TMath::Log10(s);
  Double_t logMax = TMath::Log10(e);
  Double_t binwidth = (logMax-logMin)/n;
  
  Double_t *logBinning = new Double_t[n+1];

  logBinning[0] = s;
  for (Int_t ii = 1; ii <= n; ii++)
    logBinning[ii] = s + TMath::Power(10, logMin + ii*binwidth);

  // -- Create Histograms
  // ----------------------
  for ( Int_t idx = 0; idx < fgkNSelectionCuts; idx++ ) {
    AddPtHist(new TH1F(Form("fHistOFFPt_%d",idx), 
		       Form("OFF P_{T} distribution - %s", fgkSelectionCuts[idx]), n, logBinning));
    AddPtHist(new TH1F(Form("fHistHLTPt_%d",idx), 
		       Form("HLT P_{T} distribution - %s", fgkSelectionCuts[idx]), n, logBinning));
    AddPtHist(new TH1F(Form("fHistMCPt_%d",idx), 
		       Form("MC P_{T} distribution - %s", fgkSelectionCuts[idx]), n, logBinning));
  }
  
  for ( Int_t idx = 0; idx < fgkNTrigger; idx++ ) {
    AddPtHist(new TH1F(Form("fHistOFF_PtTriggered_OFF_%d",idx), 
		       Form("OFF Triggered - OFF P_{T} distribution - %s", fgkTrigger[idx]), n, logBinning));
    AddPtHist(new TH1F(Form("fHistOFF_PtTriggered_HLT_%d",idx), 
		       Form("HLT Triggered - OFF P_{T} distribution - %s", fgkTrigger[idx]), n, logBinning));
    AddPtHist(new TH1F(Form("fHistOFF_PtTriggered_MC_%d",idx), 
		       Form("MC Triggered - OFF P_{T} distribution - %s", fgkTrigger[idx]), n, logBinning));
    AddPtHist(new TH1F(Form("fHistMC_PtTriggered_OFF_%d",idx), 
		       Form("OFF Triggered - MC P_{T} distribution - %s", fgkTrigger[idx]), n, logBinning));
    AddPtHist(new TH1F(Form("fHistMC_PtTriggered_HLT_%d",idx), 
		       Form("HLT Triggered - MC P_{T} distribution - %s", fgkTrigger[idx]), n, logBinning));
    AddPtHist(new TH1F(Form("fHistMC_PtTriggered_MC_%d",idx), 
		       Form("MC Triggered - MC P_{T} distribution - %s", fgkTrigger[idx]), n, logBinning));
  }
  
  delete[] logBinning;
  logBinning = NULL;

  return;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::SetupMultHistograms() {

  Int_t    n = 51;
  Double_t s = 0.;
  Double_t e = 50.;

  for ( Int_t idx = 0; idx < fgkNSelectionCuts; idx++ ) {
    AddMultHist(new TH1F(Form("fHistOFFMult_%d",idx), 
			 Form("Multiplicity distribution - %s", fgkSelectionCuts[idx]), n, s ,e));
    AddMultHist(new TH1F(Form("fHistHLTMult_%d",idx), 
			 Form("HLT Multiplicity distribution - %s", fgkSelectionCuts[idx]), n, s ,e));
    AddMultHist(new TH1F(Form("fHistMCMult_%d",idx), 
			 Form("MC Multiplicity distribution - %s", fgkSelectionCuts[idx]), n, s ,e));
  }
  
  for ( Int_t idx = 0; idx < fgkNTrigger; idx++ ) {
    AddMultHist(new TH1F(Form("fHistOFF_MultTriggered_OFF_%d",idx), 
			 Form("OFF Triggered - OFF Multiplicty distribution - %s", fgkTrigger[idx]), n, s ,e));
    AddMultHist(new TH1F(Form("fHistOFF_MultTriggered_HLT_%d",idx), 
			 Form("HLT Triggered - OFF Multiplicty distribution - %s", fgkTrigger[idx]), n, s ,e));
    AddMultHist(new TH1F(Form("fHistOFF_MultTriggered_MC_%d",idx), 
			 Form("MC Triggered - OFF Multiplicty distribution - %s", fgkTrigger[idx]), n, s ,e));
    AddMultHist(new TH1F(Form("fHistMC_MultTriggered_OFF_%d",idx), 
			 Form("OFF Triggered - MC Multiplicty distribution - %s", fgkTrigger[idx]), n, s ,e));
    AddMultHist(new TH1F(Form("fHistMC_MultTriggered_HLT_%d",idx), 
			 Form("HLT Triggered - MC Multiplicty distribution - %s", fgkTrigger[idx]), n, s ,e));
    AddMultHist(new TH1F(Form("fHistMC_MultTriggered_MC_%d",idx), 
			 Form("MC Triggered - MC Multiplicty distribution - %s", fgkTrigger[idx]), n, s ,e));
  }

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Helper Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
TParticle* AliAnalysisTasktrigger::GetChargedPhysicalPrimary( AliStack* stack, Int_t idx ) {
  // return charged physical primary particle

  TParticle *particle = stack->Particle(idx);
  if (!particle)
    return NULL;

  // -- primary
  if ( !(stack->IsPhysicalPrimary(idx)) )
    return NULL;
      
  // -- charged only
  if (!particle->GetPDG())
    return NULL;
  
  if (particle->GetPDG()->Charge() == 0.0) 
    return NULL;

  return particle;
}

//________________________________________________________________________
Bool_t AliAnalysisTasktrigger::IsFindableMC(Int_t idx, Float_t length) {
  
  //  return kTRUE;

  // Ok, if track longer 60cm
  AliMCParticle *mcParticle = dynamic_cast<AliMCParticle*> (fMC->GetTrack(idx));
  if(!mcParticle)
    return kFALSE;

  Int_t counter; 
  Float_t tpcTrackLength = mcParticle->GetTPCTrackLength(AliTracker::GetBz(),0.05,counter,3.0); 

  if ( tpcTrackLength > length )
    return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::AddTriggerHist(TH1F* hist) {
  // add histogram to output
  fOutputContainer->Add(hist);
  return;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::AddPtHist(TH1F* hist) {
  // add histogram to output
  hist->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  hist->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fOutputContainer->Add(hist);
  return;
}

//________________________________________________________________________
void AliAnalysisTasktrigger::AddMultHist(TH1F* hist) {
  // add histogram to output
  hist->GetXaxis()->SetTitle("#ESD tracks");
  fOutputContainer->Add(hist);
  return;
}

