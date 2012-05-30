#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TRandom.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMixInputEventHandler.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"

#include "TH2D.h"                  
#include "AliPID.h"                
#include "AliPIDResponse.h"        
#include "AliPIDCombined.h"        

#include "AliAnalysisTaskEventMixingBF.h"
#include "AliBalanceEventMixing.h"


// Analysis task for the EventMixingBF code
// Authors: Panos.Christakoglou@nikhef.nl, m.weber@cern.ch

ClassImp(AliAnalysisTaskEventMixingBF)

//________________________________________________________________________
AliAnalysisTaskEventMixingBF::AliAnalysisTaskEventMixingBF(const char *name) 
: AliAnalysisTaskSE(name), 
  fBalance(0),
  fRunShuffling(kFALSE),
  fShuffledBalance(0),
  fList(0),
  fListEventMixingBF(0),
  fListEventMixingBFS(0),
  fHistListPIDQA(0),
  fHistEventStats(0),
  fHistCentStats(0),
  fHistTriggerStats(0),
  fHistTrackStats(0),
  fHistVx(0),
  fHistVy(0),
  fHistVz(0),
  fHistClus(0),
  fHistDCA(0),
  fHistChi2(0),
  fHistPt(0),
  fHistEta(0),
  fHistPhi(0),
  fHistPhiBefore(0),
  fHistPhiAfter(0),
  fHistV0M(0),
  fHistRefTracks(0),
  fHistdEdxVsPTPCbeforePID(NULL),
  fHistBetavsPTOFbeforePID(NULL), 
  fHistProbTPCvsPtbeforePID(NULL), 
  fHistProbTOFvsPtbeforePID(NULL), 
  fHistProbTPCTOFvsPtbeforePID(NULL),
  fHistNSigmaTPCvsPtbeforePID(NULL), 
  fHistNSigmaTOFvsPtbeforePID(NULL), 
  fHistdEdxVsPTPCafterPID(NULL),
  fHistBetavsPTOFafterPID(NULL), 
  fHistProbTPCvsPtafterPID(NULL), 
  fHistProbTOFvsPtafterPID(NULL), 
  fHistProbTPCTOFvsPtafterPID(NULL),
  fHistNSigmaTPCvsPtafterPID(NULL), 
  fHistNSigmaTOFvsPtafterPID(NULL),  
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fParticleOfInterest(kPion),
  fPidDetectorConfig(kTPCTOF),
  fUsePID(kFALSE),
  fUsePIDnSigma(kTRUE),
  fUsePIDPropabilities(kFALSE), 
  fPIDNSigma(3.),
  fMinAcceptedPIDProbability(0.8),
  fESDtrackCuts(0),
  fCentralityEstimator("V0M"),
  fUseCentrality(kFALSE),
  fCentralityPercentileMin(0.), 
  fCentralityPercentileMax(5.),
  fImpactParameterMin(0.),
  fImpactParameterMax(20.),
  fUseMultiplicity(kFALSE),
  fNumberOfAcceptedTracksMin(0),
  fNumberOfAcceptedTracksMax(10000),
  fHistNumberOfAcceptedTracks(0),
  fUseOfflineTrigger(kFALSE),
  fVxMax(0.3),
  fVyMax(0.3),
  fVzMax(10.),
  nAODtrackCutBit(128),
  fPtMin(0.3),
  fPtMax(1.5),
  fEtaMin(-0.8),
  fEtaMax(-0.8),
  fDCAxyCut(-1),
  fDCAzCut(-1),
  fTPCchi2Cut(-1),
  fNClustersTPCCut(-1),
  fAcceptanceParameterization(0),
  fDifferentialV2(0),
  fUseFlowAfterBurner(kFALSE),
  fExcludeResonancesInMC(kFALSE),
  fUseMCPdgCode(kFALSE),
  fPDGCodeToBeAnalyzed(-1),
  fMainEvent(0x0),
  fMixEvent(0x0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEventMixingBF::~AliAnalysisTaskEventMixingBF() {

  // delete fBalance; 
  // delete fShuffledBalance; 
  // delete fList;
  // delete fListEventMixingBF; 
  // delete fListEventMixingBFS;

  // delete fHistEventStats; 
  // delete fHistTrackStats; 
  // delete fHistVx; 
  // delete fHistVy; 
  // delete fHistVz; 

  // delete fHistClus;
  // delete fHistDCA;
  // delete fHistChi2;
  // delete fHistPt;
  // delete fHistEta;
  // delete fHistPhi;
  // delete fHistV0M;
}

//________________________________________________________________________
void AliAnalysisTaskEventMixingBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once
  if(!fBalance) {
    fBalance = new AliBalanceEventMixing();
    fBalance->SetAnalysisLevel("ESD");
    //fBalance->SetNumberOfBins(-1,16);
    fBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6);
  }
  if(fRunShuffling) {
    if(!fShuffledBalance) {
      fShuffledBalance = new AliBalanceEventMixing();
      fShuffledBalance->SetAnalysisLevel("ESD");
      //fShuffledBalance->SetNumberOfBins(-1,16);
      fShuffledBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6);
    }
  }

  //QA list
  fList = new TList();
  fList->SetName("listQA");
  fList->SetOwner();

  //Balance Function list
  fListEventMixingBF = new TList();
  fListEventMixingBF->SetName("listEventMixingBF");
  fListEventMixingBF->SetOwner();

  if(fRunShuffling) {
    fListEventMixingBFS = new TList();
    fListEventMixingBFS->SetName("listEventMixingBFShuffled");
    fListEventMixingBFS->SetOwner();
  }

  //PID QA list
  if(fUsePID) {
    fHistListPIDQA = new TList();
    fHistListPIDQA->SetName("listQAPID");
    fHistListPIDQA->SetOwner();
  }

  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  TString gCentName[9] = {"V0M","FMD","TRK","TKL","CL0","CL1","V0MvsFMD","TKLvsV0M","ZEMvsZDC"};
  fHistCentStats = new TH2F("fHistCentStats",
                             "Centrality statistics;;Cent percentile",
			    9,-0.5,8.5,220,-5,105);
  for(Int_t i = 1; i <= 9; i++)
    fHistCentStats->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
  fList->Add(fHistCentStats);

  fHistTriggerStats = new TH1F("fHistTriggerStats","Trigger statistics;TriggerBit;N_{events}",130,0,130);
  fList->Add(fHistTriggerStats);

  fHistTrackStats = new TH1F("fHistTrackStats","Event statistics;TrackFilterBit;N_{events}",130,0,130);
  fList->Add(fHistTrackStats);

  fHistNumberOfAcceptedTracks = new TH1F("fHistNumberOfAcceptedTracks",";N_{acc.};Entries",4001,-0.5,4000.5);
  fList->Add(fHistNumberOfAcceptedTracks);

  // Vertex distributions
  fHistVx = new TH1F("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVx);
  fHistVy = new TH1F("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVy);
  fHistVz = new TH1F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",100,-20.,20.);
  fList->Add(fHistVz);

  // QA histograms
  fHistClus = new TH2F("fHistClus","# Cluster (TPC vs. ITS)",10,0,10,200,0,200);
  fList->Add(fHistClus);
  fHistChi2 = new TH1F("fHistChi2","Chi2/NDF distribution",200,0,10);
  fList->Add(fHistChi2);
  fHistDCA  = new TH2F("fHistDCA","DCA (xy vs. z)",400,-5,5,400,-5,5); 
  fList->Add(fHistDCA);
  fHistPt   = new TH1F("fHistPt","p_{T} distribution",200,0,10);
  fList->Add(fHistPt);
  fHistEta  = new TH1F("fHistEta","#eta distribution",200,-2,2);
  fList->Add(fHistEta);
  fHistPhi  = new TH1F("fHistPhi","#phi distribution",200,-20,380);
  fList->Add(fHistPhi);
  fHistPhiBefore  = new TH1F("fHistPhiBefore","#phi distribution",200,0.,2*TMath::Pi());
  fList->Add(fHistPhiBefore);
  fHistPhiAfter  = new TH1F("fHistPhiAfter","#phi distribution",200,0.,2*TMath::Pi());
  fList->Add(fHistPhiAfter);
  fHistV0M  = new TH2F("fHistV0M","V0 Multiplicity C vs. A",500, 0, 20000, 500, 0, 20000);
  fList->Add(fHistV0M);
  TString gRefTrackName[6] = {"tracks","tracksPos","tracksNeg","tracksTPConly","clusITS0","clusITS1"};
  fHistRefTracks  = new TH2F("fHistRefTracks","Nr of Ref tracks/event vs. ref track estimator;;Nr of tracks",6, 0, 6, 400, 0, 20000);
  for(Int_t i = 1; i <= 6; i++)
    fHistRefTracks->GetXaxis()->SetBinLabel(i,gRefTrackName[i-1].Data());
  fList->Add(fHistRefTracks);

  // Balance function histograms
  // Initialize histograms if not done yet
  if(!fBalance->GetHistNp(0)){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    fBalance->InitHistograms();
  }

  if(fRunShuffling) {
    if(!fShuffledBalance->GetHistNp(0)) {
      AliWarning("Histograms (shuffling) not yet initialized! --> Will be done now");
      AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
      fShuffledBalance->InitHistograms();
    }
  }

  for(Int_t a = 0; a < ANALYSIS_TYPES; a++){
    fListEventMixingBF->Add(fBalance->GetHistNp(a));
    fListEventMixingBF->Add(fBalance->GetHistNn(a));
    fListEventMixingBF->Add(fBalance->GetHistNpn(a));
    fListEventMixingBF->Add(fBalance->GetHistNnn(a));
    fListEventMixingBF->Add(fBalance->GetHistNpp(a));
    fListEventMixingBF->Add(fBalance->GetHistNnp(a));

    if(fRunShuffling) {
      fListEventMixingBFS->Add(fShuffledBalance->GetHistNp(a));
      fListEventMixingBFS->Add(fShuffledBalance->GetHistNn(a));
      fListEventMixingBFS->Add(fShuffledBalance->GetHistNpn(a));
      fListEventMixingBFS->Add(fShuffledBalance->GetHistNnn(a));
      fListEventMixingBFS->Add(fShuffledBalance->GetHistNpp(a));
      fListEventMixingBFS->Add(fShuffledBalance->GetHistNnp(a));
    }  
  }

  if(fESDtrackCuts) fList->Add(fESDtrackCuts);

  //====================PID========================//
  if(fUsePID) {
    fPIDCombined = new AliPIDCombined();
    fPIDCombined->SetDefaultTPCPriors();

    fHistdEdxVsPTPCbeforePID = new TH2D ("dEdxVsPTPCbefore","dEdxVsPTPCbefore", 1000, -10.0, 10.0, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCbeforePID); //addition 
    
    fHistBetavsPTOFbeforePID = new TH2D ("BetavsPTOFbefore","BetavsPTOFbefore", 1000, -10.0, 10., 1000, 0, 1.2); 
    fHistListPIDQA->Add(fHistBetavsPTOFbeforePID); //addition
    
    fHistProbTPCvsPtbeforePID = new TH2D ("ProbTPCvsPtbefore","ProbTPCvsPtbefore", 1000, -10.0,10.0, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCvsPtbeforePID); //addition 
    
    fHistProbTOFvsPtbeforePID = new TH2D ("ProbTOFvsPtbefore","ProbTOFvsPtbefore", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTOFvsPtbeforePID); //addition 

    fHistProbTPCTOFvsPtbeforePID =new TH2D ("ProbTPCTOFvsPtbefore","ProbTPCTOFvsPtbefore", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCTOFvsPtbeforePID); //addition 
    
    fHistNSigmaTPCvsPtbeforePID = new TH2D ("NSigmaTPCvsPtbefore","NSigmaTPCvsPtbefore", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtbeforePID); //addition 
    
    fHistNSigmaTOFvsPtbeforePID = new TH2D ("NSigmaTOFvsPtbefore","NSigmaTOFvsPtbefore", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTOFvsPtbeforePID); //addition 
    
    fHistdEdxVsPTPCafterPID = new TH2D ("dEdxVsPTPCafter","dEdxVsPTPCafter", 1000, -10, 10, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCafterPID); //addition 
    
    fHistBetavsPTOFafterPID = new TH2D ("BetavsPTOFafter","BetavsPTOFafter", 1000, -10, 10, 1000, 0, 1.2); 
    fHistListPIDQA->Add(fHistBetavsPTOFafterPID); //addition 
    
    fHistProbTPCvsPtafterPID = new TH2D ("ProbTPCvsPtafter","ProbTPCvsPtafter", 1000, -10, 10, 1000, 0, 2); 
    fHistListPIDQA->Add(fHistProbTPCvsPtafterPID); //addition 
  
    fHistProbTOFvsPtafterPID = new TH2D ("ProbTOFvsPtafter","ProbTOFvsPtafter", 1000,  -10, 10, 1000, 0, 2); 
    fHistListPIDQA->Add(fHistProbTOFvsPtafterPID); //addition  
    
    fHistProbTPCTOFvsPtafterPID =new TH2D ("ProbTPCTOFvsPtafter","ProbTPCTOFvsPtafter", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCTOFvsPtafterPID); //addition 

    fHistNSigmaTPCvsPtafterPID = new TH2D ("NSigmaTPCvsPtafter","NSigmaTPCvsPtafter", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtafterPID); //addition  
    
    fHistNSigmaTOFvsPtafterPID = new TH2D ("NSigmaTOFvsPtafter","NSigmaTOFvsPtafter", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTOFvsPtafterPID); //addition 
  }
  //====================PID========================//

  // Post output data.
  PostData(1, fList);
  PostData(2, fListEventMixingBF);
  if(fRunShuffling) PostData(3, fListEventMixingBFS);
  if(fUsePID) PostData(4, fHistListPIDQA);       //PID
}

//________________________________________________________________________
void AliAnalysisTaskEventMixingBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  // NOTHING TO DO for event mixing!
}      

//________________________________________________________________________
void  AliAnalysisTaskEventMixingBF::FinishTaskOutput(){
  //Printf("END EventMixingBF");

  if (!fBalance) {
    Printf("ERROR: fBalance not available");
    return;
  }  
  if(fRunShuffling) {
    if (!fShuffledBalance) {
      Printf("ERROR: fShuffledBalance not available");
      return;
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskEventMixingBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}

void AliAnalysisTaskEventMixingBF::UserExecMix(Option_t *)
{

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  AliMixInputEventHandler *mixIEH = SetupEventsForMixing();

  Float_t fCentrality           = 0.;
  
  // vector holding the charges/kinematics of all tracks (charge,y,eta,phi,p0,p1,p2,pt,E)
  vector<Double_t> *chargeVector[9];          // original charge
  for(Int_t i = 0; i < 9; i++){
    chargeVector[i]        = new vector<Double_t>;
  }
  
  Double_t v_charge;
  Double_t v_y;
  Double_t v_eta;
  Double_t v_phi;
  Double_t v_p[3];
  Double_t v_pt;
  Double_t v_E;

  Int_t iMainTrackUsed = -1;

  // -------------------------------------------------------------		     
  // At the moment MIXING only for AODs
  if(mixIEH){

    //AOD analysis (vertex and track cuts also here!!!!)
    if(gAnalysisLevel == "AOD") {
      AliAODEvent* aodEventMain = dynamic_cast<AliAODEvent*>(fMainEvent); 
      if(!aodEventMain) {
  	Printf("ERROR: aodEventMain not available");
  	return;
      }
      AliAODEvent *aodEventMix  = dynamic_cast<AliAODEvent *>(fMixEvent); 
     if(!aodEventMix) {
  	Printf("ERROR: aodEventMix not available");
  	return;
      }
      
     AliAODHeader *aodHeaderMain = aodEventMain->GetHeader();
     AliAODHeader *aodHeaderMix  = aodEventMix->GetHeader();    
  

      // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
      fHistEventStats->Fill(1); //all events

      // this is not needed (checked in mixing handler!)
      Bool_t isSelectedMain = kTRUE;
      Bool_t isSelectedMix = kTRUE;
      
      if(fUseOfflineTrigger){
       	isSelectedMain = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	isSelectedMix = ((AliInputEventHandler*)((AliMultiInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetFirstMultiInputHandler())->IsEventSelected();
      }
      
      if(isSelectedMain && isSelectedMix) {
       	fHistEventStats->Fill(2); //triggered events
	
      	//Centrality stuff (centrality in AOD header)
      	if(fUseCentrality) {
      	  fCentrality = aodHeaderMain->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
	  
      	  // QA for centrality estimators
      	  fHistCentStats->Fill(0.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("V0M"));
      	  fHistCentStats->Fill(1.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("FMD"));
      	  fHistCentStats->Fill(2.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("TRK"));
      	  fHistCentStats->Fill(3.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("TKL"));
      	  fHistCentStats->Fill(4.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("CL0"));
      	  fHistCentStats->Fill(5.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("CL1"));
      	  fHistCentStats->Fill(6.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("V0MvsFMD"));
      	  fHistCentStats->Fill(7.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("TKLvsV0M"));
      	  fHistCentStats->Fill(8.,aodHeaderMain->GetCentralityP()->GetCentralityPercentile("ZEMvsZDC"));
	  
      	  // take only events inside centrality class
      	  if((fCentrality < fCentralityPercentileMin) || (fCentrality > fCentralityPercentileMax)) 
      	    return;
	  
      	  // centrality QA (V0M)
      	  fHistV0M->Fill(aodEventMain->GetVZEROData()->GetMTotV0A(), aodEventMain->GetVZEROData()->GetMTotV0C());
	  
      	  // centrality QA (reference tracks)
      	  fHistRefTracks->Fill(0.,aodHeaderMain->GetRefMultiplicity());
      	  fHistRefTracks->Fill(1.,aodHeaderMain->GetRefMultiplicityPos());
      	  fHistRefTracks->Fill(2.,aodHeaderMain->GetRefMultiplicityNeg());
      	  fHistRefTracks->Fill(3.,aodHeaderMain->GetTPConlyRefMultiplicity());
      	  fHistRefTracks->Fill(4.,aodHeaderMain->GetNumberOfITSClusters(0));
      	  fHistRefTracks->Fill(5.,aodHeaderMain->GetNumberOfITSClusters(1));
      	  fHistRefTracks->Fill(6.,aodHeaderMain->GetNumberOfITSClusters(2));
      	  fHistRefTracks->Fill(7.,aodHeaderMain->GetNumberOfITSClusters(3));
      	  fHistRefTracks->Fill(8.,aodHeaderMain->GetNumberOfITSClusters(4));
      	}
	
	// // this is crashing (Bug in ROOT (to be solved)) but not needed (checked in mixing handler)
	// const AliAODVertex *vertexMain = aodEventMain->GetPrimaryVertex();
	// const AliAODVertex *vertexMix  = aodEventMix->GetPrimaryVertex();
      	
	// if(vertexMain && vertexMix) {
      	//    Double32_t fCovMain[6];
      	//    Double32_t fCovMix[6];
      	//    vertexMain->GetCovarianceMatrix(fCovMain);
      	//    vertexMix->GetCovarianceMatrix(fCovMix);
	  
	//    if(vertexMain->GetNContributors() > 0 && vertexMix->GetNContributors() > 0) {
      	//     if(fCovMain[5] != 0 && fCovMix[5] != 0) {
      	//       fHistEventStats->Fill(3); //events with a proper vertex
      	//       if(TMath::Abs(vertexMain->GetX()) < fVxMax && TMath::Abs(vertexMix->GetX()) < fVxMax ) {
      	// 	if(TMath::Abs(vertexMain->GetY()) < fVyMax && TMath::Abs(vertexMix->GetY()) < fVyMax) {
      	// 	  if(TMath::Abs(vertexMain->GetZ()) < fVzMax && TMath::Abs(vertexMix->GetZ()) < fVzMax) {
      	// 	    fHistEventStats->Fill(4); //analyzed events
      	// 	    fHistVx->Fill(vertexMain->GetX());
      	// 	    fHistVy->Fill(vertexMain->GetY());
      	// 	    fHistVz->Fill(vertexMain->GetZ());

  		    // Loop over tracks in main event
  		    for (Int_t iTracksMain = 0; iTracksMain < aodEventMain->GetNumberOfTracks(); iTracksMain++) {
  		      AliAODTrack* aodTrackMain = dynamic_cast<AliAODTrack *>(aodEventMain->GetTrack(iTracksMain));
  		      if (!aodTrackMain) {
  			Printf("ERROR: Could not receive track %d", iTracksMain);
  			continue;
  		      }
		      
  		      // AOD track cuts
		      
  		      // For ESD Filter Information: ANALYSIS/macros/AddTaskESDfilter.C
  		      // take only TPC only tracks 
  		      fHistTrackStats->Fill(aodTrackMain->GetFilterMap());
  		      if(!aodTrackMain->TestFilterBit(nAODtrackCutBit)) continue;
		      
  		      v_charge = aodTrackMain->Charge();
  		      v_y      = aodTrackMain->Y();
  		      v_eta    = aodTrackMain->Eta();
  		      v_phi    = aodTrackMain->Phi() * TMath::RadToDeg();
  		      v_E      = aodTrackMain->E();
  		      v_pt     = aodTrackMain->Pt();
  		      aodTrackMain->PxPyPz(v_p);
		      
  		      Float_t DCAxyMain = aodTrackMain->DCA();      // this is the DCA from global track (not exactly what is cut on)
  		      Float_t DCAzMain  = aodTrackMain->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
		      
		      
  		      // Kinematics cuts from ESD track cuts
  		      if( v_pt < fPtMin || v_pt > fPtMax)      continue;
  		      if( v_eta < fEtaMin || v_eta > fEtaMax)  continue;
		      
  		      // Extra DCA cuts (for systematic studies [!= -1])
  		      if( fDCAxyCut != -1 && fDCAzCut != -1){
  			if(TMath::Sqrt((DCAxyMain*DCAxyMain)/(fDCAxyCut*fDCAxyCut)+(DCAzMain*DCAzMain)/(fDCAzCut*fDCAzCut)) > 1 ){
  			  continue;  // 2D cut
  			}
  		      }
		      
  		      // Extra TPC cuts (for systematic studies [!= -1])
  		      if( fTPCchi2Cut != -1 && aodTrackMain->Chi2perNDF() > fTPCchi2Cut){
  			continue;
  		      }
  		      if( fNClustersTPCCut != -1 && aodTrackMain->GetTPCNcls() < fNClustersTPCCut){
  			continue;
  		      }
		      
  		      // fill QA histograms
  		      fHistClus->Fill(aodTrackMain->GetITSNcls(),aodTrackMain->GetTPCNcls());
  		      fHistDCA->Fill(DCAzMain,DCAxyMain);
  		      fHistChi2->Fill(aodTrackMain->Chi2perNDF());
  		      fHistPt->Fill(v_pt);
  		      fHistEta->Fill(v_eta);
  		      fHistPhi->Fill(v_phi);
		      
  		      // fill charge vector
  		      chargeVector[0]->push_back(v_charge);
  		      chargeVector[1]->push_back(v_y);
  		      chargeVector[2]->push_back(v_eta);
  		      chargeVector[3]->push_back(v_phi);
  		      chargeVector[4]->push_back(v_p[0]);
  		      chargeVector[5]->push_back(v_p[1]);
  		      chargeVector[6]->push_back(v_p[2]);
  		      chargeVector[7]->push_back(v_pt);
  		      chargeVector[8]->push_back(v_E);

  		      // -------------------------------------------------------------		     
  		      // for each track in main event loop over all tracks in mix event
  		      for (Int_t iTracksMix = 0; iTracksMix < aodEventMix->GetNumberOfTracks(); iTracksMix++) {

  			AliAODTrack* aodTrackMix = dynamic_cast<AliAODTrack *>(aodEventMix->GetTrack(iTracksMix));
  			if (!aodTrackMix) {
  			  Printf("ERROR: Could not receive track %d", iTracksMix);
  			  continue;
  			}

  			// AOD track cuts
			
  			// For ESD Filter Information: ANALYSIS/macros/AddTaskESDfilter.C
  			// take only TPC only tracks 
  			fHistTrackStats->Fill(aodTrackMix->GetFilterMap());
  			if(!aodTrackMix->TestFilterBit(nAODtrackCutBit)) continue;
			
  			v_charge = aodTrackMix->Charge();
  			v_y      = aodTrackMix->Y();
  			v_eta    = aodTrackMix->Eta();
  			v_phi    = aodTrackMix->Phi() * TMath::RadToDeg();
  			v_E      = aodTrackMix->E();
  			v_pt     = aodTrackMix->Pt();
  			aodTrackMix->PxPyPz(v_p);
		      
  			Float_t DCAxyMix = aodTrackMix->DCA();      // this is the DCA from global track (not exactly what is cut on)
  			Float_t DCAzMix  = aodTrackMix->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
			
			
  			// Kinematics cuts from ESD track cuts
  			if( v_pt < fPtMin || v_pt > fPtMax)      continue;
  			if( v_eta < fEtaMin || v_eta > fEtaMax)  continue;
			
  			// Extra DCA cuts (for systematic studies [!= -1])
  			if( fDCAxyCut != -1 && fDCAxyCut != -1){
  			  if(TMath::Sqrt((DCAxyMix*DCAxyMix)/(fDCAxyCut*fDCAxyCut)+(DCAzMix*DCAzMix)/(fDCAzCut*fDCAzCut)) > 1 ){
  			    continue;  // 2D cut
  			  }
  			}
			
  			// Extra TPC cuts (for systematic studies [!= -1])
  			if( fTPCchi2Cut != -1 && aodTrackMix->Chi2perNDF() > fTPCchi2Cut){
  			  continue;
  			}
  			if( fNClustersTPCCut != -1 && aodTrackMix->GetTPCNcls() < fNClustersTPCCut){
  			  continue;
  			}
			
  			// fill QA histograms
  			fHistClus->Fill(aodTrackMix->GetITSNcls(),aodTrackMix->GetTPCNcls());
  			fHistDCA->Fill(DCAzMix,DCAxyMix);
  			fHistChi2->Fill(aodTrackMix->Chi2perNDF());
  			fHistPt->Fill(v_pt);
  			fHistEta->Fill(v_eta);
  			fHistPhi->Fill(v_phi);
			
  			// fill charge vector
  			chargeVector[0]->push_back(v_charge);
  			chargeVector[1]->push_back(v_y);
  			chargeVector[2]->push_back(v_eta);
  			chargeVector[3]->push_back(v_phi);
  			chargeVector[4]->push_back(v_p[0]);
  			chargeVector[5]->push_back(v_p[1]);
  			chargeVector[6]->push_back(v_p[2]);
  			chargeVector[7]->push_back(v_pt);
  			chargeVector[8]->push_back(v_E);
			
			
			
  		      } //mix track loop

  		      // calculate balance function for each track in main event
		      iMainTrackUsed++; // is needed to do no double counting in Balance Function calculation   
		      if(iMainTrackUsed >= (Int_t)chargeVector[0]->size()) break; //do not allow more tracks than in mixed event!
		      fBalance->CalculateBalance(fCentrality,chargeVector,iMainTrackUsed);
  		      // clean charge vector afterwards
  		      for(Int_t i = 0; i < 9; i++){		       
  			chargeVector[i]->clear();
  		      }
		      

  		    } //main track loop
      // 		  }//Vz cut
      // 		}//Vy cut
      // 	      }//Vx cut
      // 	    }//proper vertexresolution
      //    }//proper number of contributors
      // }//vertex object valid
      }//triggered event 
    }//AOD analysis
  }
}

AliMixInputEventHandler *AliAnalysisTaskEventMixingBF::SetupEventsForMixing() {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliMultiInputEventHandler *inEvHMain = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
  if (inEvHMain) {

      AliMixInputEventHandler *mixEH = dynamic_cast<AliMixInputEventHandler *>(inEvHMain->GetFirstMultiInputHandler());
      if (!mixEH) return kFALSE;
      if (mixEH->CurrentBinIndex() < 0) {
         AliDebug(AliLog::kDebug + 1, "Current event mixEH->CurrentEntry() == -1");
         return kFALSE ;
      }
      AliDebug(AliLog::kDebug, Form("Mixing %lld %d [%lld,%lld] %d", mixEH->CurrentEntry(), mixEH->CurrentBinIndex(), mixEH->CurrentEntryMain(), mixEH->CurrentEntryMix(), mixEH->NumberMixed()));

      AliInputEventHandler      *ihMainCurrent     = inEvHMain->GetFirstInputEventHandler();
      fMainEvent = ihMainCurrent->GetEvent();

      AliMultiInputEventHandler *inEvHMixedCurrent = mixEH->GetFirstMultiInputHandler(); // for buffer = 1
      AliInputEventHandler      *ihMixedCurrent    = inEvHMixedCurrent->GetFirstInputEventHandler();
      fMixEvent                                    = ihMixedCurrent->GetEvent();
      
      return mixEH;
  }
  return NULL;
} 
