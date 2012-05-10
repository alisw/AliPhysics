#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"                  
#include "TH3D.h"
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
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "AliTHn.h"              

#include "AliPID.h"                
#include "AliPIDResponse.h"        
#include "AliPIDCombined.h"        

#include "AliAnalysisTaskBFPsi.h"
#include "AliBalancePsi.h"


// Analysis task for the BF vs Psi code
// Authors: Panos.Christakoglou@nikhef.nl

ClassImp(AliAnalysisTaskBFPsi)

//________________________________________________________________________
AliAnalysisTaskBFPsi::AliAnalysisTaskBFPsi(const char *name) 
: AliAnalysisTaskSE(name), 
  fBalance(0),
  fRunShuffling(kFALSE),
  fShuffledBalance(0),
  fList(0),
  fListBF(0),
  fListBFS(0),
  fHistListPIDQA(0),
  fHistEventStats(0),
  fHistCentStats(0),
  fHistTriggerStats(0),
  fHistTrackStats(0),
  fHistVx(0),
  fHistVy(0),
  fHistVz(0),
  fHistEventPlane(0),
  fHistClus(0),
  fHistDCA(0),
  fHistChi2(0),
  fHistPt(0),
  fHistEta(0),
  fHistRapidity(0),
  fHistPhi(0),
  fHistPhiBefore(0),
  fHistPhiAfter(0),
  fHistPhiPos(0),
  fHistPhiNeg(0),
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
  fPDGCodeToBeAnalyzed(-1) {
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
AliAnalysisTaskBFPsi::~AliAnalysisTaskBFPsi() {

  // delete fBalance; 
  // delete fShuffledBalance; 
  // delete fList;
  // delete fListBF; 
  // delete fListBFS;

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
void AliAnalysisTaskBFPsi::UserCreateOutputObjects() {
  // Create histograms
  // Called once
  if(!fBalance) {
    fBalance = new AliBalancePsi();
    fBalance->SetAnalysisLevel("ESD");
    //fBalance->SetNumberOfBins(-1,16);
    //fBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6,15.);
  }
  if(fRunShuffling) {
    if(!fShuffledBalance) {
      fShuffledBalance = new AliBalancePsi();
      fShuffledBalance->SetAnalysisLevel("ESD");
      //fShuffledBalance->SetNumberOfBins(-1,16);
      //fShuffledBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6,15.);
    }
  }

  //QA list
  fList = new TList();
  fList->SetName("listQA");
  fList->SetOwner();

  //Balance Function list
  fListBF = new TList();
  fListBF->SetName("listBF");
  fListBF->SetOwner();

  if(fRunShuffling) {
    fListBFS = new TList();
    fListBFS->SetName("listBFShuffled");
    fListBFS->SetOwner();
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
  fHistEventStats = new TH2F("fHistEventStats",
                             "Event statistics;;Centrality percentile;N_{events}",
                             4,0.5,4.5,220,-5,105);
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

  fHistNumberOfAcceptedTracks = new TH2F("fHistNumberOfAcceptedTracks",";N_{acc.};Centrality percentile;Entries",4001,-0.5,4000.5,220,-5,105);
  fList->Add(fHistNumberOfAcceptedTracks);

  // Vertex distributions
  fHistVx = new TH1F("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVx);
  fHistVy = new TH1F("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVy);
  fHistVz = new TH2F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Centrality percentile;Entries",100,-20.,20.,220,-5,105);
  fList->Add(fHistVz);

  //Event plane
  fHistEventPlane = new TH2F("fHistEventPlane",";#Psi_{2} [deg.];Centrality percentile;Counts",100,0,360.,220,-5,105);
  fList->Add(fHistEventPlane);

  // QA histograms
  fHistClus = new TH2F("fHistClus","# Cluster (TPC vs. ITS)",10,0,10,200,0,200);
  fList->Add(fHistClus);
  fHistChi2 = new TH2F("fHistChi2","Chi2/NDF distribution;#chi^{2}/ndf;Centrality percentile",200,0,10,220,-5,105);
  fList->Add(fHistChi2);
  fHistDCA  = new TH2F("fHistDCA","DCA (xy vs. z)",400,-5,5,400,-5,5); 
  fList->Add(fHistDCA);
  fHistPt   = new TH2F("fHistPt","p_{T} distribution;p_{T} (GeV/c);Centrality percentile",200,0,10,220,-5,105);
  fList->Add(fHistPt);
  fHistEta  = new TH2F("fHistEta","#eta distribution;#eta;Centrality percentile",200,-2,2,220,-5,105);
  fList->Add(fHistEta);
  fHistRapidity  = new TH2F("fHistRapidity","y distribution;y;Centrality percentile",200,-2,2,220,-5,105);
  fList->Add(fHistRapidity);
  fHistPhi  = new TH2F("fHistPhi","#phi distribution;#phi;Centrality percentile",200,-20,380,220,-5,105);
  fList->Add(fHistPhi);
  fHistPhiBefore  = new TH2F("fHistPhiBefore","#phi distribution;#phi;Centrality percentile",200,0.,2*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhiBefore);
  fHistPhiAfter  = new TH2F("fHistPhiAfter","#phi distribution;#phi;Centrality percentile",200,0.,2*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhiAfter);
  fHistPhiPos  = new TH2F("fHistPhiPos","#phi distribution for positive particles;#phi;Centrality percentile",200,-20,380,220,-5,105);
  fList->Add(fHistPhiPos);
  fHistPhiNeg  = new TH2F("fHistPhiNeg","#phi distribution for negative particles;#phi;Centrality percentile",200,-20,380,220,-5,105);
  fList->Add(fHistPhiNeg);
  fHistV0M  = new TH2F("fHistV0M","V0 Multiplicity C vs. A",500, 0, 20000, 500, 0, 20000);
  fList->Add(fHistV0M);
  TString gRefTrackName[6] = {"tracks","tracksPos","tracksNeg","tracksTPConly","clusITS0","clusITS1"};
  fHistRefTracks  = new TH2F("fHistRefTracks","Nr of Ref tracks/event vs. ref track estimator;;Nr of tracks",6, 0, 6, 400, 0, 20000);
  for(Int_t i = 1; i <= 6; i++)
    fHistRefTracks->GetXaxis()->SetBinLabel(i,gRefTrackName[i-1].Data());
  fList->Add(fHistRefTracks);

  // Balance function histograms
  // Initialize histograms if not done yet
  if(!fBalance->GetHistNp()){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    fBalance->InitHistograms();
  }

  if(fRunShuffling) {
    if(!fShuffledBalance->GetHistNp()) {
      AliWarning("Histograms (shuffling) not yet initialized! --> Will be done now");
      AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
      fShuffledBalance->InitHistograms();
    }
  }

  //for(Int_t a = 0; a < ANALYSIS_TYPES; a++){
  fListBF->Add(fBalance->GetHistNp());
  fListBF->Add(fBalance->GetHistNn());
  fListBF->Add(fBalance->GetHistNpn());
  fListBF->Add(fBalance->GetHistNnn());
  fListBF->Add(fBalance->GetHistNpp());
  fListBF->Add(fBalance->GetHistNnp());

  if(fRunShuffling) {
    fListBFS->Add(fShuffledBalance->GetHistNp());
    fListBFS->Add(fShuffledBalance->GetHistNn());
    fListBFS->Add(fShuffledBalance->GetHistNpn());
    fListBFS->Add(fShuffledBalance->GetHistNnn());
    fListBFS->Add(fShuffledBalance->GetHistNpp());
    fListBFS->Add(fShuffledBalance->GetHistNnp());
  }  
  //}

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
  PostData(2, fListBF);
  if(fRunShuffling) PostData(3, fListBFS);
  if(fUsePID) PostData(4, fHistListPIDQA);       //PID
}

//________________________________________________________________________
void AliAnalysisTaskBFPsi::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  AliESDtrack *track_TPC   = NULL;

  Int_t gNumberOfAcceptedTracks = 0;
  Float_t fCentrality           = 0.;
  Double_t gReactionPlane       = 0.;

  // vector holding the charges/kinematics of all tracks (charge,y,eta,phi,p0,p1,p2,pt,E)
  vector<Double_t> *chargeVectorShuffle[9];   // this will be shuffled
  vector<Double_t> *chargeVector[9];          // original charge
  for(Int_t i = 0; i < 9; i++){
    chargeVectorShuffle[i] = new vector<Double_t>;
    chargeVector[i]        = new vector<Double_t>;
  }

  Double_t v_charge;
  Double_t v_y;
  Double_t v_eta;
  Double_t v_phi;
  Double_t v_p[3];
  Double_t v_pt;
  Double_t v_E;

  if(fUsePID) {
    fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
  }
 
  //ESD analysis
  if(gAnalysisLevel == "ESD") {
    AliESDEvent* gESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from TaskSE
    if (!gESD) {
      Printf("ERROR: gESD not available");
      return;
    }

    // store offline trigger bits
    fHistTriggerStats->Fill(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());

    // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
    fHistEventStats->Fill(1,fCentrality); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2,fCentrality); //triggered events

      if(fUseCentrality) {
	//Centrality stuff
	AliCentrality *centrality = gESD->GetCentrality();
	
	fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
	
	// take only events inside centrality class
	if(!centrality->IsEventInCentralityClass(fCentralityPercentileMin,
						 fCentralityPercentileMax,
						 fCentralityEstimator.Data()))
	  return;
	
	// centrality QA (V0M)
	fHistV0M->Fill(gESD->GetVZEROData()->GetMTotV0A(), gESD->GetVZEROData()->GetMTotV0C());
      }
	
      const AliESDVertex *vertex = gESD->GetPrimaryVertex();
      if(vertex) {
	if(vertex->GetNContributors() > 0) {
	  if(vertex->GetZRes() != 0) {
	    fHistEventStats->Fill(3,fCentrality); //events with a proper vertex
	    if(TMath::Abs(vertex->GetXv()) < fVxMax) {
	      if(TMath::Abs(vertex->GetYv()) < fVyMax) {
		if(TMath::Abs(vertex->GetZv()) < fVzMax) {
		  fHistEventStats->Fill(4,fCentrality); //analayzed events
		  fHistVx->Fill(vertex->GetXv());
		  fHistVy->Fill(vertex->GetYv());
		  fHistVz->Fill(vertex->GetZv(),fCentrality);
		  
		  //========Get the VZERO event plane========//
		  Double_t gVZEROEventPlane = -10.0;
		  Double_t qxTot = 0.0, qyTot = 0.0;
		  AliEventplane *ep = gESD->GetEventplane();
		  if(ep) 
		    gVZEROEventPlane = ep->CalculateVZEROEventPlane(gESD,10,2,qxTot,qyTot);
		  if(gVZEROEventPlane < 0.) gVZEROEventPlane += TMath::Pi();
		  gReactionPlane = gVZEROEventPlane*TMath::RadToDeg();
		  fHistEventPlane->Fill(gReactionPlane,fCentrality);
		  //========Get the VZERO event plane========//

		  //Printf("There are %d tracks in this event", gESD->GetNumberOfTracks());
		  for (Int_t iTracks = 0; iTracks < gESD->GetNumberOfTracks(); iTracks++) {
		    AliESDtrack* track = dynamic_cast<AliESDtrack *>(gESD->GetTrack(iTracks));
		    if (!track) {
		      Printf("ERROR: Could not receive track %d", iTracks);
		      continue;
		    }	
		    
		    // take only TPC only tracks
		    track_TPC   = new AliESDtrack();
		    if(!track->FillTPCOnlyTrack(*track_TPC)) continue;
		    
		    //ESD track cuts
		    if(fESDtrackCuts) 
		      if(!fESDtrackCuts->AcceptTrack(track_TPC)) continue;
		    
		    // fill QA histograms
		    Float_t b[2];
		    Float_t bCov[3];
		    track_TPC->GetImpactParameters(b,bCov);
		    if (bCov[0]<=0 || bCov[2]<=0) {
		      AliDebug(1, "Estimated b resolution lower or equal zero!");
		      bCov[0]=0; bCov[2]=0;
		    }
		    
		    Int_t nClustersTPC = -1;
		    nClustersTPC = track_TPC->GetTPCNclsIter1();   // TPC standalone
		    //nClustersTPC = track->GetTPCclusters(0);   // global track
		    Float_t chi2PerClusterTPC = -1;
		    if (nClustersTPC!=0) {
		      chi2PerClusterTPC = track_TPC->GetTPCchi2Iter1()/Float_t(nClustersTPC);      // TPC standalone
		      //chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);     // global track
		    }

		    //===========================PID===============================//		    
		    if(fUsePID) {
		      Double_t prob[AliPID::kSPECIES]={0.};
		      Double_t probTPC[AliPID::kSPECIES]={0.};
		      Double_t probTOF[AliPID::kSPECIES]={0.};
		      Double_t probTPCTOF[AliPID::kSPECIES]={0.};

		      Double_t nSigma = 0.;
                      UInt_t detUsedTPC = 0;
		      UInt_t detUsedTOF = 0;
                      UInt_t detUsedTPCTOF = 0;

		      //Decide what detector configuration we want to use
		      switch(fPidDetectorConfig) {
		      case kTPCpid:
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
			nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)fParticleOfInterest));
			detUsedTPC = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
			for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
			  prob[iSpecies] = probTPC[iSpecies];
			break;
		      case kTOFpid:
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
			nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)fParticleOfInterest));
			detUsedTOF = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
			for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
			  prob[iSpecies] = probTOF[iSpecies];
			break;
		      case kTPCTOF:
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
			detUsedTPCTOF = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPCTOF);
			for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
			  prob[iSpecies] = probTPCTOF[iSpecies];
			break;
		      default:
			break;
		      }//end switch: define detector mask
		      
		      //Filling the PID QA
		      Double_t tofTime = -999., length = 999., tof = -999.;
		      Double_t c = TMath::C()*1.E-9;// m/ns
		      Double_t beta = -999.;
		      Double_t  nSigmaTOFForParticleOfInterest = -999.;
		      if ( (track->IsOn(AliESDtrack::kTOFin)) &&
			   (track->IsOn(AliESDtrack::kTIME))  ) { 
			tofTime = track->GetTOFsignal();//in ps
			length = track->GetIntegratedLength();
			tof = tofTime*1E-3; // ns	
			
			if (tof <= 0) {
			  //Printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
			  continue;
			}
			if (length <= 0){
			  //printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
			  continue;
			}
			
			length = length*0.01; // in meters
			tof = tof*c;
			beta = length/tof;
			
			nSigmaTOFForParticleOfInterest = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)fParticleOfInterest);
			fHistBetavsPTOFbeforePID ->Fill(track->P()*track->Charge(),beta);
			fHistProbTOFvsPtbeforePID ->Fill(track->Pt(),probTOF[fParticleOfInterest]);
			fHistNSigmaTOFvsPtbeforePID ->Fill(track->Pt(),nSigmaTOFForParticleOfInterest);
		      }//TOF signal 
		      
		      
		      Double_t  nSigmaTPCForParticleOfInterest = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)fParticleOfInterest);
		      fHistdEdxVsPTPCbeforePID -> Fill(track->P()*track->Charge(),track->GetTPCsignal());
		      fHistProbTPCvsPtbeforePID -> Fill(track->Pt(),probTPC[fParticleOfInterest]); 
		      fHistNSigmaTPCvsPtbeforePID -> Fill(track->Pt(),nSigmaTPCForParticleOfInterest); 
		      fHistProbTPCTOFvsPtbeforePID -> Fill(track->Pt(),probTPCTOF[fParticleOfInterest]);
		      //end of QA-before pid
		      
		      if ((detUsedTPC != 0)||(detUsedTOF != 0)||(detUsedTPCTOF != 0)) {
			//Make the decision based on the n-sigma
			if(fUsePIDnSigma) {
			  if(nSigma > fPIDNSigma) continue;}
			
			//Make the decision based on the bayesian
			else if(fUsePIDPropabilities) {
			  if(fParticleOfInterest != TMath::LocMax(AliPID::kSPECIES,prob)) continue;
			  if (prob[fParticleOfInterest] < fMinAcceptedPIDProbability) continue;      
			}
			
			//Fill QA after the PID
			fHistBetavsPTOFafterPID ->Fill(track->P()*track->Charge(),beta);
			fHistProbTOFvsPtafterPID ->Fill(track->Pt(),probTOF[fParticleOfInterest]);
			fHistNSigmaTOFvsPtafterPID ->Fill(track->Pt(),nSigmaTOFForParticleOfInterest);
			
			fHistdEdxVsPTPCafterPID -> Fill(track->P()*track->Charge(),track->GetTPCsignal());
			fHistProbTPCvsPtafterPID -> Fill(track->Pt(),probTPC[fParticleOfInterest]); 
			fHistProbTPCTOFvsPtafterPID -> Fill(track->Pt(),probTPCTOF[fParticleOfInterest]);
			fHistNSigmaTPCvsPtafterPID -> Fill(track->Pt(),nSigmaTPCForParticleOfInterest); 
		      }
		      
		      PostData(4, fHistListPIDQA);
		    }
                    //===========================PID===============================//
		    v_charge = track_TPC->Charge();
		    v_y      = track_TPC->Y();
		    v_eta    = track_TPC->Eta();
		    v_phi    = track_TPC->Phi() * TMath::RadToDeg();
		    v_E      = track_TPC->E();
		    v_pt     = track_TPC->Pt();
		    track_TPC->PxPyPz(v_p);
		    fHistClus->Fill(track_TPC->GetITSclusters(0),nClustersTPC);
		    fHistDCA->Fill(b[1],b[0]);
		    fHistChi2->Fill(chi2PerClusterTPC,fCentrality);
		    fHistPt->Fill(v_pt,fCentrality);
		    fHistEta->Fill(v_eta,fCentrality);
		    fHistPhi->Fill(v_phi,fCentrality);
		    fHistRapidity->Fill(v_y,fCentrality);
		    if(v_charge > 0) fHistPhiPos->Fill(v_phi,fCentrality);
		    else if(v_charge < 0) fHistPhiNeg->Fill(v_phi,fCentrality);

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

		    if(fRunShuffling) {
		      chargeVectorShuffle[0]->push_back(v_charge);
		      chargeVectorShuffle[1]->push_back(v_y);
		      chargeVectorShuffle[2]->push_back(v_eta);
		      chargeVectorShuffle[3]->push_back(v_phi);
		      chargeVectorShuffle[4]->push_back(v_p[0]);
		      chargeVectorShuffle[5]->push_back(v_p[1]);
		      chargeVectorShuffle[6]->push_back(v_p[2]);
		      chargeVectorShuffle[7]->push_back(v_pt);
		      chargeVectorShuffle[8]->push_back(v_E);
		    }
		    
		    delete track_TPC;
		    
		  } //track loop
		}//Vz cut
	      }//Vy cut
	    }//Vx cut
	  }//proper vertex resolution
	}//proper number of contributors
      }//vertex object valid
    }//triggered event 
  }//ESD analysis
  
  //AOD analysis (vertex and track cuts also here!!!!)
  else if(gAnalysisLevel == "AOD") {
    AliAODEvent* gAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
    if(!gAOD) {
      Printf("ERROR: gAOD not available");
      return;
    }

    AliAODHeader *aodHeader = gAOD->GetHeader();

    // store offline trigger bits
    fHistTriggerStats->Fill(aodHeader->GetOfflineTrigger());
    
    // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
    fHistEventStats->Fill(1,fCentrality); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2,fCentrality); //triggered events
		  
      //Centrality stuff (centrality in AOD header)
      if(fUseCentrality) {
    	fCentrality = aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
 	
    	// QA for centrality estimators
    	fHistCentStats->Fill(0.,aodHeader->GetCentralityP()->GetCentralityPercentile("V0M"));
    	fHistCentStats->Fill(1.,aodHeader->GetCentralityP()->GetCentralityPercentile("FMD"));
    	fHistCentStats->Fill(2.,aodHeader->GetCentralityP()->GetCentralityPercentile("TRK"));
    	fHistCentStats->Fill(3.,aodHeader->GetCentralityP()->GetCentralityPercentile("TKL"));
    	fHistCentStats->Fill(4.,aodHeader->GetCentralityP()->GetCentralityPercentile("CL0"));
    	fHistCentStats->Fill(5.,aodHeader->GetCentralityP()->GetCentralityPercentile("CL1"));
    	fHistCentStats->Fill(6.,aodHeader->GetCentralityP()->GetCentralityPercentile("V0MvsFMD"));
    	fHistCentStats->Fill(7.,aodHeader->GetCentralityP()->GetCentralityPercentile("TKLvsV0M"));
    	fHistCentStats->Fill(8.,aodHeader->GetCentralityP()->GetCentralityPercentile("ZEMvsZDC"));
	
    	// take only events inside centrality class
    	if((fCentrality < fCentralityPercentileMin) || (fCentrality > fCentralityPercentileMax)) 
    	  return;
	
    	// centrality QA (V0M)
    	fHistV0M->Fill(gAOD->GetVZEROData()->GetMTotV0A(), gAOD->GetVZEROData()->GetMTotV0C());
	
    	// centrality QA (reference tracks)
    	fHistRefTracks->Fill(0.,aodHeader->GetRefMultiplicity());
    	fHistRefTracks->Fill(1.,aodHeader->GetRefMultiplicityPos());
    	fHistRefTracks->Fill(2.,aodHeader->GetRefMultiplicityNeg());
    	fHistRefTracks->Fill(3.,aodHeader->GetTPConlyRefMultiplicity());
    	fHistRefTracks->Fill(4.,aodHeader->GetNumberOfITSClusters(0));
    	fHistRefTracks->Fill(5.,aodHeader->GetNumberOfITSClusters(1));
    	fHistRefTracks->Fill(6.,aodHeader->GetNumberOfITSClusters(2));
    	fHistRefTracks->Fill(7.,aodHeader->GetNumberOfITSClusters(3));
    	fHistRefTracks->Fill(8.,aodHeader->GetNumberOfITSClusters(4));
      }

      const AliAODVertex *vertex = gAOD->GetPrimaryVertex();
      
      if(vertex) {
      	Double32_t fCov[6];
      	vertex->GetCovarianceMatrix(fCov);
	
      	if(vertex->GetNContributors() > 0) {
      	  if(fCov[5] != 0) {
      	    fHistEventStats->Fill(3,fCentrality); //events with a proper vertex
      	    if(TMath::Abs(vertex->GetX()) < fVxMax) {
      	      if(TMath::Abs(vertex->GetY()) < fVyMax) {
      		if(TMath::Abs(vertex->GetZ()) < fVzMax) {
      		  fHistEventStats->Fill(4,fCentrality); //analyzed events
      		  fHistVx->Fill(vertex->GetX());
      		  fHistVy->Fill(vertex->GetY());
      		  fHistVz->Fill(vertex->GetZ(),fCentrality);
	
		  //========Get the VZERO event plane========//
		  Double_t gVZEROEventPlane = -10.0;
		  Double_t qxTot = 0.0, qyTot = 0.0;
		  AliEventplane *ep = gAOD->GetEventplane();
		  if(ep) 
		    gVZEROEventPlane = ep->CalculateVZEROEventPlane(gAOD,10,2,qxTot,qyTot);
		  if(gVZEROEventPlane < 0.) gVZEROEventPlane += TMath::Pi();
		  gReactionPlane = gVZEROEventPlane*TMath::RadToDeg();
		  fHistEventPlane->Fill(gReactionPlane,fCentrality);
		  //========Get the VZERO event plane========//

      		  //Printf("There are %d tracks in this event", gAOD->GetNumberOfTracks());
      		  for (Int_t iTracks = 0; iTracks < gAOD->GetNumberOfTracks(); iTracks++) {
      		    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(gAOD->GetTrack(iTracks));
      		    if (!aodTrack) {
      		      Printf("ERROR: Could not receive track %d", iTracks);
      		      continue;
      		    }
		    
      		    // AOD track cuts
		    
      		    // For ESD Filter Information: ANALYSIS/macros/AddTaskESDfilter.C
      		    // take only TPC only tracks 
      		    fHistTrackStats->Fill(aodTrack->GetFilterMap());
      		    if(!aodTrack->TestFilterBit(nAODtrackCutBit)) continue;
		    
      		    v_charge = aodTrack->Charge();
      		    v_y      = aodTrack->Y();
      		    v_eta    = aodTrack->Eta();
      		    v_phi    = aodTrack->Phi() * TMath::RadToDeg();
      		    v_E      = aodTrack->E();
      		    v_pt     = aodTrack->Pt();
      		    aodTrack->PxPyPz(v_p);
		    
      		    Float_t DCAxy = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
      		    Float_t DCAz  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
		    
		    
      		    // Kinematics cuts from ESD track cuts
      		    if( v_pt < fPtMin || v_pt > fPtMax)      continue;
		    if (!fUsePID) {
		      if( v_eta < fEtaMin || v_eta > fEtaMax)  continue;
		    }
		    else if (fUsePID){
		      if( v_y < fEtaMin || v_y > fEtaMax)  continue;
		    }

      		    // Extra DCA cuts (for systematic studies [!= -1])
      		    if( fDCAxyCut != -1 && fDCAzCut != -1){
      		      if(TMath::Sqrt((DCAxy*DCAxy)/(fDCAxyCut*fDCAxyCut)+(DCAz*DCAz)/(fDCAzCut*fDCAzCut)) > 1 ){
      			continue;  // 2D cut
      		      }
      		    }
		    
      		    // Extra TPC cuts (for systematic studies [!= -1])
      		    if( fTPCchi2Cut != -1 && aodTrack->Chi2perNDF() > fTPCchi2Cut){
      		      continue;
      		    }
      		    if( fNClustersTPCCut != -1 && aodTrack->GetTPCNcls() < fNClustersTPCCut){
      		      continue;
      		    }
		    
      		    // fill QA histograms
      		    fHistClus->Fill(aodTrack->GetITSNcls(),aodTrack->GetTPCNcls());
      		    fHistDCA->Fill(DCAz,DCAxy);
      		    fHistChi2->Fill(aodTrack->Chi2perNDF(),fCentrality);
      		    fHistPt->Fill(v_pt,fCentrality);
      		    fHistEta->Fill(v_eta,fCentrality);
      		    fHistPhi->Fill(v_phi,fCentrality);
		    fHistRapidity->Fill(v_y,fCentrality);
		    if(v_charge > 0) fHistPhiPos->Fill(v_phi,fCentrality);
		    else if(v_charge < 0) fHistPhiNeg->Fill(v_phi,fCentrality);

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

      		    if(fRunShuffling) {
      		      chargeVectorShuffle[0]->push_back(v_charge);
      		      chargeVectorShuffle[1]->push_back(v_y);
      		      chargeVectorShuffle[2]->push_back(v_eta);
      		      chargeVectorShuffle[3]->push_back(v_phi);
      		      chargeVectorShuffle[4]->push_back(v_p[0]);
      		      chargeVectorShuffle[5]->push_back(v_p[1]);
      		      chargeVectorShuffle[6]->push_back(v_p[2]);
      		      chargeVectorShuffle[7]->push_back(v_pt);
      		      chargeVectorShuffle[8]->push_back(v_E);
      		    }
		    		    
      		    gNumberOfAcceptedTracks += 1;
		    
      		  } //track loop
      		}//Vz cut
      	      }//Vy cut
      	    }//Vx cut
      	  }//proper vertex resolution
	}//proper number of contributors
      }//vertex object valid
    }//triggered event 
  }//AOD analysis

  //MC-ESD analysis
  if(gAnalysisLevel == "MCESD") {
    AliMCEvent*  mcEvent = MCEvent(); 
    if (!mcEvent) {
      Printf("ERROR: mcEvent not available");
      return;
    }

    AliESDEvent* gESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from TaskSE
    if (!gESD) {
      Printf("ERROR: gESD not available");
      return;
    }

    // store offline trigger bits
    fHistTriggerStats->Fill(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());

    // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
    fHistEventStats->Fill(1,fCentrality); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2,fCentrality); //triggered events

      if(fUseCentrality) {
	//Centrality stuff
	AliCentrality *centrality = gESD->GetCentrality();

	fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
	
	// take only events inside centrality class
	if(!centrality->IsEventInCentralityClass(fCentralityPercentileMin,
						 fCentralityPercentileMax,
						 fCentralityEstimator.Data()))
	  return;
	
	// centrality QA (V0M)
	fHistV0M->Fill(gESD->GetVZEROData()->GetMTotV0A(), gESD->GetVZEROData()->GetMTotV0C());
      }
	
      const AliESDVertex *vertex = gESD->GetPrimaryVertex();
      if(vertex) {
	if(vertex->GetNContributors() > 0) {
	  if(vertex->GetZRes() != 0) {
	    fHistEventStats->Fill(3,fCentrality); //events with a proper vertex
	    if(TMath::Abs(vertex->GetXv()) < fVxMax) {
	      if(TMath::Abs(vertex->GetYv()) < fVyMax) {
		if(TMath::Abs(vertex->GetZv()) < fVzMax) {
		  fHistEventStats->Fill(4,fCentrality); //analayzed events
		  fHistVx->Fill(vertex->GetXv());
		  fHistVy->Fill(vertex->GetYv());
		  fHistVz->Fill(vertex->GetZv(),fCentrality);
		  
		  //========Get the VZERO event plane========//
		  Double_t gVZEROEventPlane = -10.0;
		  Double_t qxTot = 0.0, qyTot = 0.0;
		  AliEventplane *ep = gESD->GetEventplane();
		  if(ep) 
		    gVZEROEventPlane = ep->CalculateVZEROEventPlane(gESD,10,2,qxTot,qyTot);
		  if(gVZEROEventPlane < 0.) gVZEROEventPlane += TMath::Pi();
		  gReactionPlane = gVZEROEventPlane*TMath::RadToDeg();
		  fHistEventPlane->Fill(gReactionPlane,fCentrality);
		  //========Get the VZERO event plane========//

		  //Printf("There are %d tracks in this event", gESD->GetNumberOfTracks());
		  for (Int_t iTracks = 0; iTracks < gESD->GetNumberOfTracks(); iTracks++) {
		    AliESDtrack* track = dynamic_cast<AliESDtrack *>(gESD->GetTrack(iTracks));
		    if (!track) {
		      Printf("ERROR: Could not receive track %d", iTracks);
		      continue;
		    }	
		    
		    Int_t label = TMath::Abs(track->GetLabel());
		    if(label > mcEvent->GetNumberOfTracks()) continue;
		    if(label > mcEvent->GetNumberOfPrimaries()) continue;
		    
		    AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(label));
		    if(!mcTrack) continue;

		    // take only TPC only tracks
		    track_TPC   = new AliESDtrack();
		    if(!track->FillTPCOnlyTrack(*track_TPC)) continue;
		    
		    //ESD track cuts
		    if(fESDtrackCuts) 
		      if(!fESDtrackCuts->AcceptTrack(track_TPC)) continue;
		    
		    // fill QA histograms
		    Float_t b[2];
		    Float_t bCov[3];
		    track_TPC->GetImpactParameters(b,bCov);
		    if (bCov[0]<=0 || bCov[2]<=0) {
		      AliDebug(1, "Estimated b resolution lower or equal zero!");
		      bCov[0]=0; bCov[2]=0;
		    }
		    
		    Int_t nClustersTPC = -1;
		    nClustersTPC = track_TPC->GetTPCNclsIter1();   // TPC standalone
		    //nClustersTPC = track->GetTPCclusters(0);   // global track
		    Float_t chi2PerClusterTPC = -1;
		    if (nClustersTPC!=0) {
		      chi2PerClusterTPC = track_TPC->GetTPCchi2Iter1()/Float_t(nClustersTPC);      // TPC standalone
		      //chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);     // global track
		    }
		    
		    v_charge = track_TPC->Charge();
		    v_y      = track_TPC->Y();
		    v_eta    = track_TPC->Eta();
		    v_phi    = track_TPC->Phi() * TMath::RadToDeg();
		    v_E      = track_TPC->E();
		    v_pt     = track_TPC->Pt();
		    track_TPC->PxPyPz(v_p);

		    fHistClus->Fill(track_TPC->GetITSclusters(0),nClustersTPC);
		    fHistDCA->Fill(b[1],b[0]);
		    fHistChi2->Fill(chi2PerClusterTPC,fCentrality);
		    fHistPt->Fill(v_pt,fCentrality);
		    fHistEta->Fill(v_eta,fCentrality);
		    fHistPhi->Fill(v_phi,fCentrality);
		    fHistRapidity->Fill(v_y,fCentrality);
		    if(v_charge > 0) fHistPhiPos->Fill(v_phi,fCentrality);
		    else if(v_charge < 0) fHistPhiNeg->Fill(v_phi,fCentrality);

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

		    if(fRunShuffling) {
		      chargeVectorShuffle[0]->push_back(v_charge);
		      chargeVectorShuffle[1]->push_back(v_y);
		      chargeVectorShuffle[2]->push_back(v_eta);
		      chargeVectorShuffle[3]->push_back(v_phi);
		      chargeVectorShuffle[4]->push_back(v_p[0]);
		      chargeVectorShuffle[5]->push_back(v_p[1]);
		      chargeVectorShuffle[6]->push_back(v_p[2]);
		      chargeVectorShuffle[7]->push_back(v_pt);
		      chargeVectorShuffle[8]->push_back(v_E);
		    }
		    
		    delete track_TPC;
      		    gNumberOfAcceptedTracks += 1;
		    
		  } //track loop
		}//Vz cut
	      }//Vy cut
	    }//Vx cut
	  }//proper vertex resolution
	}//proper number of contributors
      }//vertex object valid
    }//triggered event 
  }//MC-ESD analysis

  //MC analysis
  else if(gAnalysisLevel == "MC") {
    AliMCEvent*  mcEvent = MCEvent(); 
    if (!mcEvent) {
      Printf("ERROR: mcEvent not available");
      return;
    }
    fHistEventStats->Fill(1,fCentrality); //total events
    fHistEventStats->Fill(2,fCentrality); //offline trigger

    Double_t gImpactParameter = 0.;
    if(fUseCentrality) {
      //Get the MC header
      AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
      if (headerH) {
	//Printf("=====================================================");
	//Printf("Reaction plane angle: %lf",headerH->ReactionPlaneAngle());
	//Printf("Impact parameter: %lf",headerH->ImpactParameter());
	//Printf("=====================================================");
	gReactionPlane = headerH->ReactionPlaneAngle();
	gImpactParameter = headerH->ImpactParameter();
	fCentrality = gImpactParameter;
      }
      fCentrality = gImpactParameter;
      fHistEventPlane->Fill(gReactionPlane*TMath::RadToDeg(),fCentrality);

      // take only events inside centrality class (DIDN'T CHANGE THIS UP TO NOW)
      if((fImpactParameterMin > gImpactParameter) || (fImpactParameterMax < gImpactParameter))
	return;
    }
    
    AliGenEventHeader *header = mcEvent->GenEventHeader();
    if(!header) return;
    
    TArrayF gVertexArray;
    header->PrimaryVertex(gVertexArray);
    //Printf("Vertex: %lf (x) - %lf (y) - %lf (z)",
    //gVertexArray.At(0),
    //gVertexArray.At(1),
    //gVertexArray.At(2));
    fHistEventStats->Fill(3,fCentrality); //events with a proper vertex
    if(TMath::Abs(gVertexArray.At(0)) < fVxMax) {
      if(TMath::Abs(gVertexArray.At(1)) < fVyMax) {
	if(TMath::Abs(gVertexArray.At(2)) < fVzMax) {
	  fHistEventStats->Fill(4,fCentrality); //analayzed events
	  fHistVx->Fill(gVertexArray.At(0));
	  fHistVy->Fill(gVertexArray.At(1));
	  fHistVz->Fill(gVertexArray.At(2),fCentrality);
	  
	  Printf("There are %d tracks in this event", mcEvent->GetNumberOfPrimaries());
	  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfPrimaries(); iTracks++) {
	    AliMCParticle* track = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iTracks));
	    if (!track) {
	      Printf("ERROR: Could not receive particle %d", iTracks);
	      continue;
	    }
	    
	    //exclude non stable particles
	    if(!mcEvent->IsPhysicalPrimary(iTracks)) continue;

	    v_eta    = track->Eta();
	    v_pt     = track->Pt();
	    v_y      = track->Y();

	    if( v_pt < fPtMin || v_pt > fPtMax)      
	      continue;
	    if (!fUsePID) {
	      if( v_eta < fEtaMin || v_eta > fEtaMax)  continue;
	    }
	    else if (fUsePID){
	      if( v_y < fEtaMin || v_y > fEtaMax)  continue;
	    }

	    //analyze one set of particles
	    if(fUseMCPdgCode) {
	      TParticle *particle = track->Particle();
	      if(!particle) continue;
	      
	      Int_t gPdgCode = particle->GetPdgCode();
	      if(TMath::Abs(fPDGCodeToBeAnalyzed) != TMath::Abs(gPdgCode)) 
		continue;
	    }
	    
	    //Use the acceptance parameterization
	    if(fAcceptanceParameterization) {
	      Double_t gRandomNumber = gRandom->Rndm();
	      if(gRandomNumber > fAcceptanceParameterization->Eval(track->Pt())) 
		continue;
	    }
	    
	    //Exclude resonances
	    if(fExcludeResonancesInMC) {
	      TParticle *particle = track->Particle();
	      if(!particle) continue;
	      
	      Bool_t kExcludeParticle = kFALSE;
	      Int_t gMotherIndex = particle->GetFirstMother();
	      if(gMotherIndex != -1) {
		AliMCParticle* motherTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(gMotherIndex));
		if(motherTrack) {
		  TParticle *motherParticle = motherTrack->Particle();
		  if(motherParticle) {
		    Int_t pdgCodeOfMother = motherParticle->GetPdgCode();
		    //if((pdgCodeOfMother == 113)||(pdgCodeOfMother == 213)||(pdgCodeOfMother == 221)||(pdgCodeOfMother == 223)||(pdgCodeOfMother == 331)||(pdgCodeOfMother == 333)) {
		    if(pdgCodeOfMother == 113) {
		      kExcludeParticle = kTRUE;
		    }
		  }
		}
	      }
	      
	      //Exclude from the analysis decay products of rho0, rho+, eta, eta' and phi
	      if(kExcludeParticle) continue;
	    }

	    v_charge = track->Charge();
	    v_phi    = track->Phi();
	    v_E      = track->E();
	    track->PxPyPz(v_p);
	    //Printf("phi (before): %lf",v_phi);

	    fHistPt->Fill(v_pt,fCentrality);
	    fHistEta->Fill(v_eta,fCentrality);
	    fHistPhi->Fill(v_phi*TMath::RadToDeg(),fCentrality);
	    fHistRapidity->Fill(v_y,fCentrality);
	    if(v_charge > 0) fHistPhiPos->Fill(v_phi*TMath::RadToDeg(),fCentrality);
	    else if(v_charge < 0) fHistPhiNeg->Fill(v_phi*TMath::RadToDeg(),fCentrality);

	    //Flow after burner
	    if(fUseFlowAfterBurner) {
	      Double_t precisionPhi = 0.001;
	      Int_t maxNumberOfIterations = 100;

	      Double_t phi0 = v_phi;
	      Double_t gV2 = fDifferentialV2->Eval(v_pt);

	      for (Int_t j = 0; j < maxNumberOfIterations; j++) {
		Double_t phiprev = v_phi;
		Double_t fl = v_phi - phi0 + gV2*TMath::Sin(2.*(v_phi - gReactionPlane));
		Double_t fp = 1.0 + 2.0*gV2*TMath::Cos(2.*(v_phi - gReactionPlane)); 
		v_phi -= fl/fp;
		if (TMath::AreEqualAbs(phiprev,v_phi,precisionPhi)) break;
	      }
	      //Printf("phi (after): %lf\n",v_phi);
	      	      Double_t v_DeltaphiBefore = phi0 - gReactionPlane;
	      if(v_DeltaphiBefore < 0) v_DeltaphiBefore += 2*TMath::Pi();
	      fHistPhiBefore->Fill(v_DeltaphiBefore,fCentrality);

	      Double_t v_DeltaphiAfter = v_phi - gReactionPlane;
	      if(v_DeltaphiAfter < 0) v_DeltaphiAfter += 2*TMath::Pi();
	      fHistPhiAfter->Fill(v_DeltaphiAfter,fCentrality);
	    }
	    
	    v_phi *= TMath::RadToDeg();

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
	    
	    if(fRunShuffling) {
	      chargeVectorShuffle[0]->push_back(v_charge);
	      chargeVectorShuffle[1]->push_back(v_y);
	      chargeVectorShuffle[2]->push_back(v_eta);
	      chargeVectorShuffle[3]->push_back(v_phi);
	      chargeVectorShuffle[4]->push_back(v_p[0]);
	      chargeVectorShuffle[5]->push_back(v_p[1]);
	      chargeVectorShuffle[6]->push_back(v_p[2]);
	      chargeVectorShuffle[7]->push_back(v_pt);
	      chargeVectorShuffle[8]->push_back(v_E);
	    }
	    gNumberOfAcceptedTracks += 1;
		    
	  } //track loop
	  gReactionPlane *= TMath::RadToDeg();
	}//Vz cut
      }//Vy cut
    }//Vx cut
  }//MC analysis
  
  //multiplicity cut (used in pp)
  if(fUseMultiplicity) {
    if((gNumberOfAcceptedTracks < fNumberOfAcceptedTracksMin)||(gNumberOfAcceptedTracks > fNumberOfAcceptedTracksMax))
      return;
  }
  fHistNumberOfAcceptedTracks->Fill(gNumberOfAcceptedTracks,fCentrality);
  
  // calculate balance function
  if(fUseMultiplicity) 
    fBalance->CalculateBalance(gNumberOfAcceptedTracks,gReactionPlane,chargeVector);
  else                 
    fBalance->CalculateBalance(fCentrality,gReactionPlane,chargeVector);

  if(fRunShuffling) {
    // shuffle charges
    random_shuffle( chargeVectorShuffle[0]->begin(), chargeVectorShuffle[0]->end() );
    if(fUseMultiplicity) 
      fShuffledBalance->CalculateBalance(gNumberOfAcceptedTracks,gReactionPlane,chargeVectorShuffle);
    else                 
      fShuffledBalance->CalculateBalance(fCentrality,gReactionPlane,chargeVectorShuffle);
  }
}      

//________________________________________________________________________
void  AliAnalysisTaskBFPsi::FinishTaskOutput(){
  //Printf("END BF");

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
void AliAnalysisTaskBFPsi::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}
