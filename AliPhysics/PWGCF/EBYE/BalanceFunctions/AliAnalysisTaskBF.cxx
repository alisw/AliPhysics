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
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"

#include "TH2D.h"                  
#include "AliPID.h"                
#include "AliPIDResponse.h"        
#include "AliPIDCombined.h"        

#include "AliAnalysisTaskBF.h"
#include "AliBalance.h"

#include <random>


// Analysis task for the BF code
// Authors: Panos.Christakoglou@nikhef.nl

ClassImp(AliAnalysisTaskBF)

//________________________________________________________________________
AliAnalysisTaskBF::AliAnalysisTaskBF(const char *name) 
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
  fAODtrackCutBit(128),
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
AliAnalysisTaskBF::~AliAnalysisTaskBF() {

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
void AliAnalysisTaskBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  if(!fBalance) {
    fBalance = new AliBalance();
    fBalance->SetAnalysisLevel("ESD");
    //fBalance->SetNumberOfBins(-1,16);
    fBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6);
  }
  if(fRunShuffling) {
    if(!fShuffledBalance) {
      fShuffledBalance = new AliBalance();
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

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  
  if ((gAnalysisLevel == "ESD") || (gAnalysisLevel == "AOD") || (gAnalysisLevel == "MCESD")) {
    fHistEventStats = new TH2D("fHistEventStats",
			       "Event statistics;;Centrality",
			       4,0.5,4.5, 100,0,100);
  }
  
  if (gAnalysisLevel == "MC"){
    fHistEventStats = new TH2D("fHistEventStats",
			       "Event statistics;;Centrality",
			       4,0.5,4.5, 10000,0,15);
  }


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

  fHistNumberOfAcceptedTracks = new TH2D("fHistNumberOfAcceptedTracks",";N_{acc.};;Centrality",4001,-0.5,4000.5,100,0,100);
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
  fHistRapidity  = new TH1F("fHistRapidity","y distribution",200,-2,2);
  fList->Add(fHistRapidity);
  fHistPhi  = new TH1F("fHistPhi","#phi distribution",200,-20,380);
  fList->Add(fHistPhi);
  fHistPhiBefore  = new TH1F("fHistPhiBefore","#phi distribution",200,0.,2*TMath::Pi());
  fList->Add(fHistPhiBefore);
  fHistPhiAfter  = new TH1F("fHistPhiAfter","#phi distribution",200,0.,2*TMath::Pi());
  fList->Add(fHistPhiAfter);
  fHistPhiPos  = new TH1F("fHistPhiPos","#phi distribution for positive particles",200,-20,380);
  fList->Add(fHistPhiPos);
  fHistPhiNeg  = new TH1F("fHistPhiNeg","#phi distribution for negative particles",200,-20,380);
  fList->Add(fHistPhiNeg);
  fHistV0M  = new TH2F("fHistV0M","V0 Multiplicity C vs. A",500, 0, 20000, 500, 0, 20000);
  fList->Add(fHistV0M);
  TString gRefTrackName[6] = {"tracks","tracksPos","tracksNeg","tracksTPConly","clusITS0","clusITS1"};
  fHistRefTracks  = new TH2F("fHistRefTracks","Nr of Ref tracks/event vs. ref track estimator;;Nr of tracks",6, 0, 6, 400, 0, 20000);
  for(Int_t i = 1; i <= 6; i++)
    fHistRefTracks->GetXaxis()->SetBinLabel(i,gRefTrackName[i-1].Data());
  fList->Add(fHistRefTracks);

  // QA histograms for HBTinspired and Conversion cuts
  fList->Add(fBalance->GetQAHistHBTbefore());
  fList->Add(fBalance->GetQAHistHBTafter());
  fList->Add(fBalance->GetQAHistConversionbefore());
  fList->Add(fBalance->GetQAHistConversionafter());

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
    fListBF->Add(fBalance->GetHistNp(a));
    fListBF->Add(fBalance->GetHistNn(a));
    fListBF->Add(fBalance->GetHistNpn(a));
    fListBF->Add(fBalance->GetHistNnn(a));
    fListBF->Add(fBalance->GetHistNpp(a));
    fListBF->Add(fBalance->GetHistNnp(a));

    if(fRunShuffling) {
      fListBFS->Add(fShuffledBalance->GetHistNp(a));
      fListBFS->Add(fShuffledBalance->GetHistNn(a));
      fListBFS->Add(fShuffledBalance->GetHistNpn(a));
      fListBFS->Add(fShuffledBalance->GetHistNnn(a));
      fListBFS->Add(fShuffledBalance->GetHistNpp(a));
      fListBFS->Add(fShuffledBalance->GetHistNnp(a));
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
  PostData(2, fListBF);
  if(fRunShuffling) PostData(3, fListBFS);
  if(fUsePID) PostData(4, fHistListPIDQA);       //PID

  TH1::AddDirectory(oldStatus);

}

//________________________________________________________________________
void AliAnalysisTaskBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  AliESDtrack *trackTPC   = NULL;

  Int_t gNumberOfAcceptedTracks = 0;
  Float_t fCentrality           = -999.;

  // for HBT like cuts need magnetic field sign
  Float_t bSign = 0; // only used in AOD so far

  // vector holding the charges/kinematics of all tracks (charge,y,eta,phi,p0,p1,p2,pt,E)
  vector<Double_t> *chargeVectorShuffle[9];   // this will be shuffled
  vector<Double_t> *chargeVector[9];          // original charge
  for(Int_t i = 0; i < 9; i++){
    chargeVectorShuffle[i] = new vector<Double_t>;
    chargeVector[i]        = new vector<Double_t>;
  }

  Double_t vCharge;
  Double_t vY;
  Double_t vEta;
  Double_t vPhi;
  Double_t vP[3];
  Double_t vPt;
  Double_t vE;

  if(fUsePID) {
    fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
  }
 
  //ESD analysis
  if(gAnalysisLevel == "ESD") {
    AliESDEvent* gESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from TaskSE
    if (!gESD) {
      AliError("ERROR: gESD not available");
      return;
    }

    // store offline trigger bits
    fHistTriggerStats->Fill(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());

    AliCentrality *centrality = 0x0; 
    if(fUseCentrality) {
      //Centrality stuff
      centrality = gESD->GetCentrality();
      fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
    }

    // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
    fHistEventStats->Fill(1,fCentrality); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2,fCentrality); //triggered events

      if(fUseCentrality) {
	//Centrality stuff
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
	    if(TMath::Abs(vertex->GetX()) < fVxMax) {
	      if(TMath::Abs(vertex->GetY()) < fVyMax) {
		if(TMath::Abs(vertex->GetZ()) < fVzMax) {
		  fHistEventStats->Fill(4,fCentrality); //analayzed events
		  fHistVx->Fill(vertex->GetX());
		  fHistVy->Fill(vertex->GetY());
		  fHistVz->Fill(vertex->GetZ());
		  
		  //Printf("There are %d tracks in this event", gESD->GetNumberOfTracks());
		  for (Int_t iTracks = 0; iTracks < gESD->GetNumberOfTracks(); iTracks++) {
		    AliESDtrack* track = dynamic_cast<AliESDtrack *>(gESD->GetTrack(iTracks));
		    if (!track) {
		      AliError(Form("ERROR: Could not receive track %d", iTracks));
		      continue;
		    }	
		    
		    // take only TPC only tracks
		    trackTPC   = new AliESDtrack();
		    if(!track->FillTPCOnlyTrack(*trackTPC)) continue;
		    
		    //ESD track cuts
		    if(fESDtrackCuts) 
		      if(!fESDtrackCuts->AcceptTrack(trackTPC)) continue;
		    
		    // fill QA histograms
		    Float_t b[2];
		    Float_t bCov[3];
		    trackTPC->GetImpactParameters(b,bCov);
		    if (bCov[0]<=0 || bCov[2]<=0) {
		      AliDebug(1, "Estimated b resolution lower or equal zero!");
		      bCov[0]=0; bCov[2]=0;
		    }
		    
		    Int_t nClustersTPC = -1;
		    nClustersTPC = trackTPC->GetTPCNclsIter1();   // TPC standalone
		    //nClustersTPC = track->GetTPCclusters(0);   // global track
		    Float_t chi2PerClusterTPC = -1;
		    if (nClustersTPC!=0) {
		      chi2PerClusterTPC = trackTPC->GetTPCchi2Iter1()/Float_t(nClustersTPC);      // TPC standalone
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
		    vCharge = trackTPC->Charge();
		    vY      = trackTPC->Y();
		    vEta    = trackTPC->Eta();
		    vPhi    = trackTPC->Phi() * TMath::RadToDeg();
		    vE      = trackTPC->E();
		    vPt     = trackTPC->Pt();
		    trackTPC->PxPyPz(vP);
		    fHistClus->Fill(trackTPC->GetITSclusters(0),nClustersTPC);
		    fHistDCA->Fill(b[1],b[0]);
		    fHistChi2->Fill(chi2PerClusterTPC);
		    fHistPt->Fill(vPt);
		    fHistEta->Fill(vEta);
		    fHistPhi->Fill(vPhi);
		    fHistRapidity->Fill(vY);
		    if(vCharge > 0) fHistPhiPos->Fill(vPhi);
		    else if(vCharge < 0) fHistPhiNeg->Fill(vPhi);

		    // fill charge vector
		    chargeVector[0]->push_back(vCharge);
		    chargeVector[1]->push_back(vY);
		    chargeVector[2]->push_back(vEta);
		    chargeVector[3]->push_back(vPhi);
		    chargeVector[4]->push_back(vP[0]);
		    chargeVector[5]->push_back(vP[1]);
		    chargeVector[6]->push_back(vP[2]);
		    chargeVector[7]->push_back(vPt);
		    chargeVector[8]->push_back(vE);

		    if(fRunShuffling) {
		      chargeVectorShuffle[0]->push_back(vCharge);
		      chargeVectorShuffle[1]->push_back(vY);
		      chargeVectorShuffle[2]->push_back(vEta);
		      chargeVectorShuffle[3]->push_back(vPhi);
		      chargeVectorShuffle[4]->push_back(vP[0]);
		      chargeVectorShuffle[5]->push_back(vP[1]);
		      chargeVectorShuffle[6]->push_back(vP[2]);
		      chargeVectorShuffle[7]->push_back(vPt);
		      chargeVectorShuffle[8]->push_back(vE);
		    }
		    
		    delete trackTPC;
		    gNumberOfAcceptedTracks += 1;
		  } //track loop
		  // cout<<"Centrality: "<<fCentrality<<" - Accepted tracks: "<<gNumberOfAcceptedTracks<<endl;
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
      AliError("ERROR: gAOD not available");
      return;
    }

    // for HBT like cuts need magnetic field sign
    bSign = (gAOD->GetMagneticField() > 0) ? 1 : -1;

    AliAODHeader *aodHeader = dynamic_cast<AliAODHeader*>(gAOD->GetHeader());
    if(!aodHeader) AliFatal("Not a standard AOD");

    // store offline trigger bits
    fHistTriggerStats->Fill(aodHeader->GetOfflineTrigger());

    if(fUseCentrality) {
      fCentrality = aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
    }
    
    //event selection done in AliAnalysisTaskSE::Exec() --> this is not used
    fHistEventStats->Fill(1,fCentrality); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2,fCentrality); //triggered events
		  
      //Centrality stuff (centrality in AOD header)
      if(fUseCentrality) {
    	//fCentrality = aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
	// in OLD AODs (i.e. AOD049) fCentrality can be == 0
	if(fCentrality == 0) 
	  return;
 	
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
      		  fHistVz->Fill(vertex->GetZ());
		  
		  //===========================================//
		  TExMap *trackMap = new TExMap();
		  for (Int_t iTracks = 0; iTracks < gAOD->GetNumberOfTracks(); iTracks++) {
		    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(gAOD->GetTrack(iTracks));
		    if (!aodTrack) {
		      AliError(Form("ERROR: Could not receive track %d", iTracks));
		      continue;
		    }
		    Int_t gID = aodTrack->GetID();
		    //if (!aodTrack->TestFilterBit(fAODtrackCutBit)) trackMap->Add(gID, iTracks);
		    if (aodTrack->TestFilterBit(1)) trackMap->Add(gID, iTracks);
		  }
		  AliAODTrack* newAodTrack; 
		  //===========================================//

      		  //Printf("There are %d tracks in this event", gAOD->GetNumberOfTracks());
      		  for (Int_t iTracks = 0; iTracks < gAOD->GetNumberOfTracks(); iTracks++) {
      		    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(gAOD->GetTrack(iTracks));
      		    if (!aodTrack) {
      		      AliError(Form("ERROR: Could not receive track %d", iTracks));
      		      continue;
      		    }
		    
      		    // AOD track cuts		    
      		    // For ESD Filter Information: ANALYSIS/macros/AddTaskESDfilter.C
		    //===========================================//
		    // take only TPC only tracks 
		    fHistTrackStats->Fill(aodTrack->GetFilterMap());
		    if(!aodTrack->TestFilterBit(fAODtrackCutBit)) continue;

		    Int_t gID = aodTrack->GetID();
		    newAodTrack = gID >= 0 ? aodTrack : dynamic_cast<AliAODTrack*>(gAOD->GetTrack(trackMap->GetValue(-1-gID)));
                    if(!newAodTrack) AliFatal("Not a standard AOD");
		    //Printf("Label: %d - Pt: %lf (old) - %d - Pt: %lf(new)",gID,aodTrack->Pt(), newAodTrack->GetID(), newAodTrack->Pt());
                    //===========================================//

		    //fHistTrackStats->Fill(aodTrack->GetFilterMap());
      		    //if(!aodTrack->TestFilterBit(fAODtrackCutBit)) continue;
		    
      		    vCharge = aodTrack->Charge();
      		    vY      = aodTrack->Y();
      		    vEta    = aodTrack->Eta();
      		    vPhi    = aodTrack->Phi() * TMath::RadToDeg();
      		    vE      = aodTrack->E();
      		    vPt     = aodTrack->Pt();
      		    aodTrack->PxPyPz(vP);
		    
      		    Float_t dcaXY = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
      		    Float_t dcaZ  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
		    
		    
      		    // Kinematics cuts from ESD track cuts
      		    if( vPt < fPtMin || vPt > fPtMax)      continue;

		    if (!fUsePID) {
		      if( vEta < fEtaMin || vEta > fEtaMax)  continue;
		    }

		    else if (fUsePID){
		      if( vY < fEtaMin || vY > fEtaMax)  continue;
		    }

      		    // Extra DCA cuts (for systematic studies [!= -1])
      		    if( fDCAxyCut != -1 && fDCAzCut != -1){
      		      if(TMath::Sqrt((dcaXY*dcaXY)/(fDCAxyCut*fDCAxyCut)+(dcaZ*dcaZ)/(fDCAzCut*fDCAzCut)) > 1 ){
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

		   //===============================================PID==================================//		    		   
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
			nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(newAodTrack,(AliPID::EParticleType)fParticleOfInterest));
			//detUsedTPC = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPC);
                        detUsedTPC = (AliPIDResponse::kDetTPC);
			for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
			  prob[iSpecies] = probTPC[iSpecies];
			break;
		      case kTOFpid:
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
			nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(newAodTrack,(AliPID::EParticleType)fParticleOfInterest));
			//detUsedTOF = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTOF);
                        detUsedTPC = (AliPIDResponse::kDetTPC);
			for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
			  prob[iSpecies] = probTOF[iSpecies];
			break;
		      case kTPCTOF:
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
			//detUsedTPCTOF = fPIDCombined->ComputeProbabilities(newAodTrack, fPIDResponse, probTPCTOF);
                        detUsedTPC = (AliPIDResponse::kDetTPC);
			for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
			  prob[iSpecies] = probTPCTOF[iSpecies];
			break;
		      default:
			break;
		      }//end switch: define detector mask
		      
		      //Filling the PID QA
		      Double_t tofTime = -999., tof = -999.; //length = 999., tof = -999.;
		      //Double_t c = TMath::C()*1.E-9;// m/ns
		      //Double_t beta = -999.;
		      Double_t  nSigmaTOFForParticleOfInterest = -999.;
		      if ( (newAodTrack->IsOn(AliAODTrack::kTOFin)) &&
			   (newAodTrack->IsOn(AliAODTrack::kTIME))  ) { 
			tofTime = newAodTrack->GetTOFsignal();//in ps
			//length = newAodTrack->GetIntegratedLength();
			tof = tofTime*1E-3; // ns	
			
			if (tof <= 0) {
			  //Printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
			  continue;
			}
			//if (length <= 0){
			  //printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
			  //continue;
			//}
			
			//length = length*0.01; // in meters
			//tof = tof*c;
			//beta = length/tof;
			
			nSigmaTOFForParticleOfInterest = fPIDResponse->NumberOfSigmasTOF(newAodTrack,(AliPID::EParticleType)fParticleOfInterest);
			//fHistBetavsPTOFbeforePID ->Fill(aodTrack->P()*aodTrack->Charge(),beta);
			fHistProbTOFvsPtbeforePID ->Fill(newAodTrack->Pt(),probTOF[fParticleOfInterest]);
			fHistNSigmaTOFvsPtbeforePID ->Fill(newAodTrack->Pt(),nSigmaTOFForParticleOfInterest);
		      }//TOF signal 
		      
		      
		      Double_t  nSigmaTPCForParticleOfInterest = fPIDResponse->NumberOfSigmasTPC(newAodTrack,(AliPID::EParticleType)fParticleOfInterest);
		      fHistdEdxVsPTPCbeforePID -> Fill(newAodTrack->P()*newAodTrack->Charge(),newAodTrack->GetTPCsignal());
		      fHistProbTPCvsPtbeforePID -> Fill(newAodTrack->Pt(),probTPC[fParticleOfInterest]); 
		      fHistNSigmaTPCvsPtbeforePID -> Fill(newAodTrack->Pt(),nSigmaTPCForParticleOfInterest); 
		      fHistProbTPCTOFvsPtbeforePID -> Fill(newAodTrack->Pt(),probTPCTOF[fParticleOfInterest]);
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
			//fHistBetavsPTOFafterPID ->Fill(newAodTrack->P()*newAodTrack->Charge(),beta);
			fHistProbTOFvsPtafterPID ->Fill(newAodTrack->Pt(),probTOF[fParticleOfInterest]);
			fHistNSigmaTOFvsPtafterPID ->Fill(newAodTrack->Pt(),nSigmaTOFForParticleOfInterest);
			
			fHistdEdxVsPTPCafterPID -> Fill(newAodTrack->P()*newAodTrack->Charge(),newAodTrack->GetTPCsignal());
			fHistProbTPCvsPtafterPID -> Fill(newAodTrack->Pt(),probTPC[fParticleOfInterest]); 
			fHistProbTPCTOFvsPtafterPID -> Fill(newAodTrack->Pt(),probTPCTOF[fParticleOfInterest]);
			fHistNSigmaTPCvsPtafterPID -> Fill(newAodTrack->Pt(),nSigmaTPCForParticleOfInterest); 
		      }
		      
		      PostData(4, fHistListPIDQA);
		    }

		    //=========================================================PID=================================================================//
		    		    
      		    // fill QA histograms
      		    fHistClus->Fill(aodTrack->GetITSNcls(),aodTrack->GetTPCNcls());
      		    fHistDCA->Fill(dcaZ,dcaXY);
      		    fHistChi2->Fill(aodTrack->Chi2perNDF());
      		    fHistPt->Fill(vPt);
      		    fHistEta->Fill(vEta);
      		    fHistPhi->Fill(vPhi);
		    fHistRapidity->Fill(vY);
		    if(vCharge > 0) fHistPhiPos->Fill(vPhi);
		    else if(vCharge < 0) fHistPhiNeg->Fill(vPhi);

      		    // fill charge vector
      		    chargeVector[0]->push_back(vCharge);
      		    chargeVector[1]->push_back(vY);
      		    chargeVector[2]->push_back(vEta);
      		    chargeVector[3]->push_back(vPhi);
      		    chargeVector[4]->push_back(vP[0]);
      		    chargeVector[5]->push_back(vP[1]);
      		    chargeVector[6]->push_back(vP[2]);
      		    chargeVector[7]->push_back(vPt);
      		    chargeVector[8]->push_back(vE);

      		    if(fRunShuffling) {
      		      chargeVectorShuffle[0]->push_back(vCharge);
      		      chargeVectorShuffle[1]->push_back(vY);
      		      chargeVectorShuffle[2]->push_back(vEta);
      		      chargeVectorShuffle[3]->push_back(vPhi);
      		      chargeVectorShuffle[4]->push_back(vP[0]);
      		      chargeVectorShuffle[5]->push_back(vP[1]);
      		      chargeVectorShuffle[6]->push_back(vP[2]);
      		      chargeVectorShuffle[7]->push_back(vPt);
      		      chargeVectorShuffle[8]->push_back(vE);
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
      AliError("ERROR: mcEvent not available");
      return;
    }

    AliESDEvent* gESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from TaskSE
    if (!gESD) {
      AliError("ERROR: gESD not available");
      return;
    }

    // store offline trigger bits
    fHistTriggerStats->Fill(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());

    AliCentrality *centrality = 0x0; 
    if(fUseCentrality) {
	centrality = gESD->GetCentrality();
	fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
    }

    // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
    fHistEventStats->Fill(1,fCentrality); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2,fCentrality); //triggered events

      if(fUseCentrality) {
	//Centrality stuff
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
	    if(TMath::Abs(vertex->GetX()) < fVxMax) {
	      if(TMath::Abs(vertex->GetY()) < fVyMax) {
		if(TMath::Abs(vertex->GetZ()) < fVzMax) {
		  fHistEventStats->Fill(4,fCentrality); //analayzed events
		  fHistVx->Fill(vertex->GetX());
		  fHistVy->Fill(vertex->GetY());
		  fHistVz->Fill(vertex->GetZ());
		  
		  //Printf("There are %d tracks in this event", gESD->GetNumberOfTracks());
		  for (Int_t iTracks = 0; iTracks < gESD->GetNumberOfTracks(); iTracks++) {
		    AliESDtrack* track = dynamic_cast<AliESDtrack *>(gESD->GetTrack(iTracks));
		    if (!track) {
		      AliError(Form("ERROR: Could not receive track %d", iTracks));
		      continue;
		    }	
		    
		    Int_t label = TMath::Abs(track->GetLabel());
		    if(label > mcEvent->GetNumberOfTracks()) continue;
		    if(label > mcEvent->GetNumberOfPrimaries()) continue;
		    
		    AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(label));
		    if(!mcTrack) continue;

		    // take only TPC only tracks
		    trackTPC   = new AliESDtrack();
		    if(!track->FillTPCOnlyTrack(*trackTPC)) continue;
		    
		    //ESD track cuts
		    if(fESDtrackCuts) 
		      if(!fESDtrackCuts->AcceptTrack(trackTPC)) continue;
		    
		    // fill QA histograms
		    Float_t b[2];
		    Float_t bCov[3];
		    trackTPC->GetImpactParameters(b,bCov);
		    if (bCov[0]<=0 || bCov[2]<=0) {
		      AliDebug(1, "Estimated b resolution lower or equal zero!");
		      bCov[0]=0; bCov[2]=0;
		    }
		    
		    Int_t nClustersTPC = -1;
		    nClustersTPC = trackTPC->GetTPCNclsIter1();   // TPC standalone
		    //nClustersTPC = track->GetTPCclusters(0);   // global track
		    Float_t chi2PerClusterTPC = -1;
		    if (nClustersTPC!=0) {
		      chi2PerClusterTPC = trackTPC->GetTPCchi2Iter1()/Float_t(nClustersTPC);      // TPC standalone
		      //chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);     // global track
		    }
		    
		    vCharge = trackTPC->Charge();
		    vY      = trackTPC->Y();
		    vEta    = trackTPC->Eta();
		    vPhi    = trackTPC->Phi() * TMath::RadToDeg();
		    vE      = trackTPC->E();
		    vPt     = trackTPC->Pt();
		    trackTPC->PxPyPz(vP);

		    fHistClus->Fill(trackTPC->GetITSclusters(0),nClustersTPC);
		    fHistDCA->Fill(b[1],b[0]);
		    fHistChi2->Fill(chi2PerClusterTPC);
		    fHistPt->Fill(vPt);
		    fHistEta->Fill(vEta);
		    fHistPhi->Fill(vPhi);
		    fHistRapidity->Fill(vY);
		    if(vCharge > 0) fHistPhiPos->Fill(vPhi);
		    else if(vCharge < 0) fHistPhiNeg->Fill(vPhi);

		    // fill charge vector
		    chargeVector[0]->push_back(vCharge);
		    chargeVector[1]->push_back(vY);
		    chargeVector[2]->push_back(vEta);
		    chargeVector[3]->push_back(vPhi);
		    chargeVector[4]->push_back(vP[0]);
		    chargeVector[5]->push_back(vP[1]);
		    chargeVector[6]->push_back(vP[2]);
		    chargeVector[7]->push_back(vPt);
		    chargeVector[8]->push_back(vE);

		    if(fRunShuffling) {
		      chargeVectorShuffle[0]->push_back(vCharge);
		      chargeVectorShuffle[1]->push_back(vY);
		      chargeVectorShuffle[2]->push_back(vEta);
		      chargeVectorShuffle[3]->push_back(vPhi);
		      chargeVectorShuffle[4]->push_back(vP[0]);
		      chargeVectorShuffle[5]->push_back(vP[1]);
		      chargeVectorShuffle[6]->push_back(vP[2]);
		      chargeVectorShuffle[7]->push_back(vPt);
		      chargeVectorShuffle[8]->push_back(vE);
		    }
		    
		    delete trackTPC;
      		    gNumberOfAcceptedTracks += 1;
		    
		  } //track loop
		  //cout<<"Centrality: "<<fCentrality<<" - Accepted tracks: "<<gNumberOfAcceptedTracks<<endl;
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
      AliError("ERROR: mcEvent not available");
      return;
    }

    //fHistEventStats->Fill(1,fCentrality); //total events
    //fHistEventStats->Fill(2,fCentrality); //offline trigger

    Double_t gReactionPlane = 0., gImpactParameter = 0.;
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

      // take only events inside centrality class (DIDN'T CHANGE THIS UP TO NOW)
      if((fImpactParameterMin > gImpactParameter) || (fImpactParameterMax < gImpactParameter))
	return;
    }

    fHistEventStats->Fill(1,fCentrality); //total events
    fHistEventStats->Fill(2,fCentrality); //offline trigger
    
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
	  fHistVz->Fill(gVertexArray.At(2));
	  
	  AliInfo(Form("There are %d tracks in this event", mcEvent->GetNumberOfPrimaries()));
	  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfPrimaries(); iTracks++) {
	    AliMCParticle* track = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iTracks));
	    if (!track) {
	      AliError(Form("ERROR: Could not receive particle %d", iTracks));
	      continue;
	    }
	    
	    //exclude non stable particles
	    if(!mcEvent->IsPhysicalPrimary(iTracks)) continue;

	    vEta    = track->Eta();
	    vPt     = track->Pt();
	    vY      = track->Y();

	    if( vPt < fPtMin || vPt > fPtMax)      
	      continue;
	    if (!fUsePID) {
	      if( vEta < fEtaMin || vEta > fEtaMax)  continue;
	    }
	    else if (fUsePID){
	      if( vY < fEtaMin || vY > fEtaMax)  continue;
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

	    vCharge = track->Charge();
	    vPhi    = track->Phi();
	    vE      = track->E();
	    track->PxPyPz(vP);
	    //Printf("phi (before): %lf",vPhi);

	    fHistPt->Fill(vPt);
	    fHistEta->Fill(vEta);
	    fHistPhi->Fill(vPhi);
	    fHistRapidity->Fill(vY);
	    if(vCharge > 0) fHistPhiPos->Fill(vPhi);
	    else if(vCharge < 0) fHistPhiNeg->Fill(vPhi);

	    //Flow after burner
	    if(fUseFlowAfterBurner) {
	      Double_t precisionPhi = 0.001;
	      Int_t maxNumberOfIterations = 100;

	      Double_t phi0 = vPhi;
	      Double_t gV2 = fDifferentialV2->Eval(vPt);

	      for (Int_t j = 0; j < maxNumberOfIterations; j++) {
		Double_t phiprev = vPhi;
		Double_t fl = vPhi - phi0 + gV2*TMath::Sin(2.*(vPhi - gReactionPlane));
		Double_t fp = 1.0 + 2.0*gV2*TMath::Cos(2.*(vPhi - gReactionPlane)); 
		vPhi -= fl/fp;
		if (TMath::AreEqualAbs(phiprev,vPhi,precisionPhi)) break;
	      }
	      //Printf("phi (after): %lf\n",vPhi);
	      Double_t vDeltaphiBefore = phi0 - gReactionPlane;
	      if(vDeltaphiBefore < 0) vDeltaphiBefore += 2*TMath::Pi();
	      fHistPhiBefore->Fill(vDeltaphiBefore);

	      Double_t vDeltaphiAfter = vPhi - gReactionPlane;
	      if(vDeltaphiAfter < 0) vDeltaphiAfter += 2*TMath::Pi();
	      fHistPhiAfter->Fill(vDeltaphiAfter);
	    }
	    
	    vPhi *= TMath::RadToDeg();

	    // fill charge vector
	    chargeVector[0]->push_back(vCharge);
	    chargeVector[1]->push_back(vY);
	    chargeVector[2]->push_back(vEta);
	    chargeVector[3]->push_back(vPhi);
	    chargeVector[4]->push_back(vP[0]);
	    chargeVector[5]->push_back(vP[1]);
	    chargeVector[6]->push_back(vP[2]);
	    chargeVector[7]->push_back(vPt);
	    chargeVector[8]->push_back(vE);
	    
	    if(fRunShuffling) {
	      chargeVectorShuffle[0]->push_back(vCharge);
	      chargeVectorShuffle[1]->push_back(vY);
	      chargeVectorShuffle[2]->push_back(vEta);
	      chargeVectorShuffle[3]->push_back(vPhi);
	      chargeVectorShuffle[4]->push_back(vP[0]);
	      chargeVectorShuffle[5]->push_back(vP[1]);
	      chargeVectorShuffle[6]->push_back(vP[2]);
	      chargeVectorShuffle[7]->push_back(vPt);
	      chargeVectorShuffle[8]->push_back(vE);
	    }
	    gNumberOfAcceptedTracks += 1;
		    
	  } //track loop
	}//Vz cut
      }//Vy cut
    }//Vx cut
  }//MC analysis
  
  //multiplicity cut (used in pp)
  if(fUseMultiplicity) {
    if((gNumberOfAcceptedTracks < fNumberOfAcceptedTracksMin)||(gNumberOfAcceptedTracks > fNumberOfAcceptedTracksMax))
      return;
  }
  fHistNumberOfAcceptedTracks->Fill(gNumberOfAcceptedTracks, fCentrality);

  // calculate balance function
  if(fUseMultiplicity) 
    fBalance->CalculateBalance(gNumberOfAcceptedTracks,chargeVector,bSign);
  else                 
    fBalance->CalculateBalance(fCentrality,chargeVector,bSign);

  if(fRunShuffling) {
    // shuffle charges
    std::random_device rd;
    std::default_random_engine engine{rd()};
    std::shuffle( chargeVectorShuffle[0]->begin(), chargeVectorShuffle[0]->end(), engine );
    if(fUseMultiplicity) 
      fShuffledBalance->CalculateBalance(gNumberOfAcceptedTracks,chargeVectorShuffle,bSign);
    else                 
      fShuffledBalance->CalculateBalance(fCentrality,chargeVectorShuffle,bSign);
  }
}      

//________________________________________________________________________
void  AliAnalysisTaskBF::FinishTaskOutput(){
  //Printf("END BF");

  if (!fBalance) {
    AliError("ERROR: fBalance not available");
    return;
  }  
  if(fRunShuffling) {
    if (!fShuffledBalance) {
      AliError("ERROR: fShuffledBalance not available");
      return;
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}
