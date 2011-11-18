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

#include "AliAnalysisTaskBF.h"
#include "AliBalance.h"


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
  fHistV0M(0),
  fHistRefTracks(0),
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

  // Post output data.
  PostData(1, fList);
  PostData(2, fListBF);
  if(fRunShuffling) PostData(3, fListBFS);
}

//________________________________________________________________________
void AliAnalysisTaskBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  AliESDtrack *track_TPC   = NULL;

  Int_t gNumberOfAcceptedTracks = 0;
  Float_t fCentrality           = 0.;

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
    fHistEventStats->Fill(1); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2); //triggered events

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
	    fHistEventStats->Fill(3); //events with a proper vertex
	    if(TMath::Abs(vertex->GetXv()) < fVxMax) {
	      if(TMath::Abs(vertex->GetYv()) < fVyMax) {
		if(TMath::Abs(vertex->GetZv()) < fVzMax) {
		  fHistEventStats->Fill(4); //analayzed events
		  fHistVx->Fill(vertex->GetXv());
		  fHistVy->Fill(vertex->GetYv());
		  fHistVz->Fill(vertex->GetZv());
		  
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

		    v_charge = track_TPC->Charge();
		    v_y      = track_TPC->Y();
		    v_eta    = track_TPC->Eta();
		    v_phi    = track_TPC->Phi() * TMath::RadToDeg();
		    v_E      = track_TPC->E();
		    v_pt     = track_TPC->Pt();
		    track_TPC->PxPyPz(v_p);

		    fHistClus->Fill(track_TPC->GetITSclusters(0),nClustersTPC);
		    fHistDCA->Fill(b[1],b[0]);
		    fHistChi2->Fill(chi2PerClusterTPC);
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
    fHistEventStats->Fill(1); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2); //triggered events
		  
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
      	    fHistEventStats->Fill(3); //events with a proper vertex
      	    if(TMath::Abs(vertex->GetX()) < fVxMax) {
      	      if(TMath::Abs(vertex->GetY()) < fVyMax) {
      		if(TMath::Abs(vertex->GetZ()) < fVzMax) {
      		  fHistEventStats->Fill(4); //analyzed events
      		  fHistVx->Fill(vertex->GetX());
      		  fHistVy->Fill(vertex->GetY());
      		  fHistVz->Fill(vertex->GetZ());
		  
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
      		    if( v_eta < fEtaMin || v_eta > fEtaMax)  continue;
		    
      		    // Extra DCA cuts (for systematic studies [!= -1])
      		    if( fDCAxyCut != -1 && fDCAxyCut != -1){
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
      		    fHistChi2->Fill(aodTrack->Chi2perNDF());
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
    fHistEventStats->Fill(1); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
      isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
      fHistEventStats->Fill(2); //triggered events

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
	    fHistEventStats->Fill(3); //events with a proper vertex
	    if(TMath::Abs(vertex->GetXv()) < fVxMax) {
	      if(TMath::Abs(vertex->GetYv()) < fVyMax) {
		if(TMath::Abs(vertex->GetZv()) < fVzMax) {
		  fHistEventStats->Fill(4); //analayzed events
		  fHistVx->Fill(vertex->GetXv());
		  fHistVy->Fill(vertex->GetYv());
		  fHistVz->Fill(vertex->GetZv());
		  
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
		    fHistChi2->Fill(chi2PerClusterTPC);
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
  }//MC-ESD analysis

  //MC analysis
  else if(gAnalysisLevel == "MC") {
    AliMCEvent*  mcEvent = MCEvent(); 
    if (!mcEvent) {
      Printf("ERROR: mcEvent not available");
      return;
    }
    fHistEventStats->Fill(1); //total events
    fHistEventStats->Fill(2); //offline trigger

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
      }
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
    fHistEventStats->Fill(3); //events with a proper vertex
    if(TMath::Abs(gVertexArray.At(0)) < fVxMax) {
      if(TMath::Abs(gVertexArray.At(1)) < fVyMax) {
	if(TMath::Abs(gVertexArray.At(2)) < fVzMax) {
	  fHistEventStats->Fill(4); //analayzed events
	  fHistVx->Fill(gVertexArray.At(0));
	  fHistVy->Fill(gVertexArray.At(1));
	  fHistVz->Fill(gVertexArray.At(2));
	  
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
	    
	    if( v_pt < fPtMin || v_pt > fPtMax)      
	      continue;
	    if( v_eta < fEtaMin || v_eta > fEtaMax)  
	      continue;
	    
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
		    if((pdgCodeOfMother == 113)||(pdgCodeOfMother == 213)||(pdgCodeOfMother == 221)||(pdgCodeOfMother == 223)||(pdgCodeOfMother == 331)||(pdgCodeOfMother == 333)) {
		      kExcludeParticle = kTRUE;
		    }
		  }
		}
	      }
	      
	      //Exclude from the analysis decay products of rho0, rho+, eta, eta' and phi
	      if(kExcludeParticle) continue;
	    }

	    v_charge = track->Charge();
	    v_y      = track->Y();
	    v_phi    = track->Phi() * TMath::RadToDeg();
	    v_E      = track->E();
	    track->PxPyPz(v_p);

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
  fHistNumberOfAcceptedTracks->Fill(gNumberOfAcceptedTracks);
  
  // calculate balance function
  fBalance->CalculateBalance(fCentrality,chargeVector);

  if(fRunShuffling) {
    // shuffle charges
    random_shuffle( chargeVectorShuffle[0]->begin(), chargeVectorShuffle[0]->end() );
    fShuffledBalance->CalculateBalance(fCentrality,chargeVectorShuffle);
  }
}      

//________________________________________________________________________
void  AliAnalysisTaskBF::FinishTaskOutput(){
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
void AliAnalysisTaskBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}
