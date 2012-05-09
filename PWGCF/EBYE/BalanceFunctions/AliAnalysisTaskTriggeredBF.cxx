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

#include "AliLog.h"

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

#include "TH2D.h"    
#include "AliTHn.h"              

#include "AliAnalysisTaskTriggeredBF.h"
#include "AliBalanceTriggered.h"


// Analysis task for the TriggeredBF code
// Authors: Panos.Christakoglou@nikhef.nl, m.weber@cern.ch

ClassImp(AliAnalysisTaskTriggeredBF)

//________________________________________________________________________
AliAnalysisTaskTriggeredBF::AliAnalysisTaskTriggeredBF(const char *name) 
: AliAnalysisTaskSE(name), 
  fBalance(0),
  fRunShuffling(kFALSE),
  fShuffledBalance(0),
  fList(0),
  fListTriggeredBF(0),
  fListTriggeredBFS(0),
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
  fNClustersTPCCut(-1)
{
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
AliAnalysisTaskTriggeredBF::~AliAnalysisTaskTriggeredBF() {

  // Destructor

}

//________________________________________________________________________
void AliAnalysisTaskTriggeredBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once
  if(!fBalance) {
    fBalance = new AliBalanceTriggered();
    fBalance->SetAnalysisLevel("AOD");
  }
  if(fRunShuffling) {
    if(!fShuffledBalance) {
      fShuffledBalance = new AliBalanceTriggered();
      fShuffledBalance->SetAnalysisLevel("AOD");
    }
  }

  //QA list
  fList = new TList();
  fList->SetName("listQA");
  fList->SetOwner();

  //Balance Function list
  fListTriggeredBF = new TList();
  fListTriggeredBF->SetName("listTriggeredBF");
  fListTriggeredBF->SetOwner();

  if(fRunShuffling) {
    fListTriggeredBFS = new TList();
    fListTriggeredBFS->SetName("listTriggeredBFShuffled");
    fListTriggeredBFS->SetOwner();
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

  fListTriggeredBF->Add(fBalance->GetHistNp());
  fListTriggeredBF->Add(fBalance->GetHistNn());
  fListTriggeredBF->Add(fBalance->GetHistNpn());
  fListTriggeredBF->Add(fBalance->GetHistNnn());
  fListTriggeredBF->Add(fBalance->GetHistNpp());
  fListTriggeredBF->Add(fBalance->GetHistNnp());
  
  if(fRunShuffling) {
    fListTriggeredBFS->Add(fShuffledBalance->GetHistNp());
    fListTriggeredBFS->Add(fShuffledBalance->GetHistNn());
    fListTriggeredBFS->Add(fShuffledBalance->GetHistNpn());
    fListTriggeredBFS->Add(fShuffledBalance->GetHistNnn());
    fListTriggeredBFS->Add(fShuffledBalance->GetHistNpp());
    fListTriggeredBFS->Add(fShuffledBalance->GetHistNnp());
  }  


 

  // Post output data.
  PostData(1, fList);
  PostData(2, fListTriggeredBF);
  if(fRunShuffling) PostData(3, fListTriggeredBFS);
}

//________________________________________________________________________
void AliAnalysisTaskTriggeredBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  Float_t fCentrality           = 0.;
  
  // vector holding the charges/kinematics of all tracks (charge,y,eta,phi,p0,p1,p2,pt,E)
  vector<Double_t> *chargeVector[9];          // original charge
  vector<Double_t> *chargeVectorShuffled[9];  // shuffled charge

  for(Int_t i = 0; i < 9; i++){
    chargeVector[i]         = new vector<Double_t>;
    chargeVectorShuffled[i] = new vector<Double_t>;
  }
  
  Double_t v_charge;
  Double_t v_y;
  Double_t v_eta;
  Double_t v_phi;
  Double_t v_p[3];
  Double_t v_pt;
  Double_t v_E;

  // -------------------------------------------------------------		     
  // AOD analysis (vertex and track cuts also here!!!!)
  if(gAnalysisLevel == "AOD") {
    AliAODEvent* aodEventMain = dynamic_cast<AliAODEvent*>(InputEvent()); 
    if(!aodEventMain) {
      AliError("aodEventMain not available");
      return;
    }
    
    
    AliAODHeader *aodHeaderMain = aodEventMain->GetHeader();
    
    // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
    fHistEventStats->Fill(1); //all events
    
    Bool_t isSelectedMain = kTRUE;
    
    if(fUseOfflineTrigger)
      isSelectedMain = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    
    if(isSelectedMain) {
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
      
      const AliAODVertex *vertexMain = aodEventMain->GetPrimaryVertex();
      
      if(vertexMain) {
	Double32_t fCovMain[6];
	vertexMain->GetCovarianceMatrix(fCovMain);
	
	if(vertexMain->GetNContributors() > 0) {
	  if(fCovMain[5] != 0) {
	    fHistEventStats->Fill(3); //events with a proper vertex
	    if(TMath::Abs(vertexMain->GetX()) < fVxMax) {
	      if(TMath::Abs(vertexMain->GetY()) < fVyMax) {
		if(TMath::Abs(vertexMain->GetZ()) < fVzMax) {
		  fHistEventStats->Fill(4); //analyzed events
		  fHistVx->Fill(vertexMain->GetX());
		  fHistVy->Fill(vertexMain->GetY());
		  fHistVz->Fill(vertexMain->GetZ());
		  
		  // Loop over tracks in main event
		  for (Int_t iTracksMain = 0; iTracksMain < aodEventMain->GetNumberOfTracks(); iTracksMain++) {
		    AliAODTrack* aodTrackMain = dynamic_cast<AliAODTrack *>(aodEventMain->GetTrack(iTracksMain));
		    if (!aodTrackMain) {
		      AliError(Form("Could not receive track %d", iTracksMain));
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
		    
		    Float_t DCAxy = aodTrackMain->DCA();      // this is the DCA from global track (not exactly what is cut on)
		    Float_t DCAz  = aodTrackMain->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
		    
		    
		    // Kinematics cuts from ESD track cuts
		    if( v_pt < fPtMin || v_pt > fPtMax)      continue;
		    if( v_eta < fEtaMin || v_eta > fEtaMax)  continue;
		    
		    // Extra DCA cuts (for systematic studies [!= -1])
		    if( fDCAxyCut != -1 && fDCAzCut != -1){
		      if(TMath::Sqrt((DCAxy*DCAxy)/(fDCAxyCut*fDCAxyCut)+(DCAz*DCAz)/(fDCAzCut*fDCAzCut)) > 1 ){
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
		    fHistDCA->Fill(DCAz,DCAxy);
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

		    if(fRunShuffling) {
		      chargeVectorShuffled[0]->push_back(v_charge);
		      chargeVectorShuffled[1]->push_back(v_y);
		      chargeVectorShuffled[2]->push_back(v_eta);
		      chargeVectorShuffled[3]->push_back(v_phi);
		      chargeVectorShuffled[4]->push_back(v_p[0]);
		      chargeVectorShuffled[5]->push_back(v_p[1]);
		      chargeVectorShuffled[6]->push_back(v_p[2]);
		      chargeVectorShuffled[7]->push_back(v_pt);
		      chargeVectorShuffled[8]->push_back(v_E);
		    }
		    
		  } //track loop
		  
		  // calculate balance function
		  fBalance->FillBalance(fCentrality,chargeVector);
		  
		  // calculate shuffled balance function
		  if(fRunShuffling) {
		    random_shuffle(chargeVectorShuffled[0]->begin(), chargeVectorShuffled[0]->end());
		    fShuffledBalance->FillBalance(fCentrality,chargeVectorShuffled);
		  }

		  // clean charge vector afterwards
		  for(Int_t i = 0; i < 9; i++){		       
		    chargeVector[i]->clear();
		    chargeVectorShuffled[i]->clear();
		  }

		}//Vz cut
	      }//Vy cut
	    }//Vx cut
	  }//proper vertex resolution
	}//proper number of contributors
      }//vertex object valid
    }//triggered event 
  }//AOD analysis
  else{
    AliError("Triggered Balance Function analysis only for AODs!");
  }
}     

//________________________________________________________________________
void  AliAnalysisTaskTriggeredBF::FinishTaskOutput(){

  if (!fBalance) {
    AliError("fBalance not available");
    return;
  }  
  if(fRunShuffling) {
    if (!fShuffledBalance) {
      AliError("fShuffledBalance not available");
      return;
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskTriggeredBF::Terminate(Option_t *) {
  // Called once at the end of the query

  // not implemented ...

}

void AliAnalysisTaskTriggeredBF::UserExecMix(Option_t *)
{

  // not yet done for event mixing!
  return;

}

