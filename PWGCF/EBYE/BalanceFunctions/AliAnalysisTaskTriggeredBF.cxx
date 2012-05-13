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

#include "AliEventPoolManager.h" 

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
  fRunMixing(kFALSE),
  fMixingTracks(50000),
  fMixedBalance(0),
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
  if(fRunMixing) {
    if(!fMixedBalance) {
      fMixedBalance = new AliBalanceTriggered();
      fMixedBalance->SetAnalysisLevel("AOD");
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
  if(fRunMixing) {
    fListTriggeredBFM = new TList();
    fListTriggeredBFM->SetName("listTriggeredBFMixed");
    fListTriggeredBFM->SetOwner();
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

  if(fRunMixing) {
    if(!fMixedBalance->GetHistNp()) {
      AliWarning("Histograms (mixing) not yet initialized! --> Will be done now");
      AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
      fMixedBalance->InitHistograms();
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

  if(fRunMixing) {
    fListTriggeredBFM->Add(fMixedBalance->GetHistNp());
    fListTriggeredBFM->Add(fMixedBalance->GetHistNn());
    fListTriggeredBFM->Add(fMixedBalance->GetHistNpn());
    fListTriggeredBFM->Add(fMixedBalance->GetHistNnn());
    fListTriggeredBFM->Add(fMixedBalance->GetHistNpp());
    fListTriggeredBFM->Add(fMixedBalance->GetHistNnp());
  }  


  // Event Mixing
  Int_t trackDepth = fMixingTracks; 
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
   
  Double_t centralityBins[] = {0,1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
  Double_t* centbins = centralityBins;
  Int_t nCentralityBins  = 26;

  
  // bins for second buffer are shifted by 100 cm
  Double_t vertexBins[] = {-10, -7, -5, -3, -1, 1, 3, 5, 7, 10}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
  Double_t* vtxbins = vertexBins;
  Int_t nVertexBins  = 9;

  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins);


  // Post output data.
  PostData(1, fList);
  PostData(2, fListTriggeredBF);
  if(fRunShuffling) PostData(3, fListTriggeredBFS);
  if(fRunMixing) PostData(4, fListTriggeredBFM);
}

//________________________________________________________________________
void AliAnalysisTaskTriggeredBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  Float_t fCentrality = 0.;  

  // -------------------------------------------------------------		     
  // AOD analysis (vertex and track cuts also here!!!!)
  if(gAnalysisLevel == "AOD") {
    AliVEvent* eventMain = dynamic_cast<AliVEvent*>(InputEvent()); 
    if(!eventMain) {
      AliError("eventMain not available");
      return;
    }

    // check event cuts and fill event histograms
    if((fCentrality = IsEventAccepted(eventMain)) < 0){
      return;
    }
    
    // get the accepted tracks in main event
    TObjArray *tracksMain = GetAcceptedTracks(eventMain);

    // store charges of all accepted tracks, shuffle and reassign (two extra loops!)
    TObjArray* tracksShuffled = NULL;
    if(fRunShuffling){
      tracksShuffled = GetShuffledTracks(tracksMain);
    }
    
    // Event mixing --> UPDATE POOL IS MISSING!!!
    if (fRunMixing)
      {
        // 1. First get an event pool corresponding in mult (cent) and
        //    zvertex to the current event. Once initialized, the pool
        //    should contain nMix (reduced) events. This routine does not
        //    pre-scan the chain. The first several events of every chain
        //    will be skipped until the needed pools are filled to the
        //    specified depth. If the pool categories are not too rare, this
        //    should not be a problem. If they are rare, you could lose`
	//    statistics.
	
	// 2. Collect the whole pool's content of tracks into one TObjArray
	//    (bgTracks), which is effectively a single background super-event.
	
	// 3. The reduced and bgTracks arrays must both be passed into
	//    FillCorrelations(). Also nMix should be passed in, so a weight
	//    of 1./nMix can be applied.
	
	AliEventPool* pool = fPoolMgr->GetEventPool(fCentrality, eventMain->GetPrimaryVertex()->GetZ());
	
	if (!pool)
	  AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality, eventMain->GetPrimaryVertex()->GetZ()));
	
	//pool->SetDebug(1);
	
	if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5) 
	  {
	    
	    Int_t nMix = pool->GetCurrentNEvents();
	    cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
	    
	    //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
	    //((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
	    //if (pool->IsReady())
	    //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(3);
	    
	    // Fill mixed-event histos here  
	    for (Int_t jMix=0; jMix<nMix; jMix++) 
	      {
		TObjArray* tracksMixed = pool->GetEvent(jMix);
		fMixedBalance->FillBalance(fCentrality,tracksMixed); //NOW ONLY THE MIXED EVENT ITSELF IS FILLED --> DO ONE TRACK OF MAIN AND ONE OF MIXED (LIKE UEHISTOGRAMS!!!!)
	      }
	  }
      }
    
    // calculate balance function
    fBalance->FillBalance(fCentrality,tracksMain);//,chargeVectorMixed); // here comes the mixing... in some time
    
    // calculate shuffled balance function
    if(fRunShuffling && tracksShuffled != NULL) {
       fShuffledBalance->FillBalance(fCentrality,tracksShuffled);
    }
    
  }//AOD analysis
  else{
    AliError("Triggered Balance Function analysis only for AODs!");
  }
}     

//________________________________________________________________________
Float_t AliAnalysisTaskTriggeredBF::IsEventAccepted(AliVEvent *event){
  // Checks the Event cuts
  // Fills Event statistics histograms
  
  // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
  fHistEventStats->Fill(1); //all events

  Bool_t isSelectedMain = kTRUE;
  Float_t fCentrality = -1.;
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  
  if(fUseOfflineTrigger)
    isSelectedMain = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  if(isSelectedMain) {
    fHistEventStats->Fill(2); //triggered events
    
    //Centrality stuff 
    if(fUseCentrality) {
      if(gAnalysisLevel == "AOD") { //centrality in AOD header
	AliAODHeader *header = (AliAODHeader*) event->GetHeader();
	fCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());

	// QA for centrality estimators
	fHistCentStats->Fill(0.,header->GetCentralityP()->GetCentralityPercentile("V0M"));
	fHistCentStats->Fill(1.,header->GetCentralityP()->GetCentralityPercentile("FMD"));
	fHistCentStats->Fill(2.,header->GetCentralityP()->GetCentralityPercentile("TRK"));
	fHistCentStats->Fill(3.,header->GetCentralityP()->GetCentralityPercentile("TKL"));
	fHistCentStats->Fill(4.,header->GetCentralityP()->GetCentralityPercentile("CL0"));
	fHistCentStats->Fill(5.,header->GetCentralityP()->GetCentralityPercentile("CL1"));
	fHistCentStats->Fill(6.,header->GetCentralityP()->GetCentralityPercentile("V0MvsFMD"));
	fHistCentStats->Fill(7.,header->GetCentralityP()->GetCentralityPercentile("TKLvsV0M"));
	fHistCentStats->Fill(8.,header->GetCentralityP()->GetCentralityPercentile("ZEMvsZDC"));
	
	// centrality QA (V0M)
	fHistV0M->Fill(event->GetVZEROData()->GetMTotV0A(), event->GetVZEROData()->GetMTotV0C());
	
	// centrality QA (reference tracks)
	fHistRefTracks->Fill(0.,header->GetRefMultiplicity());
	fHistRefTracks->Fill(1.,header->GetRefMultiplicityPos());
	fHistRefTracks->Fill(2.,header->GetRefMultiplicityNeg());
	fHistRefTracks->Fill(3.,header->GetTPConlyRefMultiplicity());
	fHistRefTracks->Fill(4.,header->GetNumberOfITSClusters(0));
	fHistRefTracks->Fill(5.,header->GetNumberOfITSClusters(1));
	fHistRefTracks->Fill(6.,header->GetNumberOfITSClusters(2));
	fHistRefTracks->Fill(7.,header->GetNumberOfITSClusters(3));
	fHistRefTracks->Fill(8.,header->GetNumberOfITSClusters(4));
      }
    }
    
    
    const AliVVertex *vertex = event->GetPrimaryVertex();
    
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

		// take only events inside centrality class
		if((fCentrality > fCentralityPercentileMin) && (fCentrality < fCentralityPercentileMax)){
		  return fCentrality;		
		}//centrality class
	      }//Vz cut
	    }//Vy cut
	  }//Vx cut
	}//proper vertex resolution
      }//proper number of contributors
    }//vertex object valid
  }//triggered event 
  
  // in all other cases return -1 (event not accepted)
  return -1;
}

//________________________________________________________________________
TObjArray* AliAnalysisTaskTriggeredBF::GetAcceptedTracks(AliVEvent *event){
  // Returns TObjArray with tracks after all track cuts (only for AOD!)
  // Fills QA histograms

  //output TObjArray holding all good tracks
  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted->SetOwner(kTRUE);

  Double_t v_charge;
  Double_t v_eta;
  Double_t v_phi;
  Double_t v_pt;
  
  // Loop over tracks in event
  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
    if (!aodTrack) {
      AliError(Form("Could not receive track %d", iTracks));
      continue;
    }
    
    // AOD track cuts
    
    // For ESD Filter Information: ANALYSIS/macros/AddTaskESDfilter.C
    // take only TPC only tracks 
    fHistTrackStats->Fill(aodTrack->GetFilterMap());
    if(!aodTrack->TestFilterBit(nAODtrackCutBit)) continue;
    
    v_charge = aodTrack->Charge();
    v_eta    = aodTrack->Eta();
    v_phi    = aodTrack->Phi() * TMath::RadToDeg();
    v_pt     = aodTrack->Pt();
    
    Float_t DCAxy = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
    Float_t DCAz  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
    
    
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
    
    // add the track to the TObjArray
    tracksAccepted->Add(new AliBFBasicParticle(v_eta, v_phi, v_pt, v_charge));
  }

  return tracksAccepted;
}

//________________________________________________________________________
TObjArray* AliAnalysisTaskTriggeredBF::GetShuffledTracks(TObjArray *tracks){
  // Clones TObjArray and returns it with tracks after shuffling the charges

  TObjArray* tracksShuffled = new TObjArray;
  tracksShuffled->SetOwner(kTRUE);

  vector<Short_t> *chargeVector = new vector<Short_t>;   //original charge of accepted tracks 

  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* track = (AliVParticle*) tracks->At(i);
    chargeVector->push_back(track->Charge());
  }  
 
  random_shuffle(chargeVector->begin(), chargeVector->end());
  
  for(Int_t i = 0; i < tracks->GetEntriesFast(); i++){
    AliVParticle* track = (AliVParticle*) tracks->At(i);
    tracksShuffled->Add(new AliBFBasicParticle(track->Eta(), track->Phi(), track->Pt(),chargeVector->at(i)));
  }
   
  return tracksShuffled;
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

