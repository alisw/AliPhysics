#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskBF.h"
#include "AliBalance.h"


// Analysis task for the BF code
// Authors: Panos Cristakoglou@cern.ch

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
  fNClustersTPCCut(-1){
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

  TObjArray *array         = new TObjArray();
  AliESDtrack *track_TPC   = NULL;


  // vector holding the charges of all tracks
  vector<Int_t> chargeVectorShuffle;   // this will be shuffled
  vector<Int_t> chargeVector;          // original charge 
  
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
	//Int_t nCentrality = 0;
	//nCentrality = centrality->GetCentralityClass5(fCentralityEstimator.Data());
	//cout<<nCentrality<<" "<<centrality->IsEventInCentralityClass(fCentralityPercentileMin,fCentralityPercentileMax,fCentralityEstimator.Data())<<endl;
	
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
		    
		    fHistClus->Fill(track_TPC->GetITSclusters(0),nClustersTPC);
		    fHistDCA->Fill(b[1],b[0]);
		    fHistChi2->Fill(chi2PerClusterTPC);
		    fHistPt->Fill(track_TPC->Pt());
		    fHistEta->Fill(track_TPC->Eta());
		    fHistPhi->Fill(track_TPC->Phi()*TMath::RadToDeg());
		    
		    // fill BF array
		    array->Add(track_TPC);
		    
		    // fill charge vector
		    chargeVector.push_back(track_TPC->Charge());
		    if(fRunShuffling){
		      chargeVectorShuffle.push_back(track_TPC->Charge());
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
	Float_t fCentrality     = aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
	// cout<<fCentralityEstimator.Data()<<" = "<<fCentrality<<" ,  others are V0M =  "
	// 	  << aodHeader->GetCentralityP()->GetCentralityPercentile("V0M")
	// 	  <<"  FMD = "<<aodHeader->GetCentralityP()->GetCentralityPercentile("FMD")
	// 	  <<"  TRK = "<<aodHeader->GetCentralityP()->GetCentralityPercentile("TRK")
	// 	  <<"  TKL = "<<aodHeader->GetCentralityP()->GetCentralityPercentile("TKL")
	// 	  <<"  CL0 ="<<aodHeader->GetCentralityP()->GetCentralityPercentile("CL0")
	// 	  <<"  CL1 ="<<aodHeader->GetCentralityP()->GetCentralityPercentile("CL1")
	// 	  <<"  V0MvsFMD = "<<aodHeader->GetCentralityP()->GetCentralityPercentile("V0MvsFMD")
	// 	  <<"  TKLvsV0M = "<<aodHeader->GetCentralityP()->GetCentralityPercentile("TKLvsV0M")
	// 	  <<"  ZEMvsZDC = "<<aodHeader->GetCentralityP()->GetCentralityPercentile("ZEMvsZDC")
	// 	  <<endl;
	
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
		    
		    Float_t pt  = aodTrack->Pt();
		    Float_t eta = aodTrack->Eta();
		    
		    Float_t DCAxy = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
		    Float_t DCAz  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
		    
		    
		    // Kinematics cuts from ESD track cuts
		    if( pt < fPtMin || pt > fPtMax)      continue;
		    if( eta < fEtaMin || eta > fEtaMax)  continue;
		    
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
		    fHistPt->Fill(pt);
		    fHistEta->Fill(eta);
		    fHistPhi->Fill(aodTrack->Phi()*TMath::RadToDeg());
		    
		    // fill BF array
		    array->Add(aodTrack);
		    
		    // fill charge vector
		    chargeVector.push_back(aodTrack->Charge());
		    if(fRunShuffling) {
		      chargeVectorShuffle.push_back(aodTrack->Charge());
		    }
		  } //track loop
		}//Vz cut
	      }//Vy cut
	    }//Vx cut
	  }//proper vertex resolution
	}//proper number of contributors
      }//vertex object valid
    }//triggered event 
  }//AOD analysis

  //MC analysis
  else if(gAnalysisLevel == "MC") {
    
    AliMCEvent*  mcEvent = MCEvent(); 
    if (!mcEvent) {
      Printf("ERROR: mcEvent not available");
      return;
    }
    
    Printf("There are %d tracks in this event", mcEvent->GetNumberOfPrimaries());
    for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfPrimaries(); iTracks++) {
      AliMCParticle* track = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iTracks));
      if (!track) {
	Printf("ERROR: Could not receive particle %d", iTracks);
	continue;
      }

      if( track->Pt() < fPtMin || track->Pt() > fPtMax)      continue;
      if( track->Eta() < fEtaMin || track->Eta() > fEtaMax)  continue;

      array->Add(track);

      // fill charge vector
      chargeVector.push_back(track->Charge());
      if(fRunShuffling) {
	chargeVectorShuffle.push_back(track->Charge());
      }

    } //track loop
  }//MC analysis
  
  // calculate balance function
  fBalance->CalculateBalance(array,chargeVector);

  if(fRunShuffling) {
    // shuffle charges
    random_shuffle( chargeVectorShuffle.begin(), chargeVectorShuffle.end() );
    fShuffledBalance->CalculateBalance(array,chargeVectorShuffle);
  }
  
  delete array;
  

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
