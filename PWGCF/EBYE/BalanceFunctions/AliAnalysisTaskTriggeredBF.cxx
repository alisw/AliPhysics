#include <vector>
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

#include <random>

// Analysis task for the TriggeredBF code
// Authors: Panos.Christakoglou@nikhef.nl, m.weber@cern.ch

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// For the V0 part:
// --> AliAnalysisTaskExtractV0AOD (by david.chinellato@gmail.com)
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

using std::cout;
using std::endl;
using std::vector;

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
  fPoolMgr(0),
  fRunV0(kFALSE),
  fPIDResponse(0),
  fPIDCombined(0),
  fList(0),
  fListTriggeredBF(0),
  fListTriggeredBFS(0),
  fListTriggeredBFM(0),
  fHistListPIDQA(0),
  fHistListV0(0),
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
  fHistV0MultiplicityBeforeTrigSel(0),
  fHistV0MultiplicityForTrigEvt(0),
  fHistV0MultiplicityForSelEvt(0),
  fHistV0MultiplicityForSelEvtNoTPCOnly(0),
  fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup(0),
  fHistMultiplicityBeforeTrigSel(0),
  fHistMultiplicityForTrigEvt(0),
  fHistMultiplicity(0),
  fHistMultiplicityNoTPCOnly(0),
  fHistMultiplicityNoTPCOnlyNoPileup(0),
  fHistV0InvMassK0(0),
  fHistV0InvMassLambda(0),
  fHistV0InvMassAntiLambda(0),
  fHistV0Armenteros(0),
  fHistV0SelInvMassK0(0),
  fHistV0SelInvMassLambda(0),
  fHistV0SelInvMassAntiLambda(0),
  fHistV0SelArmenteros(0),
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
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskTriggeredBF::~AliAnalysisTaskTriggeredBF() {

  // Destructor

}

//________________________________________________________________________
void AliAnalysisTaskTriggeredBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

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

  //------------------------------------------------
  // V0 Multiplicity Histograms
  //------------------------------------------------
  if(fRunV0){
    fHistListV0 = new TList();
    fHistListV0->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
    
    if(! fHistV0MultiplicityBeforeTrigSel) {
      fHistV0MultiplicityBeforeTrigSel = new TH1F("fHistV0MultiplicityBeforeTrigSel", 
						  "V0s per event (before Trig. Sel.);Nbr of V0s/Evt;Events", 
						  25, 0, 25);
      fHistListV0->Add(fHistV0MultiplicityBeforeTrigSel);
    }
    
    if(! fHistV0MultiplicityForTrigEvt) {
      fHistV0MultiplicityForTrigEvt = new TH1F("fHistV0MultiplicityForTrigEvt", 
					       "V0s per event (for triggered evt);Nbr of V0s/Evt;Events", 
					       25, 0, 25);
      fHistListV0->Add(fHistV0MultiplicityForTrigEvt);
    }
    
    if(! fHistV0MultiplicityForSelEvt) {
      fHistV0MultiplicityForSelEvt = new TH1F("fHistV0MultiplicityForSelEvt", 
					      "V0s per event;Nbr of V0s/Evt;Events", 
					      25, 0, 25);
      fHistListV0->Add(fHistV0MultiplicityForSelEvt);
    }
    
    if(! fHistV0MultiplicityForSelEvtNoTPCOnly) {
    fHistV0MultiplicityForSelEvtNoTPCOnly = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnly", 
						     "V0s per event;Nbr of V0s/Evt;Events", 
						     25, 0, 25);
    fHistListV0->Add(fHistV0MultiplicityForSelEvtNoTPCOnly);
    }
    
    if(! fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup) {
      fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup", 
							     "V0s per event;Nbr of V0s/Evt;Events", 
							       25, 0, 25);
      fHistListV0->Add(fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup);
    }
    
    //------------------------------------------------
    // Track Multiplicity Histograms
    //------------------------------------------------
    
    if(! fHistMultiplicityBeforeTrigSel) {
      fHistMultiplicityBeforeTrigSel = new TH1F("fHistMultiplicityBeforeTrigSel", 
					      "Tracks per event;Nbr of Tracks;Events", 
						200, 0, 200); 		
      fHistListV0->Add(fHistMultiplicityBeforeTrigSel);
    }
    if(! fHistMultiplicityForTrigEvt) {
      fHistMultiplicityForTrigEvt = new TH1F("fHistMultiplicityForTrigEvt", 
					     "Tracks per event;Nbr of Tracks;Events", 
					     200, 0, 200); 		
      fHistListV0->Add(fHistMultiplicityForTrigEvt);
    }
    if(! fHistMultiplicity) {
      fHistMultiplicity = new TH1F("fHistMultiplicity", 
				   "Tracks per event;Nbr of Tracks;Events", 
				   200, 0, 200); 		
      fHistListV0->Add(fHistMultiplicity);
    }
    if(! fHistMultiplicityNoTPCOnly) {
      fHistMultiplicityNoTPCOnly = new TH1F("fHistMultiplicityNoTPCOnly", 
					    "Tracks per event;Nbr of Tracks;Events", 
					    200, 0, 200); 		
    fHistListV0->Add(fHistMultiplicityNoTPCOnly);
    }
    if(! fHistMultiplicityNoTPCOnlyNoPileup) {
      fHistMultiplicityNoTPCOnlyNoPileup = new TH1F("fHistMultiplicityNoTPCOnlyNoPileup", 
						    "Tracks per event;Nbr of Tracks;Events", 
						    200, 0, 200); 		
      fHistListV0->Add(fHistMultiplicityNoTPCOnlyNoPileup);
    }

    //------------------------------------------------
    // V0 selection Histograms (before)
    //------------------------------------------------
    if(!fHistV0InvMassK0) {
      fHistV0InvMassK0 = new TH1F("fHistV0InvMassK0",
				  "Invariant Mass for K0;Mass (GeV/c^{2});Events",
				  200,0,2);
      fHistListV0->Add(fHistV0InvMassK0);
    }
    if(!fHistV0InvMassLambda) {
      fHistV0InvMassLambda = new TH1F("fHistV0InvMassLambda",
				  "Invariant Mass for Lambda;Mass (GeV/c^{2});Events",
				  200,0,2);
      fHistListV0->Add(fHistV0InvMassLambda);
    }
    if(!fHistV0InvMassAntiLambda) {
      fHistV0InvMassAntiLambda = new TH1F("fHistV0InvMassAntiLambda",
				  "Invariant Mass for AntiLambda;Mass (GeV/c^{2});Events",
				  200,0,2);
      fHistListV0->Add(fHistV0InvMassAntiLambda);
    }
    if(!fHistV0Armenteros) {
      fHistV0Armenteros = new TH2F("fHistV0Armenteros",
				  "Armenteros plot;#alpha;q_{t}",
				   200,-1,1,200,0,0.5);
      fHistListV0->Add(fHistV0Armenteros);
    }
    
    //------------------------------------------------
    // V0 selection Histograms (after)
    //------------------------------------------------
    if(!fHistV0SelInvMassK0) {
      fHistV0SelInvMassK0 = new TH1F("fHistV0SelInvMassK0",
				  "Invariant Mass for K0;Mass (GeV/c^{2});Events",
				  200,0,2);
      fHistListV0->Add(fHistV0SelInvMassK0);
    }
    if(!fHistV0SelInvMassLambda) {
      fHistV0SelInvMassLambda = new TH1F("fHistV0SelInvMassLambda",
				  "Invariant Mass for Lambda;Mass (GeV/c^{2});Events",
				  200,0,2);
      fHistListV0->Add(fHistV0SelInvMassLambda);
    }
    if(!fHistV0SelInvMassAntiLambda) {
      fHistV0SelInvMassAntiLambda = new TH1F("fHistV0SelInvMassAntiLambda",
				  "Invariant Mass for AntiLambda;Mass (GeV/c^{2});Events",
				  200,0,2);
      fHistListV0->Add(fHistV0SelInvMassAntiLambda);
    }
    if(!fHistV0SelArmenteros) {
      fHistV0SelArmenteros = new TH2F("fHistV0SelArmenteros",
				  "Armenteros plot;#alpha;q_{t}",
				   200,-1,1,200,0,0.5);
      fHistListV0->Add(fHistV0SelArmenteros);
    }
  }//V0
    
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

  // PID Response task active?
  if(fRunV0) {
    fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
  }

  // Event Mixing
  Int_t trackDepth = fMixingTracks; 
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
   
  Double_t centralityBins[] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,90.,100.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
  Double_t* centbins        = centralityBins;
  Int_t nCentralityBins     = sizeof(centralityBins) / sizeof(Double_t) - 1;
  
  // bins for second buffer are shifted by 100 cm
  Double_t vertexBins[] = {-10., -7., -5., -3., -1., 1., 3., 5., 7., 10.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
  Double_t* vtxbins     = vertexBins;
  Int_t nVertexBins     = sizeof(vertexBins) / sizeof(Double_t) - 1;

  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins);


  // Post output data.
  PostData(1, fList);
  PostData(2, fListTriggeredBF);
  if(fRunShuffling) PostData(3, fListTriggeredBFS);
  if(fRunMixing) PostData(4, fListTriggeredBFM);
  if(fRunV0) PostData(5,fHistListV0);

  TH1::AddDirectory(oldStatus);

}

//________________________________________________________________________
void AliAnalysisTaskTriggeredBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  Float_t fCentrality = -1.;  

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
    TObjArray *tracksMain = NULL;
    if(fRunV0) tracksMain = GetAcceptedV0s(eventMain);
    else       tracksMain = GetAcceptedTracks(eventMain);

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
	
	if (!pool){
	  AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality, eventMain->GetPrimaryVertex()->GetZ()));
	}
	else{

	//pool->SetDebug(1);
	
	  if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5){ 
	    
	    
	    Int_t nMix = pool->GetCurrentNEvents();
	    //cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
	    
	    //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
	    //((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
	    //if (pool->IsReady())
	    //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(3);
	    
	    // Fill mixed-event histos here  
	    for (Int_t jMix=0; jMix<nMix; jMix++) 
	      {
		TObjArray* tracksMixed = pool->GetEvent(jMix);
		fMixedBalance->FillBalance(fCentrality,tracksMain,tracksMixed); 
	      }
	  }
	  
	  // Update the Event pool
	  pool->UpdatePool(tracksMain);
	  //pool->PrintInfo();
	  
	}//pool NULL check  
      }//run mixing
    
    // calculate balance function
    fBalance->FillBalance(fCentrality,tracksMain,NULL);
    
    // calculate shuffled balance function
    if(fRunShuffling && tracksShuffled != NULL) {
      fShuffledBalance->FillBalance(fCentrality,tracksShuffled,NULL);
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
  Int_t nV0s          = event->GetNumberOfV0s();
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  
  if(fUseOfflineTrigger)
    isSelectedMain = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  //V0 QA histograms (before trigger selection)
  if(fRunV0){
    fHistMultiplicityBeforeTrigSel->Fill ( -1 );
    fHistV0MultiplicityBeforeTrigSel->Fill ( nV0s );
  }
  
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

	//V0 QA histograms (after trigger selection)
	if(fRunV0){
	  fHistMultiplicityForTrigEvt->Fill ( fCentrality );
	  fHistV0MultiplicityForTrigEvt->Fill ( nV0s );
	}
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



		//V0 QA histograms (vertex Z check)
		if(fRunV0){
		  fHistV0MultiplicityForSelEvt ->Fill( nV0s );
		  fHistMultiplicity->Fill(fCentrality);

		  //V0 QA histograms (Only look at events with well-established PV)
		  const AliAODVertex *lPrimarySPDVtx = ((AliAODEvent*)event)->GetPrimaryVertexSPD();
		  if(lPrimarySPDVtx){
		    fHistMultiplicityNoTPCOnly->Fill ( fCentrality );
		    fHistV0MultiplicityForSelEvtNoTPCOnly->Fill ( nV0s );
		    
		    //V0 QA histograms (Pileup Rejection)
		    // FIXME : quality selection regarding pile-up rejection 
		    fHistMultiplicityNoTPCOnlyNoPileup->Fill(fCentrality);
		    fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup ->Fill( nV0s );

		  }
		  else{
		    return -1;
		  }
		}
		
		
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

  Short_t vCharge = 0;
  Double_t vEta    = 0.;
  Double_t vPhi    = 0.;
  Double_t vPt     = 0.;
  
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
    
    vCharge = aodTrack->Charge();
    vEta    = aodTrack->Eta();
    vPhi    = aodTrack->Phi() * TMath::RadToDeg();
    vPt     = aodTrack->Pt();
    
    Float_t dcaXY = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
    Float_t dcaZ  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
    
    
    // Kinematics cuts from ESD track cuts
    if( vPt < fPtMin || vPt > fPtMax)      continue;
    if( vEta < fEtaMin || vEta > fEtaMax)  continue;
    
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
    
    // fill QA histograms
    fHistClus->Fill(aodTrack->GetITSNcls(),aodTrack->GetTPCNcls());
    fHistDCA->Fill(dcaZ,dcaXY);
    fHistChi2->Fill(aodTrack->Chi2perNDF());
    fHistPt->Fill(vPt);
    fHistEta->Fill(vEta);
    fHistPhi->Fill(vPhi);
    
    // add the track to the TObjArray
    tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge,1.));
  }

  return tracksAccepted;
}

//________________________________________________________________________
TObjArray* AliAnalysisTaskTriggeredBF::GetAcceptedV0s(AliVEvent *event){
  // Returns TObjArray with tracks after all track cuts (only for AOD!)
  // Fills QA histograms

  //output TObjArray holding all good tracks
  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted->SetOwner(kTRUE);

  Short_t vCharge = 0;
  Double_t vEta    = 0.;
  Double_t vPhi    = 0.;
  Double_t vPt     = 0.;
  
  //------------------------------------------------
  // MAIN LAMBDA LOOP STARTS HERE (basically a copy of AliAnalysisTaskExtractV0AOD)
  //------------------------------------------------

  // parameters (for the time being hard coded here) --> from David for EbyE Lambdas
  Bool_t fkUseOnTheFly = kFALSE;
  Double_t fRapidityBoundary  = 0.5; 
  Double_t fCutDaughterEta    = 0.8;
  Double_t fCutV0Radius       = 0.9;
  Double_t fCutDCANegToPV     = 0.1;
  Double_t fCutDCAPosToPV     = 0.1;
  Double_t fCutDCAV0Daughters = 1.0;
  Double_t fCutV0CosPA        = 0.9995;
  Double_t fMassLambda        = 1.115683;
  Double_t fCutMassLambda     = 0.007;
  Double_t fCutProperLifetime = 3*7.9;
  Double_t fCutLeastNumberOfCrossedRows = 70;
  Double_t fCutLeastNumberOfCrossedRowsOverFindable = 0.8;
  Double_t fCutTPCPIDNSigmasProton  = 3.0;
  Double_t fCutTPCPIDNSigmasPion    = 5.0;


  //Variable definition
  Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
  Double_t lChi2V0 = 0;
  Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
  Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
  Double_t lV0CosineOfPointingAngle = 0;
  Double_t lV0Radius = 0, lPt = 0;
  Double_t lEta = 0, lPhi = 0;
  Double_t lRap = 0, lRapK0Short = 0, lRapLambda = 0;
  Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
  Double_t lAlphaV0 = 0, lPtArmV0 = 0;
  
  Double_t fMinV0Pt = 0; 
  Double_t fMaxV0Pt = 100; 
  

  
  // some event observables
  Int_t nv0s = event->GetNumberOfV0s();
  Double_t tPrimaryVtxPosition[3];
  const AliVVertex *primaryVtx = event->GetPrimaryVertex();
  tPrimaryVtxPosition[0] = primaryVtx->GetX();
  tPrimaryVtxPosition[1] = primaryVtx->GetY();
  tPrimaryVtxPosition[2] = primaryVtx->GetZ();


  //loop over V0s  
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) 
    {// This is the begining of the V0 loop
      AliAODv0 *v0 = ((AliAODEvent*)event)->GetV0(iV0);
      if (!v0) continue;

      //Obsolete at AOD level... 
      //---> Fix On-the-Fly candidates, count how many swapped
      //if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ){
      //  fHistSwappedV0Counter -> Fill( 1 );
      //}else{
      //  fHistSwappedV0Counter -> Fill( 0 ); 
      //}
      //if ( fkUseOnTheFly ) CheckChargeV0(v0); 
      
      Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0); 
      Double_t tV0mom[3];
      v0->GetPxPyPz( tV0mom ); 
      Double_t lV0TotalMomentum = TMath::Sqrt(
					      tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );
      
      lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
      lPt = v0->Pt();
      lEta = v0->Eta();
      lPhi = v0->Phi()*TMath::RadToDeg();
      lRapK0Short = v0->RapK0Short();
      lRapLambda  = v0->RapLambda();
      lRap        = lRapLambda;//v0->Y(); //FIXME!!!
      if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;
      
      //UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPosID());
      //UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetPosID());

      Double_t lMomPos[3]; //v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      Double_t lMomNeg[3]; //v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
      lMomPos[0] = v0->MomPosX();
      lMomPos[1] = v0->MomPosY();
      lMomPos[2] = v0->MomPosZ();
      lMomNeg[0] = v0->MomNegX();
      lMomNeg[1] = v0->MomNegY();
      lMomNeg[2] = v0->MomNegZ();
      
      AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
      AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
      if (!pTrack || !nTrack) {
	AliError("ERROR: Could not retreive one of the daughter track");
	continue;
      }

      //Daughter Eta for Eta selection, afterwards
      Double_t lNegEta = nTrack->Eta();
      Double_t lPosEta = pTrack->Eta();
      
      // Filter like-sign V0 (next: add counter and distribution)
      if ( pTrack->Charge() == nTrack->Charge()){
	continue;
      } 
      
      //Quick test this far! 
      

      //________________________________________________________________________
      // Track quality cuts 
      Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
      Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
      Float_t lLeastNbrCrossedRows =  (lPosTrackCrossedRows>lNegTrackCrossedRows) ? lNegTrackCrossedRows : lPosTrackCrossedRows;

      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
      if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
      
      if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;
      
      //Findable clusters > 0 condition
      if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;
      
      //Compute ratio Crossed Rows / Findable clusters
      //Note: above test avoids division by zero! 
      Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF())); 
      Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF())); 
      Float_t lLeastNbrCrossedRowsOverFindable = (lPosTrackCrossedRowsOverFindable>lNegTrackCrossedRowsOverFindable) ? lNegTrackCrossedRowsOverFindable : lPosTrackCrossedRowsOverFindable;

      //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
      if ( lLeastNbrCrossedRowsOverFindable < 0.8) continue;
      
      //End track Quality Cuts
      //________________________________________________________________________
      
      
      lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
      lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
          
      lOnFlyStatus = v0->GetOnFlyStatus();
      lChi2V0 = v0->Chi2V0();
      lDcaV0Daughters = v0->DcaV0Daughters();
      lDcaV0ToPrimVertex = v0->DcaV0ToPrimVertex();
      lV0CosineOfPointingAngle = v0->CosPointingAngle(tPrimaryVtxPosition);
      
      // Distance over total momentum
      Double_t lDistOverTotMom = TMath::Sqrt(
				    TMath::Power( tDecayVertexV0[0] - tPrimaryVtxPosition[0] , 2) +
				    TMath::Power( tDecayVertexV0[1] - tPrimaryVtxPosition[1] , 2) +
				    TMath::Power( tDecayVertexV0[2] - tPrimaryVtxPosition[2] , 2)
				    );
      lDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure
      
      
      // Getting invariant mass infos directly from ESD
      lInvMassK0s        = v0->MassK0Short();
      lInvMassLambda     = v0->MassLambda();
      lInvMassAntiLambda = v0->MassAntiLambda();
      lAlphaV0 = v0->AlphaV0();
      lPtArmV0 = v0->PtArmV0();

      //Official means of acquiring N-sigmas 
      Double_t lNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
      Double_t lNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
      Double_t lNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
      Double_t lNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );

      //V0 QA histograms (before V0 selection)
      fHistV0InvMassK0->Fill(lInvMassK0s);
      fHistV0InvMassLambda->Fill(lInvMassLambda);
      fHistV0InvMassAntiLambda->Fill(lInvMassAntiLambda);
      fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
      
      
      //First Selection: Reject OnFly
      if( (lOnFlyStatus == 0 && fkUseOnTheFly == kFALSE) || (lOnFlyStatus != 0 && fkUseOnTheFly == kTRUE ) ){
	

      	//Second Selection: rough 20-sigma band, parametric. 	
      	//K0Short: Enough to parametrize peak broadening with linear function.    
      	Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*lPt; 
      	Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*lPt;
	
      	//Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
      	//[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
      	Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*lPt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*lPt); 
      	Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*lPt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*lPt);
	
      	//Do Selection      
      	if( (lInvMassLambda     < lUpperLimitLambda  && lInvMassLambda     > lLowerLimitLambda     ) || 
      	    (lInvMassAntiLambda < lUpperLimitLambda  && lInvMassAntiLambda > lLowerLimitLambda     ) || 
      	    (lInvMassK0s        < lUpperLimitK0Short && lInvMassK0s        > lLowerLimitK0Short    ) ){


      // 	  //Pre-selection in case this is AA...
      // 	  //if( fkIsNuclear == kFALSE ) fTree->Fill();
      // 	  //if( fkIsNuclear == kTRUE){ 
      // 	  //If this is a nuclear collision___________________
      // 	  // ... pre-filter with TPC, daughter eta selection

	  
	  if( (lInvMassLambda     < lUpperLimitLambda  && lInvMassLambda     > lLowerLimitLambda 
      	       && TMath::Abs(lNSigmasPosProton) < 6.0 && TMath::Abs(lNSigmasNegPion) < 6.0 ) || 
      	      (lInvMassAntiLambda < lUpperLimitLambda  && lInvMassAntiLambda > lLowerLimitLambda 
      	       && TMath::Abs(lNSigmasNegProton) < 6.0 && TMath::Abs(lNSigmasPosPion) < 6.0 ) ||  
      	      (lInvMassK0s        < lUpperLimitK0Short && lInvMassK0s        > lLowerLimitK0Short 
      	       && TMath::Abs(lNSigmasNegPion)   < 6.0 && TMath::Abs(lNSigmasPosPion) < 6.0 ) ){
	    
      	    //insane test
      	    if ( TMath::Abs(lNegEta)<0.8 && TMath::Abs(lPosEta)<0.8 ){

	      // start the fine selection (usually done in post processing, but we don't have time to waste) --> Lambdas!
	      if(
		 TMath::Abs(lRap)<fRapidityBoundary &&
		 TMath::Abs(lNegEta)       <= fCutDaughterEta               &&                   
		 TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
		 lV0Radius                 >= fCutV0Radius                  &&
		 lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
		 lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
		 lDcaV0Daughters           <= fCutDCAV0Daughters            &&
		 lV0CosineOfPointingAngle  >= fCutV0CosPA                   && 
		 fMassLambda*lDistOverTotMom    <= fCutProperLifetime       &&
		 lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
		 lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable &&
		 lPtArmV0 * 5 < TMath::Abs(lAlphaV0)                        && 
		 ((TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmasPion     &&
		  TMath::Abs(lNSigmasPosProton) <= fCutTPCPIDNSigmasProton) ||
		  (TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmasPion     &&
		   TMath::Abs(lNSigmasNegProton) <= fCutTPCPIDNSigmasProton)) 		 
		 )
		{

		  //V0 QA histograms (after V0 selection)
		  fHistV0SelInvMassK0->Fill(lInvMassK0s);
		  fHistV0SelInvMassLambda->Fill(lInvMassLambda);
		  fHistV0SelInvMassAntiLambda->Fill(lInvMassAntiLambda);

		  // this means a V0 candidate is found
		  if(TMath::Abs(lInvMassLambda-fMassLambda) < fCutMassLambda ||
		     TMath::Abs(lInvMassAntiLambda-fMassLambda) < fCutMassLambda){

		    fHistV0SelArmenteros->Fill(lAlphaV0,lPtArmV0);		  

		    vEta    = lEta;
		    vPhi    = lPhi;
		    vPt     = lPt;
		    if(lAlphaV0 > 0) vCharge = 1;
		    if(lAlphaV0 < 0) vCharge = -1;

		    // fill QA histograms
		    fHistPt->Fill(vPt);
		    fHistEta->Fill(vEta);
		    fHistPhi->Fill(vPhi);
		    
		    // add the track to the TObjArray
		    tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge,1.));
		  }
		}
	      }
	  }
	  //}//end nuclear_____________________________________
	}
      }
    }//V0 loop
  
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
 
  std::random_device rd;
  std::default_random_engine engine{rd()};
  std::shuffle(chargeVector->begin(), chargeVector->end(), engine);
  
  for(Int_t i = 0; i < tracks->GetEntriesFast(); i++){
    AliVParticle* track = (AliVParticle*) tracks->At(i);
    tracksShuffled->Add(new AliBFBasicParticle(track->Eta(), track->Phi(), track->Pt(),chargeVector->at(i),1.));
  }

  delete chargeVector;
   
  return tracksShuffled;
}

//________________________________________________________________________
void  AliAnalysisTaskTriggeredBF::FinishTaskOutput(){
  //checks if Balance Function objects are there (needed to write the histograms)
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

