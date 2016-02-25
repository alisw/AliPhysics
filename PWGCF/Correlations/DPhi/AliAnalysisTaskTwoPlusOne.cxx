#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TTree.h>
#include <TRandom3.h>
//#include <TParameter.h>

#include "AliAnalysisTaskTwoPlusOne.h"
#include "AliCFParticle.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliTwoPlusOneContainer.h"
#include "AliUEHist.h"

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliVParticle.h"
#include "AliCFContainer.h"

#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"

#include "AliEventPoolManager.h"
#include <iostream>

ClassImp( AliAnalysisTaskTwoPlusOne )

//________________________________________________________________________
AliAnalysisTaskTwoPlusOne::AliAnalysisTaskTwoPlusOne(const char *name)
: AliAnalysisTaskSE(name),
  fMixingTracks(10000),
  fMode(0),
  fIsNano(0),
  fAnalyseUE(0x0),
// pointers to UE classes
  fHistos(0x0),
// handlers and events
  fAOD(0x0),
  fMcEvent(0x0),
  fMcHandler(0x0),
  fPoolMgr(0x0),
  fEventCombination(0x0),
  fUsedEvents(0),
// histogram settings
  fListOfHistos(0x0), 
// event QA
  fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default) 
  fZVertex(7.),
  fCentralityMethod("V0M"),
// track cuts
  fTrackEtaCut(0.9),
  fTrackEtaCutMin(-1.),
  fPtMin(0.5),
  fDCAXYCut(0),
  fSharedClusterCut(-1),
  fCrossedRowsCut(-1),
  fFoundFractionCut(-1),
  fFilterBit(0xFF),
  fTrackStatus(0),
  fThreeParticleMixed(0),
  fUseEventCombination(0),
  fUsePP(0),
  fUEHist_name("TwoPlusOne"),
  fCustomBinning(),
  fAlpha(0.2),
  fUseLeadingPt(1),
  fUseAllT1(1),
  fUseBackgroundSameOneSide(0),
  fUseBackgroundSameFromMixedComb(0),
  fUseSmallerPtAssoc(0),
  fRunCorrelations(1),
  fRunIfPoolReady(0),
  fRandomPosition(0),
  fSelectCentrality(1),
  fEfficiencyCorrection(0)
{

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskTwoPlusOne::~AliAnalysisTaskTwoPlusOne()
{
    // Destructor

}

//________________________________________________________________________
void AliAnalysisTaskTwoPlusOne::UserCreateOutputObjects()
{
  // Initialize class with main algorithms, event and track selection. 
  fAnalyseUE = new AliAnalyseLeadingTrackUE();
  fAnalyseUE->SetParticleSelectionCriteria(fFilterBit, kFALSE, fTrackEtaCut, fTrackEtaCutMin, fPtMin);
  fAnalyseUE->SetDCAXYCut(fDCAXYCut);
  fAnalyseUE->SetSharedClusterCut(fSharedClusterCut);
  fAnalyseUE->SetCrossedRowsCut(fCrossedRowsCut);
  fAnalyseUE->SetFoundFractionCut(fFoundFractionCut);
  fAnalyseUE->SetTrackStatus(fTrackStatus);
  fAnalyseUE->SetDebug(fDebug); 
  fAnalyseUE->DefineESDCuts(fFilterBit);

  fListOfHistos = new TList();
  fListOfHistos->SetOwner(kTRUE); 

  fHistos = new AliTwoPlusOneContainer("AliTwoPlusOneContainer", fUEHist_name, fCustomBinning, fAlpha);
  fHistos->GetData()->SetTrackEtaCut(fTrackEtaCut);
  fHistos->SetUseLeadingPt(fUseLeadingPt);
  fHistos->SetUseAllT1(fUseAllT1);
  fHistos->SetUseBackgroundSameOneSide(fUseBackgroundSameOneSide);
  fHistos->SetUseSmallerPtAssoc(fUseSmallerPtAssoc);

  if (fEfficiencyCorrection)
    fHistos->SetEfficiencyCorrection(fEfficiencyCorrection);

  fListOfHistos->Add(fHistos);

  fListOfHistos->Add(new TH1F("eventStat", ";;events", 4, -0.5, 3.5));
  fListOfHistos->Add(new TH2F("eventStatCent", ";events;centrality", 4, -0.5, 3.5, 201, 0, 100.5));
  fListOfHistos->Add(new TH2F("mixedDist", ";centrality;tracks;events", 101, 0, 101, 200, 0, fMixingTracks * 1.5));

  PostData(1,fListOfHistos);

  // Add task configuration to output list 
  AddSettingsTree();
  
  // event mixing
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemention of AliEventPoolManager
   
  Int_t nCentralityBins  = fHistos->GetData()->GetTrackHist(AliUEHist::kToward)->GetNBins(3);

  
  Double_t* centralityBins = (Double_t*) fHistos->GetData()->GetTrackHist(AliUEHist::kToward)->GetAxis(3, 0)->GetXbins()->GetArray();
  
  Int_t nZvtxBins = fHistos->GetData()->GetTrackHist(AliUEHist::kToward)->GetNBins(5);
  
  Double_t* zvtxbin = (Double_t*) fHistos->GetData()->GetTrackHist(AliUEHist::kToward)->GetAxis(5, 0)->GetXbins()->GetArray();

  fPoolMgr = new AliEventPoolManager(poolsize, fMixingTracks, nCentralityBins, centralityBins, nZvtxBins, zvtxbin);
  fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);

  fEventCombination = new TObjArray;

  // MC handler
  if(fMode)
    fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}

//________________________________________________________________________
void AliAnalysisTaskTwoPlusOne::UserExec(Option_t *)
{
  // exec (per event)
  fAnalyseUE->NextEvent();

  ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(0);

  // Support for AOD based analysis
  fAOD = dynamic_cast<AliAODEvent*> (InputEvent());

  Double_t centrality = 0;
  AliCentrality *centralityObj = 0;
  
  if (!fUsePP && fAOD){
    centralityObj = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP();
  
    if (centralityObj)
      centrality = centralityObj->GetCentralityPercentile(fCentralityMethod);
    else
      centrality = -1;

    // remove outliers
    if (centrality == 0)
      {
	if (fAOD->GetVZEROData())
	  {
	    Float_t multV0 = 0;
	    for (Int_t i=0; i<64; i++)
	      multV0 += fAOD->GetVZEROData()->GetMultiplicity(i);
	    if (multV0 < 19500)
	      {
		centrality = -1;
		AliInfo("Rejecting event due to too small V0 multiplicity");
	      }
	  }
      }
  }

  if (fMcHandler){
    fMcEvent = fMcHandler->MCEvent();
  }

  if(fMode){
      AliGenEventHeader* eventHeader = GetFirstHeader();
      if (!eventHeader)
      {
	// We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
	// (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
	AliError("Event header not found. Skipping this event.");
	//fHistos->FillEvent(0, AliUEHist::kCFStepAnaTopology);
	return;
      }
      
      AliCollisionGeometry* collGeometry = dynamic_cast<AliCollisionGeometry*> (eventHeader);
      if (!collGeometry)
      {
	eventHeader->Dump();
	AliFatal("Asking for MC_b centrality, but event header has no collision geometry information");
      }
      
      double impact_parameter = collGeometry->ImpactParameter();
      //put centrality on the middle of the bin
      if(impact_parameter<1.60)
	centrality = 0.5;
      else if(impact_parameter<2.27)
	centrality = 1.5;
      else if(impact_parameter<2.79)
	centrality = 2.5;
      else if(impact_parameter<3.22)
	centrality = 3.5;
      else if(impact_parameter<3.60)
	centrality = 4.5;
      else if(impact_parameter<5.09)
	centrality = 7.5;
      else if(impact_parameter<7.20)
	centrality = 15;
      else if(impact_parameter<8.83)
	centrality = 25;
      else if(impact_parameter<10.20)
	centrality = 35;
      else if(impact_parameter<11.40)
	centrality = 45;
      else if(impact_parameter<12.49)
	centrality = 55;
      else if(impact_parameter<13.49)
	centrality = 65;
      else if(impact_parameter<14.44)
	centrality = 75;
      else if(impact_parameter<15.46)
	centrality = 85;
      else 
	centrality = 95;
  }

  if(fIsNano){
    centrality = ((AliNanoAODHeader*) fAOD->GetHeader())->GetCentrality();
  }

  if(centrality>=0)
    ((TH1F*) fListOfHistos->FindObject("eventStatCent"))->Fill((Double_t)0, centrality);

  // Get MC primaries
  TObjArray* tracksMC;

  if(fMode==1){
    Int_t   fParticleSpeciesTrigger = -1;
    TObjArray* tmpList = fAnalyseUE->GetAcceptedParticles(fMcEvent, 0, kTRUE, fParticleSpeciesTrigger, kTRUE);
    tracksMC = CloneAndReduceTrackList(tmpList);
  }


  AliInfo(Form("Centrality is %f", centrality));


  // Vertex selection *************************************************
  if(fMode==0 && !fAnalyseUE->VertexSelection(InputEvent(), fnTracksVertex, fZVertex)) return;

  // optimization
  if (centrality < 0)
    return;

  TObjArray* tracks;
  if(fMode==0)
    tracks = fAnalyseUE->GetAcceptedParticles(InputEvent(), 0, kTRUE, -1, kTRUE);
  else
    tracks  = tracksMC;

  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = CloneAndReduceTrackList(tracks);
  delete tracks;

  Double_t zVtx = 0;

  if(fMode==0){//data mode
    const AliVVertex* vertex = InputEvent()->GetPrimaryVertex();
    zVtx = vertex->GetZ();
  }

  if(fSelectCentrality){
    if((centrality>7.5 && centrality<30)||centrality>50)
      return;
  }

  //at this point of the code the event is acctepted
  //if this run is used to add 30-50% centrality events to the multiplicity of central events this is done here, all other events are skipped. 
  if(fUseEventCombination){
    if(centrality>30&&centrality<=50){

      AddEventCombination(tracksClone);
      //do only continue if there are 4 events in the fEventCombination
      if(fUsedEvents==4)
	tracksClone = fEventCombination;
      else
	return;
    }
  }

  Bool_t applyEfficiency = kFALSE;
  if(fEfficiencyCorrection)
    applyEfficiency = kTRUE;

  fHistos->FillParticleDist(centrality, zVtx, tracksClone, 1.0, applyEfficiency);

  if(fRunCorrelations){

    AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx);

    if (fRunIfPoolReady && !pool)
	AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, zVtx));
    if (!fRunIfPoolReady || pool->IsReady()){  

    fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kSameNS, tracksClone, tracksClone, tracksClone, tracksClone, 1.0, kFALSE, kFALSE, applyEfficiency);//same event for near and away side

    fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::k1plus1, tracksClone, tracksClone, tracksClone, tracksClone, 1.0, kTRUE, kFALSE, applyEfficiency);//get number of possible away side triggers in the trigger area and outside of it
    
    if(!fUseBackgroundSameFromMixedComb)
      fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kBackgroundSameNS, tracksClone, tracksClone, tracksClone, tracksClone, 1.0, kFALSE, kTRUE, applyEfficiency);//background estimation for the same event

    ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(1);
    ((TH1F*) fListOfHistos->FindObject("eventStatCent"))->Fill((Double_t)1, centrality);
    
    // event mixing

    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.
    
    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.
    
    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.
    
    if (!pool)
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, zVtx));
    if (pool->IsReady()){    
      ((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
      Int_t nMix = pool->GetCurrentNEvents();

      ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
      ((TH1F*) fListOfHistos->FindObject("eventStatCent"))->Fill((Double_t)2, centrality);
      
      // Fill mixed-event histos here  
      for (Int_t jMix=0; jMix<nMix; jMix++){
	TObjArray* bgTracks = pool->GetEvent(jMix);
	
	//standard mixed event
	if(!fThreeParticleMixed){
	  fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kMixedNS, tracksClone, tracksClone, bgTracks, bgTracks, 1.0 / (2*nMix), kFALSE, kFALSE, applyEfficiency);
	  fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kMixedNS, bgTracks, bgTracks, tracksClone, tracksClone, 1.0 / (2*nMix), kFALSE, kFALSE, applyEfficiency);
	}
	
	//1plus1 mixed event
	fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kMixed1plus1, tracksClone, tracksClone, bgTracks, bgTracks, 1.0 / (2*nMix), kTRUE, kFALSE, applyEfficiency);
	fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kMixed1plus1, bgTracks, bgTracks, tracksClone, tracksClone, 1.0 / (2*nMix), kTRUE, kFALSE, applyEfficiency);
	
	//mixed combinatorics
	fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kMixedCombNS, tracksClone, bgTracks, tracksClone, bgTracks, 1.0 / (2*nMix), kFALSE, kFALSE, applyEfficiency);
	fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kMixedCombNS, bgTracks, tracksClone, bgTracks, tracksClone, 1.0 / (2*nMix), kFALSE, kFALSE, applyEfficiency);
	
	//background same from mixed comb
	if(fUseBackgroundSameFromMixedComb){
	  fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kBackgroundSameNS, tracksClone, bgTracks, tracksClone, bgTracks, 1.0 / (2*nMix), kFALSE, kTRUE, applyEfficiency);
	  fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kBackgroundSameNS, bgTracks, tracksClone, bgTracks, tracksClone, 1.0 / (2*nMix), kFALSE, kTRUE, applyEfficiency);
	}

	  //mixed event for background Same
	  fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kMixedBackgroundSameNS, tracksClone, tracksClone, bgTracks, bgTracks, 1.0 / nMix, kFALSE, kTRUE, applyEfficiency);
      }
        
    
    }
    }
    // ownership is with the pool now
    pool->UpdatePool(tracksClone);
  }
}

//________________________________________________________________________
void AliAnalysisTaskTwoPlusOne::Terminate(Option_t*)
{
    //terminate function is called at the end
    //can be used to draw histograms etc.
    
    
}


TObjArray* AliAnalysisTaskTwoPlusOne::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliCFParticle which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* particle = (AliVParticle*) tracks->UncheckedAt(i);
    Double_t part_eta = particle->Eta();
    Double_t part_phi = particle->Phi();

    if(fRandomPosition){
      part_phi=gRandom->Uniform(0, TMath::TwoPi());
      part_eta=gRandom->Uniform(-1*fTrackEtaCut, fTrackEtaCut);
    }

    AliCFParticle* copy = new AliCFParticle(particle->Pt(), part_eta, part_phi, particle->Charge(), 0);
    copy->SetUniqueID(particle->GetUniqueID());
    tracksClone->Add(copy);
  }
  
  return tracksClone;
}

//____________________________________________________________________
void AliAnalysisTaskTwoPlusOne::AddEventCombination(TObjArray* tracks)
{
  //if fEventCombination was full before, clear it
  if(fUsedEvents==4){
    fEventCombination->Clear();
    fUsedEvents = 0;
  }

  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliCFParticle* part = (AliCFParticle*) tracks->UncheckedAt(i);
    fEventCombination->Add(part);
  }
  fUsedEvents++;
}


//____________________________________________________________________
void  AliAnalysisTaskTwoPlusOne::AddSettingsTree()
{
  //Write settings to output list
  TTree *settingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  settingsTree->Branch("fnTracksVertex", &fnTracksVertex,"nTracksVertex/I");
  settingsTree->Branch("fZVertex", &fZVertex,"ZVertex/D");
  //settingsTree->Branch("fCentralityMethod", fCentralityMethod.Data(),"CentralityMethod/C");
  settingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  settingsTree->Branch("fTrackEtaCutMin", &fTrackEtaCutMin, "TrackEtaCutMin/D");
  settingsTree->Branch("fPtMin", &fPtMin, "PtMin/D");
  settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fSharedClusterCut", &fSharedClusterCut,"SharedClusterCut/D");
  settingsTree->Branch("fCrossedRowsCut", &fCrossedRowsCut,"CrossedRowsCut/I");
  settingsTree->Branch("fFoundFractionCut", &fFoundFractionCut,"FoundFractionCut/D");
  settingsTree->Branch("fTrackStatus", &fTrackStatus,"TrackStatus/I");
  settingsTree->Branch("fThreeParticleMixed", &fThreeParticleMixed,"fThreeParticleMixed/I");
  settingsTree->Branch("fMixingTracks", &fMixingTracks,"MixingTracks/I");
  settingsTree->Branch("fMode", &fMode,"Mode/I");
  
  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
}  


//____________________________________________________________________
AliGenEventHeader* AliAnalysisTaskTwoPlusOne::GetFirstHeader()
{
  // get first MC header from either ESD/AOD (including cocktail header if available)
  
  if (fMcEvent)
  {
    // ESD
    AliHeader* header = (AliHeader*) fMcEvent->Header();
    if (!header)
      return 0;
      
    AliGenCocktailEventHeader* cocktailHeader = dynamic_cast<AliGenCocktailEventHeader*> (header->GenEventHeader());
    if (cocktailHeader)
      return dynamic_cast<AliGenEventHeader*> (cocktailHeader->GetHeaders()->First());

    return dynamic_cast<AliGenEventHeader*> (header->GenEventHeader());
  }
}
