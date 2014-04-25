#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TTree.h>
//#include <TParameter.h>

#include "AliAnalysisTaskTwoPlusOne.h"
#include "AliCFParticle.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliTwoPlusOneContainer.h"
#include "AliUEHist.h"

#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliVParticle.h"
#include "AliCFContainer.h"

#include "AliEventPoolManager.h"


ClassImp( AliAnalysisTaskTwoPlusOne )

//________________________________________________________________________
AliAnalysisTaskTwoPlusOne::AliAnalysisTaskTwoPlusOne(const char *name)
: AliAnalysisTaskSE(name),
  fMixingTracks(10000),
  fAnalyseUE(0x0),
// pointers to UE classes
  fHistos(0x0),
// handlers and events
  fAOD(0x0),
  fPoolMgr(0x0),
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
  fSelectBit(AliVEvent::kMB|AliVEvent::kUserDefined),
  fCustomBinning(),
  fAlpha(0.2)
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
  fAnalyseUE->SetEventSelection(fSelectBit);

  fListOfHistos = new TList();
  fListOfHistos->SetOwner(kTRUE); 

  fHistos = new AliTwoPlusOneContainer("AliTwoPlusOneContainer", fCustomBinning, fAlpha);
  fHistos->GetData()->SetTrackEtaCut(fTrackEtaCut);

  fListOfHistos->Add(fHistos);

  fListOfHistos->Add(new TH1F("eventStat", ";;events", 4, -0.5, 3.5));

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

  centralityObj = fAOD->GetHeader()->GetCentralityP();


 if (centralityObj)
   centrality = centralityObj->GetCentralityPercentile(fCentralityMethod);
 else
   centrality = -1;
      
  if (fAOD){
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

  
  AliInfo(Form("Centrality is %f", centrality));


  // Vertex selection *************************************************
  if(!fAnalyseUE->VertexSelection(InputEvent(), fnTracksVertex, fZVertex)) return;

  // optimization
  if (centrality < 0)
    return;

  TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(InputEvent(), 0, kTRUE, -1, kTRUE);
  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = CloneAndReduceTrackList(tracks);
  delete tracks;

  const AliVVertex* vertex = InputEvent()->GetPrimaryVertex();
  Double_t zVtx = vertex->GetZ();

  fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kSameNS, tracksClone, tracksClone, tracksClone, tracksClone, 1.0);//same event for near and away side

  ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(1);
  
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
  
  AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx);

  if (!pool)
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, zVtx));
  if (pool->IsReady()){    
    Int_t nMix = pool->GetCurrentNEvents();

    ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
    
    // Fill mixed-event histos here  
    for (Int_t jMix=0; jMix<nMix; jMix++){
      TObjArray* bgTracks = pool->GetEvent(jMix);
      
      fHistos->FillCorrelations(centrality, zVtx, AliTwoPlusOneContainer::kMixedNS, tracksClone, bgTracks, tracksClone, bgTracks, 1.0 / nMix);
    }
  }

  // ownership is with the pool now
  pool->UpdatePool(tracksClone);
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
    AliCFParticle* copy = new AliCFParticle(particle->Pt(), particle->Eta(), particle->Phi(), particle->Charge(), 0);
    copy->SetUniqueID(particle->GetUniqueID());
    tracksClone->Add(copy);
  }
  
  return tracksClone;
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
  settingsTree->Branch("fSelectBit", &fSelectBit,"EventSelectionBit/I");
  settingsTree->Branch("fMixingTracks", &fMixingTracks,"MixingTracks/I");
  
  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
}  
