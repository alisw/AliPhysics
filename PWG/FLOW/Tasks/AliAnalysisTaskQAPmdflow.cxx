#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSeqCollection.h"
#include "TObjArray.h"
#include "TObjArray.h"
#include "TChain.h"
#include "TMCProcess.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TNtuple.h"

#include "AliLog.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliESDtrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliVertex.h"
#include "AliFlowEventSimple.h"
#include "AliFlowEvent.h"
#include "AliFlowVector.h"
#include "AliESDPmdTrack.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisTaskQAPmdflow.h"

ClassImp(AliAnalysisTaskQAPmdflow)

//________________________________________________________________________
AliAnalysisTaskQAPmdflow::AliAnalysisTaskQAPmdflow()
  : AliAnalysisTaskSE(),
    fOutput(NULL),
    fEventCuts(NULL),
    fRPTrackCuts(NULL),
    fPOITrackCuts(NULL)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskQAPmdflow::AliAnalysisTaskQAPmdflow(const char* name)
  : AliAnalysisTaskSE(name),
    fOutput(NULL),
    fEventCuts(NULL),
    fRPTrackCuts(NULL),
    fPOITrackCuts(NULL)
{
  // Constructor
  DefineOutput(1, TObjArray::Class());
}

//________________________________________________________________________
void AliAnalysisTaskQAPmdflow::UserCreateOutputObjects()
{
  // Called once at the beginning
  fOutput=new TObjArray();
  
  //define histograms
  TH1F* histB = new TH1F("PMD ADC cutB","PMD ADC cutB",500,0,10000);
  fOutput->Add(histB); 

  TH1F* histA = new TH1F("PMD ADC cutA","PMD ADC cutA",500,0,10000);
  fOutput->Add(histA); 

  TH1F* histCelB = new TH1F("PMD ncell CutB", "PMD ncell CutB",100,0,100);
  fOutput->Add(histCelB); 
  TH1F* histCelA = new TH1F("PMD ncell CutA", "PMD ncell CutA",100,0,100);
  fOutput->Add(histCelA);
  
  //post data here as it doesn't change anyway (the pointer to list anyway)

  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskQAPmdflow::UserExec(Option_t *)
{
  //get teh input data
  AliESDEvent* event = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!event)
    {
      AliFatal("no ESD event");
      return;
    }

  fRPTrackCuts->SetEvent(event);

  AliFlowTrackCuts::trackParameterType sourceRP = fRPTrackCuts->GetParamType();
  AliFlowTrackCuts::trackParameterType sourcePOI = fPOITrackCuts->GetParamType();
  Bool_t PmdTrRp = kFALSE;
  Bool_t PmdTrPoi = kFALSE;
  
  if(sourcePOI == 4) PmdTrPoi = kTRUE;
  if(sourceRP  == 4) PmdTrRp = kTRUE;
  if((!PmdTrPoi) && (!PmdTrRp))
    {
      printf("Error : PMD track is not used as POI or RP");
      return;
    }
  
  TH1F* hPmdAdcB = static_cast<TH1F*>(fOutput->At(0));
  TH1F* hPmdAdcA = static_cast<TH1F*>(fOutput->At(1));
  TH1F* hPmdNcelB = static_cast<TH1F*>(fOutput->At(2));
  TH1F* hPmdNcelA = static_cast<TH1F*>(fOutput->At(3));
  
  //Bool_t passevent = fEventCuts->IsSelected(event);
  //Bool_t isSelectedEventSelection = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);

  Int_t ntracks = 0;
  
  if(PmdTrRp) ntracks = fRPTrackCuts->GetNumberOfInputObjects();
  if(PmdTrPoi) ntracks = fPOITrackCuts->GetNumberOfInputObjects();
    
  for (Int_t i=0; i < ntracks; i++)
    {
      Bool_t pass = kFALSE;
      TObject* obj = 0x0;
      if(PmdTrPoi){
	obj = fPOITrackCuts->GetInputObject(i);
	if (!obj) continue;
	pass = fPOITrackCuts->IsSelected(obj,i);
      }
      
      if(PmdTrRp){
	obj = fRPTrackCuts->GetInputObject(i);
	if (!obj) continue;
	pass = fRPTrackCuts->IsSelected(obj,i);
      }
      AliESDPmdTrack* trackpmd = dynamic_cast<AliESDPmdTrack*>(obj);
      if (trackpmd)
	{
	  Int_t   det   = trackpmd->GetDetector();
	  Float_t adc   = trackpmd->GetClusterADC(); 
	  Float_t ncel  = trackpmd->GetClusterCells(); 
	  if(det == 0){
	  hPmdAdcB->Fill(adc); if(pass) hPmdAdcA->Fill(adc);
	  hPmdNcelB->Fill(ncel); if(pass) hPmdNcelA->Fill(ncel); 
	  }
	}
    }
}

//________________________________________________________________________
void AliAnalysisTaskQAPmdflow::Terminate(Option_t *)
{
  //terminate
  
}

//________________________________________________________________________
AliAnalysisTaskQAPmdflow::~AliAnalysisTaskQAPmdflow()
{
  //dtor
  delete fRPTrackCuts;
  delete fPOITrackCuts;
  delete fEventCuts;
}

