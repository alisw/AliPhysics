/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//------------------------------------------------------------------------------
// AliPtResolAnalysisPbPb class. 
// 
// a. functionality:
// - fills analysis control histograms
//
// b. data members:
// - control histograms
//
// Author: J.Otwinowski 04/11/2008 
//------------------------------------------------------------------------------

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "THnSparse.h"

#include "AliHeader.h"  
#include "AliInputEventHandler.h"  
#include "AliAnalysisManager.h"  
#include "AliGenEventHeader.h"  
#include "AliStack.h"  
#include "AliESDEvent.h"  
#include "AliMCEvent.h"  
#include "AliESDtrackCuts.h"  
#include "AliLog.h" 
#include "AliMultiplicity.h"
#include "AliTracker.h"
#include "AliCentrality.h"

#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AliPtResolAnalysisPbPb.h"


using namespace std;

ClassImp(AliPtResolAnalysisPbPb)

//_____________________________________________________________________________
  AliPtResolAnalysisPbPb::AliPtResolAnalysisPbPb(): AlidNdPt(),
  fAnalysisFolder(0),
  fTrackParamHist(0),
  fTrackParamHist2(0),
  fCentralityEstimator(0)
{
  // default constructor
  Init();
}

//_____________________________________________________________________________
AliPtResolAnalysisPbPb::AliPtResolAnalysisPbPb(Char_t* name, Char_t* title): AlidNdPt(name,title),
  fAnalysisFolder(0),
  fTrackParamHist(0),
  fTrackParamHist2(0),
  fCentralityEstimator(0)
{
  Init();
}

//_____________________________________________________________________________
AliPtResolAnalysisPbPb::~AliPtResolAnalysisPbPb() {
  //
  // destructor
  //
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
  if(fTrackParamHist) delete fTrackParamHist; fTrackParamHist=0;
  if(fTrackParamHist2) delete fTrackParamHist2; fTrackParamHist2=0;
}

//_____________________________________________________________________________
void AliPtResolAnalysisPbPb::Init(){
  //
  // Generic histograms to be corrected
  //
  //1/pT:#sigma(1/pT):centrality
  Int_t binsTrackParamHist[3]={400,300,11};
  Double_t minTrackParamHist[3]={0,0,0}; 
  Double_t maxTrackParamHist[3]={1,0.015,100};

  Double_t centrBins[12] = {0.0,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100}; 


  fTrackParamHist = new THnSparseF("fTrackParamHist","1/pT:#sigma(1/pT):centrality",3,binsTrackParamHist,minTrackParamHist,maxTrackParamHist);
  fTrackParamHist->SetBinEdges(2,centrBins);
  fTrackParamHist->GetAxis(0)->SetTitle("1/pT (GeV/c)^{-1}");
  fTrackParamHist->GetAxis(1)->SetTitle("#sigma(1/pT)");
  fTrackParamHist->GetAxis(2)->SetTitle("centrality");
  fTrackParamHist->Sumw2();
  
  //pt:sigma(1/pT)*pT:centrality
  const Int_t ptNbins = 73;
  Double_t bins[74] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 };

  Int_t binsTrackParamHist2[3]={ptNbins,200,11};
  Double_t minTrackParamHist2[3]={0,0,0}; 
  Double_t maxTrackParamHist2[3]={100,0.2,100};

  fTrackParamHist2 = new THnSparseF("fTrackParamHist2","pT:#sigma(1/pT)*pT:centrality",3,binsTrackParamHist2,minTrackParamHist2,maxTrackParamHist2);
  fTrackParamHist2->SetBinEdges(0,bins);
  fTrackParamHist2->GetAxis(0)->SetTitle("pT (GeV/c)");
  fTrackParamHist2->GetAxis(1)->SetTitle("#sigma(1/pT)*pT");
  fTrackParamHist2->GetAxis(2)->SetTitle("centrality");
  fTrackParamHist2->Sumw2();

  // init folder
  fAnalysisFolder = CreateFolder("folderdNdPt","Analysis dNdPt Folder");
  
}

//_____________________________________________________________________________
void AliPtResolAnalysisPbPb::Process(AliESDEvent *const esdEvent, AliMCEvent *const mcEvent)
{
  //
  // Process real and/or simulated events
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }

  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;
    //SetPhysicsTriggerSelection(physicsSelection);

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }


   // centrality determination
   Float_t centralityF = -1.;
   AliCentrality *esdCentrality = esdEvent->GetCentrality();
   centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());

  // get reconstructed vertex  
  const AliESDVertex* vtxESD = 0; 
  Bool_t isRecVertex = kFALSE;
  if(evtCuts->IsRecVertexRequired()) 
  {
    Bool_t bRedoTPCVertex = evtCuts->IsRedoTPCVertex();
    Bool_t bUseConstraints = evtCuts->IsUseBeamSpotConstraint();
    vtxESD = AlidNdPtHelper::GetVertex(esdEvent,evtCuts,accCuts,esdTrackCuts,GetAnalysisMode(),kFALSE,bRedoTPCVertex,bUseConstraints); 
    isRecVertex = AlidNdPtHelper::TestRecVertex(vtxESD, esdEvent->GetPrimaryVertexSPD(), GetAnalysisMode(), kFALSE);
  }

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD) && isRecVertex; 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());

  //
  TObjArray *allChargedTracks=0;
  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    // get all charged tracks
    allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,GetAnalysisMode());
    if(!allChargedTracks) return;

    Int_t entries = allChargedTracks->GetEntries();
    //printf("entries %d \n",entries);

    // fill histograms
    for(Int_t i=0; i<entries;++i) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
      if(!track) continue;
      if(track->Charge()==0) continue;

      // only postive charged 
      if(GetParticleMode() == AlidNdPtHelper::kPlus && track->Charge() < 0) 
        continue;
      
      // only negative charged 
      if(GetParticleMode() == AlidNdPtHelper::kMinus && track->Charge() > 0) 
        continue;

      if(esdTrackCuts->AcceptTrack(track)) 
      {
        if(accCuts->AcceptTrack(track)) 
        {
	  //Double_t x, par[5], cov[15];
	  //track->GetExternalParameters(x, p);
	  //track->GetExternalCovariance(cov);

	  Double_t v[3] = {track->OneOverPt(),TMath::Sqrt(track->GetSigma1Pt2()),centralityF};
	  fTrackParamHist->Fill(v);

	  Double_t v2[3] = {track->Pt(),track->Pt()*TMath::Sqrt(track->GetSigma1Pt2()),centralityF};
	  fTrackParamHist2->Fill(v2);
        }
      }  
    }
  }
 }
}

//_____________________________________________________________________________
Long64_t AliPtResolAnalysisPbPb::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  //
  //TList *collPhysSelection = new TList;

  // collection of generated histograms

  Int_t count=0;
  while((obj = iter->Next()) != 0) {
    AliPtResolAnalysisPbPb* entry = dynamic_cast<AliPtResolAnalysisPbPb*>(obj);
    if (entry == 0) continue; 
    
    //
    fTrackParamHist->Add(entry->fTrackParamHist);
    fTrackParamHist2->Add(entry->fTrackParamHist2);
  }

return count;
}

//_____________________________________________________________________________
void AliPtResolAnalysisPbPb::Analyse() 
{
  // Analyse histograms
  //
  TH1::AddDirectory(kFALSE);
  TObjArray *aFolderObj = new TObjArray;
  if(!aFolderObj) return;
  
  //
  // Reconstructed event vertex
  //
  
  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);
  if(!fAnalysisFolder) { 
    if(aFolderObj) delete aFolderObj;
    return;
  }

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliPtResolAnalysisPbPb::ExportToFolder(TObjArray * const array) 
{
  // recreate folder avery time and export objects to new one
  //
  if(!array) return NULL;

  AliPtResolAnalysisPbPb * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
     // get name and title from old folder
     name = folder->GetName();  
     title = folder->GetTitle();  

	 // delete old one
     delete folder;

	 // create new one
     newFolder = CreateFolder(name.Data(),title.Data());
     newFolder->SetOwner();

	 // add objects to folder
     while(i < size) {
	   newFolder->Add(array->At(i));
	   i++;
	 }
  }

return newFolder;
}

//_____________________________________________________________________________
TFolder* AliPtResolAnalysisPbPb::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
