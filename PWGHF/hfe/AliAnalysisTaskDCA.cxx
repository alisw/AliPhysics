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
//
// The analysis task:
// impact parameter resolution and pull study
// for tracks which survivied the particle cuts
// 
// 
// Authors:
//  Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//  Carlo Bombonati <carlo.bombonati@cern.ch>
//
#include <Riostream.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TString.h>

#include <TCanvas.h>

#include "AliAnalysisManager.h"

#include "AliCFManager.h"

#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliPIDResponse.h"

#include "AliHFEpid.h"
#include "AliHFEcuts.h"
#include "AliHFEdca.h"
#include "AliHFEtools.h"

#include "AliAnalysisTaskDCA.h"


//____________________________________________________________
AliAnalysisTaskDCA::AliAnalysisTaskDCA():
  AliAnalysisTaskSE("Impact Parameter Resolution and Pull Analysis")
  , fPlugins(0)
  , fCuts(0x0)
  , fHFEpid(0x0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fCFM(0x0)
  , fDCA(0x0)
  , fNclustersITS(0x0)
  , fMinNprimVtxContrbutor(0x0)
  , fNEvents(0x0)
  , fResidualList(0x0)
  , fPullList(0x0)
  , fDcaList(0x0)
  , fKfDcaList(0x0)
  , fMcVertexList(0x0)
  , fDataDcaList(0x0)
  , fDataVertexList(0x0)
  , fDataPullList(0x0)
  , fMcPidList(0x0) 
  , fDataPidList(0x0)
  , fHfeDcaList(0x0)
  , fHfeDataDcaList(0x0)
  , fOutput(0x0)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliAnalysisTaskDCA::AliAnalysisTaskDCA(const char * name):
  AliAnalysisTaskSE(name)
  , fPlugins(0)
  , fCuts(0x0)
  , fHFEpid(0x0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fCFM(0x0)
  , fDCA(0x0)
  , fNclustersITS(0x0)
  , fMinNprimVtxContrbutor(0x0)
  , fNEvents(0x0)
  , fResidualList(0x0)
  , fPullList(0x0)
  , fDcaList(0x0) 
  , fKfDcaList(0x0)
  , fMcVertexList(0x0)
  , fDataDcaList(0x0)
  , fDataVertexList(0x0)
  , fDataPullList(0x0)
  , fMcPidList(0x0)
  , fDataPidList(0x0)
  , fHfeDcaList(0x0)
  , fHfeDataDcaList(0x0)
  , fOutput(0x0)
{
  //
  // Default constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TH1I::Class());
  DefineOutput(2, TList::Class());


  //-CUTS SETTING-//
  Int_t nMinTPCcluster = 100;
  Float_t maxDcaXY = 0.5;
  Float_t maxDcaZ = 1.0;
  //--------------//
  AliHFEcuts *hfecuts = new AliHFEcuts;
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(nMinTPCcluster);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetMaxImpactParam(maxDcaXY, maxDcaZ);
  SetHFECuts(hfecuts);
 
  fHFEpid = new AliHFEpid("PIDforDCAanalysis");

}

//____________________________________________________________
AliAnalysisTaskDCA::AliAnalysisTaskDCA(const AliAnalysisTaskDCA &ref):
  AliAnalysisTaskSE(ref)
  , fPlugins(ref.fPlugins)
  , fCuts(ref.fCuts)  
  , fHFEpid(ref.fHFEpid)
  , fPIDdetectors(ref.fPIDdetectors)
  , fPIDstrategy(ref.fPIDstrategy)
  , fCFM(ref.fCFM)
  , fDCA(ref.fDCA)
  , fNclustersITS(ref.fNclustersITS)
  , fMinNprimVtxContrbutor(ref.fMinNprimVtxContrbutor)
  , fNEvents(ref.fNEvents)
  , fResidualList(ref.fResidualList)
  , fPullList(ref.fPullList)
  , fDcaList(ref.fDcaList)
  , fKfDcaList(ref.fKfDcaList)
  , fMcVertexList(ref.fMcVertexList)
  , fDataDcaList(ref.fDataDcaList)
  , fDataVertexList(ref.fDataVertexList)
  , fDataPullList(ref.fDataPullList)
  , fMcPidList(ref.fMcPidList)
  , fDataPidList(ref.fDataPidList)
  , fHfeDcaList(ref.fHfeDcaList)
  , fHfeDataDcaList(ref.fHfeDataDcaList)
  , fOutput(ref.fOutput)
{
  //
  // Copy Constructor
  //
  AliInfo("Copy Constructor");
  ref.Copy(*this);
}

//____________________________________________________________
AliAnalysisTaskDCA &AliAnalysisTaskDCA::operator=(const AliAnalysisTaskDCA &ref){
  //
  // Assignment operator
  //
  if(this == &ref) return *this;
  AliAnalysisTaskSE::operator=(ref);
  fPlugins = ref.fPlugins;
  fCuts = ref.fCuts;
  fHFEpid = ref.fHFEpid;
  fPIDdetectors = ref.fPIDdetectors;
  fPIDstrategy = ref.fPIDstrategy;
  fCFM = ref.fCFM;
  fDCA = ref.fDCA;
  fNclustersITS = ref.fNclustersITS;
  fMinNprimVtxContrbutor = ref.fMinNprimVtxContrbutor;
  fNEvents = ref.fNEvents;
  fResidualList = ref.fResidualList;
  fPullList = ref.fPullList;
  fDcaList = ref.fDcaList;
  fKfDcaList = ref.fKfDcaList;
  fMcVertexList = ref.fMcVertexList;
  fDataDcaList = ref.fDataDcaList;
  fDataVertexList = ref.fDataVertexList;
  fDataPullList = ref.fDataPullList;
  fMcPidList = ref.fMcPidList;
  fDataPidList = ref.fDataPidList;
  fHfeDcaList = ref.fHfeDcaList;    
  fHfeDataDcaList = ref.fHfeDataDcaList;
  fOutput = ref.fOutput;

  return *this;
}

//____________________________________________________________
AliAnalysisTaskDCA::~AliAnalysisTaskDCA(){
  //
  // Destructor
  //

  if(fHFEpid) delete fHFEpid;
  if(fCFM) delete fCFM;
  if(fDCA) delete fDCA;  
  if(fNEvents) delete fNEvents;
  if(fResidualList){ 
    fResidualList->Clear();
    delete fResidualList;
  }

  if(fPullList){ 
    fPullList->Clear();
    delete fPullList;
  }
  
  if(fDcaList){ 
    fDcaList->Clear();
    delete fDcaList;
  }
  if(fKfDcaList){ 
    fKfDcaList->Clear();
    delete fKfDcaList;
  }

  if(fMcVertexList){
    fMcVertexList->Clear();
    delete   fMcVertexList;
  }

  if(fDataDcaList){ 
    fDataDcaList->Clear();
    delete fDataDcaList;
  }

  if(fDataVertexList){
    fDataVertexList->Clear();
    delete   fDataVertexList;
  }
  if(fDataPullList){ 
    fDataPullList->Clear();
    delete fDataPullList;
  }

  if(fMcPidList){
    fMcPidList -> Clear();
    delete fMcPidList;
  }
  if(fDataPidList){
    fDataPidList -> Clear();
    delete fDataPidList;
  }

  if(fHfeDcaList) {
    fHfeDcaList->Clear();
    delete fHfeDcaList;
  } 
  
  if(fHfeDataDcaList) {
    fHfeDataDcaList->Clear();
    delete fHfeDataDcaList;
  } 
  
  if(fOutput){ 
    fOutput->Clear();
    delete fOutput;
  }
  
}

//____________________________________________________________
void AliAnalysisTaskDCA::UserCreateOutputObjects(){
  // create output objects
  // fNEvents
  // residual and pull
  //printf("\n=====UserCreateOutputObjects=====\n");
  
  // Automatic determination of the analysis mode
  AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!inputHandler){
    AliError("NoEvent Handler available");
    return;
  }

  if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
    SetAODAnalysis();
  } else {
    SetESDAnalysis();
    if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())
      SetHasMCData();
  }
  

  fNEvents = new TH1I("nEvents", "Number of Events in the Analysis", 5, -0.5, 4.5); // Number of Events neccessary for the analysis and not a QA histogram
  if(!fOutput) fOutput = new TList;
  // Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  MakeParticleContainer();
  // Temporary fix: Initialize particle cuts with 0x0
  for(Int_t istep = 0; istep < fCFM->GetParticleContainer()->GetNStep(); istep++)
    fCFM->SetParticleCutsList(istep, 0x0);
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  
  fCuts->Initialize(fCFM);
  
  if(!fHFEpid) AliWarning("Hello, fHFEpid is not available");
  cout<<"  Hello this is a cout "<<endl<<endl;

  if(GetPlugin(kHFEpid)) {   
    fHFEpid->SetHasMCData(HasMCData());
    fHFEpid->AddDetector("TOF", 0);
    fHFEpid->AddDetector("TPC", 1);
    cout<<endl<<" ---> TPC and TOF added to the PID"<<endl;
		fHFEpid->ConfigureTOF();
    fHFEpid->ConfigureTPCdefaultCut();
    fHFEpid->InitializePID();
  }

  // dca study----------------------------------

  
  if(!fDCA) fDCA = new AliHFEdca;
  if(!fResidualList) fResidualList = new TList();
  if(!fPullList) fPullList = new TList();
  if(!fDcaList) fDcaList = new TList();
  if(!fKfDcaList) fKfDcaList = new TList();
  if(!fMcVertexList) fMcVertexList = new TList();
  if(!fDataDcaList) fDataDcaList = new TList();
  if(!fDataVertexList) fDataVertexList = new TList();
  if(!fDataPullList) fDataPullList = new TList();
  if(!fMcPidList) fMcPidList = new TList();
  if(!fDataPidList) fDataPidList = new TList();  
  
  if(!fHfeDcaList) fHfeDcaList = new TList();
  if(!fHfeDataDcaList) fHfeDataDcaList = new TList();

  if(HasMCData()) {    
    if(GetPlugin(kImpactPar) ) {
      fDCA->CreateHistogramsResidual(fResidualList);
      fDCA->CreateHistogramsPull(fPullList);
      fDCA->CreateHistogramsDca(fDcaList);
      fOutput->AddAt(fResidualList,0);
      fOutput->AddAt(fPullList,1);
      fOutput->AddAt(fDcaList,2);
    } 
    if(GetPlugin(kKFdca)){
      fDCA->CreateHistogramsKfDca(fKfDcaList);
      fOutput->AddAt(fDcaList,3);
    }
    if(GetPlugin(kPrimVtx)){//<---
      fDCA->CreateHistogramsVertex(fMcVertexList);
      fOutput->AddAt(fMcVertexList,4);
    }
    if(GetPlugin(kCombinedPid)){//<---
      fDCA->CreateHistogramsPid(fMcPidList);
      fOutput->AddAt(fMcPidList, 5);
    }
    if(GetPlugin(kHFEpid)){//<---
      fDCA->CreateHistogramsHfeDca(fHfeDcaList);
      fOutput->AddAt(fHfeDcaList, 6);
    }
  } // mc case

  if(!HasMCData())  { 
    
    if(GetPlugin(kPrimVtx)){
      fDCA->CreateHistogramsDataVertex(fDataVertexList);  
      fOutput->AddAt(fDataVertexList,0);
    }    

    if(GetPlugin(kCombinedPid)){
      fDCA->CreateHistogramsDataDca(fDataDcaList);  
      fDCA->CreateHistogramsDataPull(fDataPullList);  
      fDCA->CreateHistogramsDataPid(fDataPidList);
      fOutput->AddAt(fDataDcaList,1);
      fOutput->AddAt(fDataPullList,2);
      fOutput->AddAt(fDataPidList, 3);
    }
    if(GetPlugin(kHFEpid)){
      fDCA->CreateHistogramsHfeDataDca(fHfeDataDcaList);
      fOutput->AddAt(fHfeDataDcaList, 4);
    }
    


  }  // data case
  
}

//____________________________________________________________
void AliAnalysisTaskDCA::UserExec(Option_t *){
  //
  // Run the analysis
  // 
  //printf("\n=====UserExec=====\n");
  if(HasMCData()) printf("WITH MC!\n");

  AliDebug(3, "Processing ESD events");

  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return;
  }
  if(HasMCData()){
    AliDebug(4, Form("MC Event: %p", fMCEvent));
    if(!fMCEvent){
      AliError("No MC Event, but MC Data required");
      return;
    }
  }

  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }
 

  // protection
  if(IsESDanalysis() && HasMCData()){
    // Protect against missing MC trees
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH){
      AliError("No MC Event Handler available");
      return;
    }
    if(!mcH->InitOk()) return;
    if(!mcH->TreeK()) return;
    if(!mcH->TreeTR()) return;
  }

  
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
    pidResponse = AliHFEtools::GetDefaultPID(HasMCData(), fInputEvent->IsA() == AliAODEvent::Class());
  }
  fHFEpid->SetPIDResponse(pidResponse);
  ProcessDcaAnalysis();

  PostData(1, fNEvents);
  PostData(2, fOutput);
}
//____________________________________________________________
void AliAnalysisTaskDCA::ProcessDcaAnalysis(){

  //printf("\n=====ProcessDcaAnalysis=====\n");

  //
  // Loop ESD
  //
  AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!fESD){
    AliError("ESD Event required for ESD Analysis");
      return;
  }

  AliMCEvent *fMC = 0x0;
  if(HasMCData()){
    fMC = dynamic_cast<AliMCEvent*>(fMCEvent);
    if(!fMC){
      AliError("MC Event required for Analysis");
      return;
    }
  }

  fNEvents->Fill(1);  // original event number before cut
  fDCA->ApplyExtraCuts(fESD,fMinNprimVtxContrbutor);  // cut on primVtx contributors
  fNEvents->Fill(3);  // events number after cut
  fCFM->SetRecEventInfo(fESD);
 
  // event cut level
  if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) return;

  AliESDtrack *track = 0x0;  
  AliMCParticle *mctrack = 0x0;
  AliESDVertex *vtxESDSkip = 0x0;
  AliHFEpidObject hfetrack;

  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    
    track = fESD->GetTrack(itrack);
    if(HasMCData()) mctrack = dynamic_cast<AliMCParticle *>(fMC->GetTrack(TMath::Abs(track->GetLabel())));

    // RecPrim: primary cuts
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
    // RecKine: ITSTPC cuts  
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
    // HFEcuts: ITS layers cuts
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
    
    if(track->GetITSclusters(0)<=fNclustersITS) continue;  // require number of ITS clusters
    
    // track accepted, do PID
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    if(HasMCData()) hfetrack.SetMCTrack(mctrack);

    //printf("Track %d passed all the cuts!\n",itrack);

    if(HasMCData()){
      if(GetPlugin(kPrimVtx))
	fDCA->FillHistogramsVtx(fESD, fMC);
      if(GetPlugin(kImpactPar)) 
	fDCA->FillHistogramsDca(fESD, track, fMC);
      if(GetPlugin(kKFdca)) 
	fDCA->FillHistogramsKfDca(fESD, track, fMC);
      if(GetPlugin(kCombinedPid)) 
	fDCA->FillHistogramsPid(track, fMC);
      if(GetPlugin(kHFEpid)) { // data-like
	if(fHFEpid->IsSelected(&hfetrack)){ 
	  
	  //	  printf("Found an electron in p+p collision! from HFE pid \n");
	  if(!vtxESDSkip){
	    // method from Andrea D 28.05.2010
	    AliVertexerTracks *vertexer = new AliVertexerTracks(fESD->GetMagneticField());
	    vertexer->SetITSMode();
	    vertexer->SetMinClusters(fNclustersITS);
	    Int_t skipped[2];
	    skipped[0] = (Int_t)track->GetID();
	    vertexer->SetSkipTracks(1,skipped);
	    vtxESDSkip = (AliESDVertex*)vertexer->FindPrimaryVertex(fESD);
	    delete vertexer; vertexer = NULL;
	    if(vtxESDSkip->GetNContributors()<fMinNprimVtxContrbutor) continue;
	  }
	  //printf("\n[ABOUT TO FILL HFE DCA: MC!]\n");
	  fDCA->FillHistogramsHfeDataDca(fESD, track, vtxESDSkip); 
	}
      } // plugin for hfepid 
    }  // MC

    if(!HasMCData()){
      if(GetPlugin(kPrimVtx))
	fDCA->FillHistogramsDataVtx(fESD);
      if(GetPlugin(kCombinedPid)) {

	// method from Andrea D 28.05.2010
	AliVertexerTracks *vertexer = new AliVertexerTracks(fESD->GetMagneticField());
	vertexer->SetITSMode();
	vertexer->SetMinClusters(fNclustersITS);
	Int_t skipped[2];
	skipped[0] = (Int_t)track->GetID();
	vertexer->SetSkipTracks(1,skipped);
	vtxESDSkip = (AliESDVertex*)vertexer->FindPrimaryVertex(fESD);
	delete vertexer; vertexer = NULL;
	if(vtxESDSkip->GetNContributors()<fMinNprimVtxContrbutor) continue;

	fDCA->FillHistogramsDataDca(fESD, track, vtxESDSkip);
	fDCA->FillHistogramsDataPid(track);
      }
      if(GetPlugin(kHFEpid)) {
	if(fHFEpid->IsSelected(&hfetrack)) {
	  //	  printf("Found an electron in p+p collision! from HFE pid \n");
	  if(!vtxESDSkip){
	    // method from Andrea D 28.05.2010
	    AliVertexerTracks *vertexer = new AliVertexerTracks(fESD->GetMagneticField());
	    vertexer->SetITSMode();
	    vertexer->SetMinClusters(fNclustersITS);
	    Int_t skipped[2];
	    skipped[0] = (Int_t)track->GetID();
	    vertexer->SetSkipTracks(1,skipped);
	    vtxESDSkip = (AliESDVertex*)vertexer->FindPrimaryVertex(fESD);
	    delete vertexer; vertexer = NULL;
	    if(vtxESDSkip->GetNContributors()<fMinNprimVtxContrbutor) continue;
	  }
	printf("\n[ABOUT TO FILL HFE DCA: DATA!]\n");
	  fDCA->FillHistogramsHfeDataDca(fESD, track,vtxESDSkip);    
	} 
      } // plugin for hfepid
    }  // data case

  } // track loop

}


//____________________________________________________________
void AliAnalysisTaskDCA::Terminate(Option_t *){
  //
  // Terminate not implemented at the moment
  //  
  //printf("\n=====Terminate=====\n");
  
  if(GetPlugin(kPostProcess)){
    fOutput = dynamic_cast<TList *>(GetOutputData(1));
    if(!fOutput){
      AliError("Results not available");
      return;
    }
    PostProcess();
  }
  
}


//____________________________________________________________
void AliAnalysisTaskDCA::Load(TString filename){

  //printf("\n=====Load=====\n");

  // no need for postprocessing for the moment
  TFile *input = TFile::Open(filename.Data());
  if(!input || input->IsZombie()){
    AliError("Cannot read file");
    return;
  }

  input->Close();
  delete input;
  
  
}

//____________________________________________________________
void AliAnalysisTaskDCA::PostProcess(){
  // do post processing
  // should do fitting here for dca resolution
  // moved to an external macro to do the job
  
  //printf("\n=====PostProcess=====\n");
  Load("HFEdca.root");
  TCanvas *c1 = new TCanvas("c1", "number of analyzed events", 300, 400);
  fNEvents->Draw();
  c1->SaveAs("temp.png");
 
}




//____________________________________________________________
void AliAnalysisTaskDCA::PrintStatus() const {
  
  //
  // Print Analysis status
  //
  printf("\n\tAnalysis Settings\n\t========================================\n");
  printf("\t Running on %s\n", !HasMCData()?"p+p collision data":"MC sample");
  printf("\t Cuts: %s\n", (fCuts != NULL) ? "YES" : "NO");
  printf("\t Impact parameter analysis is %s\n", GetPlugin(kImpactPar)?"ON":"OFF");
  printf("\t Using AliKFParticle for analysis? %s\n", GetPlugin(kKFdca)?"ON":"OFF");
  printf("\t Primary vertex analysis is %s\n", GetPlugin(kPrimVtx)?"ON":"OFF");
  printf("\t Combined pid analysis is %s\n", GetPlugin(kCombinedPid)?"ON":"OFF");
  printf("\t HFE pid analysis is %s\n", GetPlugin(kHFEpid)?"ON":"OFF");
  printf("\t Post process analysis is %s\n", GetPlugin(kPostProcess)?"ON":"OFF");
  printf("\t ");
  printf("\n");
}

//__________________________________________                                                  
void AliAnalysisTaskDCA::SwitchOnPlugin(Int_t plug){
  //                                            
  // Switch on Plugin          
  // Available:                                  
  //  - analyze impact parameter
  //  - Post Processing  

  AliDebug(2,Form("SwitchOnPlugin %d",plug));  

  switch(plug){
  case kPostProcess: 
    SETBIT(fPlugins, plug); 
    break;
  case kImpactPar: 
    SETBIT(fPlugins, plug); 
    break;
  case kPrimVtx: 
    SETBIT(fPlugins, plug); 
    break;
  case kCombinedPid:
    SETBIT(fPlugins, plug); 
    break;
  case kHFEpid:
    SETBIT(fPlugins, plug); 
    break;
  case kKFdca:
    SETBIT(fPlugins, plug); 
    break;
  default: 
    AliError("Unknown Plugin");
  };
}


//____________________________________________________________
void AliAnalysisTaskDCA::MakeParticleContainer(){

  //printf("\n=====MakeParticleContainer=====\n");
  //
  // Create the particle container (borrowed from AliAnalysisTaskHFE)
  //
  const Int_t kNvar   = 3 ; //number of variables on the grid:pt,eta, phi
  const Double_t kPtmin = 0.1, kPtmax = 10.;
  const Double_t kEtamin = -0.9, kEtamax = 0.9;
  const Double_t kPhimin = 0., kPhimax = 2. * TMath::Pi();

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 40; //bins in pt
  iBin[1] =  8; //bins in eta 
  iBin[2] = 18; // bins in phi

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    binEdges[ivar] = new Double_t[iBin[ivar] + 1];

  //values for bin lower bounds
  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)TMath::Power(10,TMath::Log10(kPtmin) + (TMath::Log10(kPtmax)-TMath::Log10(kPtmin))/iBin[0]*(Double_t)i);  
  for(Int_t i=0; i<=iBin[1]; i++) binEdges[1][i]=(Double_t)kEtamin  + (kEtamax-kEtamin)/iBin[1]*(Double_t)i;
  for(Int_t i=0; i<=iBin[2]; i++) binEdges[2][i]=(Double_t)kPhimin  + (kPhimax-kPhimin)/iBin[2]*(Double_t)i;

  //one "container" for MC
  const Int_t kNcutStepsESDtrack = AliHFEcuts::kNcutStepsRecTrack + 1;
  const Int_t kNcutStepsTrack = AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack;
  AliCFContainer* container = new AliCFContainer("container","container for tracks", (kNcutStepsTrack + 1 + 2*(kNcutStepsESDtrack + 1)), kNvar, iBin);

  //setting the bin limits
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    container -> SetBinLimits(ivar, binEdges[ivar]);
  fCFM->SetParticleContainer(container);

  //create correlation matrix for unfolding
  Int_t thnDim[2*kNvar];
  for (int k=0; k<kNvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = iBin[k];
    thnDim[k+kNvar] = iBin[k];
  }


}

//____________________________________________________________
void AliAnalysisTaskDCA::AddPIDdetector(TString detector){
  
  //
  // Adding PID detector to the task
  //
  //printf("\n=====AddPIDdetector=====\n");
  
  if(!fPIDdetectors.Length()) 
    fPIDdetectors = detector;
  else
    fPIDdetectors += ":" + detector;
}

