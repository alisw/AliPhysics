/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//*****************************************************
//   Class AliEventShape
//   author: Antonio Ortiz Velasquez
//   antonio.ortiz@nucleares.unam.mx
//*****************************************************

#include "AliAnaTransverseEventShapeTask.h"

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>
// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include "AliPPVsMultUtils.h"

#include "AliCentrality.h" 
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 

#include "AliTransverseEventShape.h"

// STL includes
#include <iostream>
using namespace std;

ClassImp(AliAnaTransverseEventShapeTask)
//_____________________________________________________________________________
AliAnaTransverseEventShapeTask::AliAnaTransverseEventShapeTask():
AliAnalysisTaskSE(),
  fESD(0x0),
  fAOD(0x0),
  fPPVsMultUtils(0),
  fMC(0x0),
  fMCStack(0x0),
  fMCArray(0x0),
    
  fUseHybrid(0x0),
  fTrackFilterHybrid1(0x0),
  fTrackFilterHybrid2(0x0),
  fTrackFilterESA(0x0),
  fMinMultESA(0x0),
  fSizeStepESA(0x0),
  fIsAbsEtaESA(0x0),
  fEtaMaxCutESA(0x0),
  fEtaMinCutESA(0x0),
  fPtMaxCutESA(0x0),
  fPtMinCutESA(0x0),
  
  fCentEst("V0M"),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),
  fAnalysisPbPb(kFALSE),
  ftrigBit(0x0),
  fRandom(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),  
  fEtaCut(0.9),  
  fMinCent(0.0),
  fMaxCent(100.0),
  fTriggeredEventMB(-999),
  fVtxStatus(-999),
  fZvtx(-999),
  fZvtxMC(-999),
  fRun(-999),
  fEventId(-999),
  fListOfObjects(0), 
  hVtxBeforeCuts(0x0), 
  hVtxAfterCuts(0x0),
  hn1(0x0),
  hso(0x0),
  hst(0x0),
  HMultRef(0x0),
  HStMultRef0(0x0),
  HSoMultRef0(0x0),
  fESASelection(0x0)

{
  //default constructor
  
  
}


AliAnaTransverseEventShapeTask::AliAnaTransverseEventShapeTask(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fAOD(0x0),
  fPPVsMultUtils(0),
  fMC(0x0),
  fMCStack(0x0),
  fMCArray(0x0),
  
  fUseHybrid(0x0),
  fTrackFilterHybrid1(0x0),
  fTrackFilterHybrid2(0x0),
  fTrackFilterESA(0x0),
  fMinMultESA(0x0),
  fSizeStepESA(0x0),
  fIsAbsEtaESA(0x0),
  fEtaMaxCutESA(0x0),
  fEtaMinCutESA(0x0),
  fPtMaxCutESA(0x0),
  fPtMinCutESA(0x0),
  
  fCentEst("V0M"),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),
  fAnalysisPbPb(kFALSE),
  ftrigBit(0x0),
  fRandom(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),  
  fEtaCut(0.9),  
  fMinCent(0.0),
  fMaxCent(100.0),
  fStoreMcIn(kFALSE),//
  fMcProcessType(-999),
  fTriggeredEventMB(-999),
  fVtxStatus(-999),
  fZvtx(-999),
  fZvtxMC(-999),
  fRun(-999),
  fEventId(-999),
  fListOfObjects(0), 
  hVtxBeforeCuts(0x0),
  hVtxAfterCuts(0x0),
  hn1(0x0),
  hso(0x0),
  hst(0x0),
  HMultRef(0x0),
  HStMultRef0(0x0),
  HSoMultRef0(0x0),
  fESASelection(0x0)
{


  // Default constructor (should not be used)
    
  DefineOutput(1, TList::Class());
}

AliAnaTransverseEventShapeTask::~AliAnaTransverseEventShapeTask() {
  //
  // Destructor
  //
  
  if (fESASelection) {
    delete fESASelection;
    fESASelection = 0x0;
  }  

}
//______________________________________________________________________________
void AliAnaTransverseEventShapeTask::UserCreateOutputObjects()
{ 
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested 
  // We also create the random generator here so it might get different seeds...
  fRandom = new TRandom(0); // 0 means random seed
 
  //Helper
  if(! fESASelection ) {
    fESASelection = new AliTransverseEventShape();
    fESASelection->SetUseHybridESA(fUseHybrid);
    fESASelection->SetTrackFilterESAHyb1(fTrackFilterHybrid1);
    fESASelection->SetTrackFilterESAHyb2(fTrackFilterHybrid2);
    fESASelection->SetTrackFilterESA(fTrackFilterESA);
    fESASelection->SetMinMultForESA(fMinMultESA);
    fESASelection->SetStepSizeESA(fSizeStepESA);
    fESASelection->SetIsEtaAbsESA(fIsAbsEtaESA);
    fESASelection->SetTrackEtaMinESA(fEtaMinCutESA);
    fESASelection->SetTrackEtaMaxESA(fEtaMaxCutESA);
    fESASelection->SetTrackPtMinESA(fPtMinCutESA);
    fESASelection->SetTrackPtMaxESA(fPtMaxCutESA);
    fESASelection->Init();

  }
  
  //OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();
  
  //
  // Histograms
  hn1=new TH1D("hn1","hn1",11,-1,10);
  fListOfObjects->Add(hn1);
  
  hVtxBeforeCuts = new TH1D("hVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(hVtxBeforeCuts);  
  hVtxAfterCuts = new TH1D("hVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(hVtxAfterCuts);
  
  hso = new TH1D("hso",";spherocity; entries",200,-1.0,1.0);
  fListOfObjects->Add(hso);
  hst = new TH1D("hst",";sphericity; entries",200,-1.0,1.0);
  fListOfObjects->Add(hst);
  HSoMultRef0 = new TH1D("HSoMultRef0",";sphericity; entries",200,-1.0,1.0);
  fListOfObjects->Add(HSoMultRef0);
  HStMultRef0 = new TH1D("HStMultRef0",";sphericity; entries",200,-1.0,1.0);
  fListOfObjects->Add(HStMultRef0); 
  HMultRef= new TH1D("HMultRef",";Multiplicity; entries",200,0.0,200);
  fListOfObjects->Add(HMultRef);
  
  // Post output data.
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnaTransverseEventShapeTask::UserExec(Option_t *) 
{
  // Main loop
  
  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }
  
  
  Bool_t isPileup = kFALSE;
  
  if (fAnalysisType == "ESD"){
    fESD = dynamic_cast<AliESDEvent*>(event);
    if(!fESD){
      Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
    
    isPileup = fESD->IsPileupFromSPD();
    if(fPileUpRej)
      if(isPileup)
	return;
  } else {
    fAOD = dynamic_cast<AliAODEvent*>(event);
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
    
    isPileup = fAOD->IsPileupFromSPD();
    if(fPileUpRej)
      if(isPileup)
	return;    
  }
  
  if (fAnalysisMC) {
    
    if (fAnalysisType == "ESD"){
      fMC = dynamic_cast<AliMCEvent*>(MCEvent());
      if(!fMC){
	Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
      }    
      
      fMCStack = fMC->Stack();
      
      if(!fMCStack){
	Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
      }   
      
      
    } else { // AOD
      
      fMC = dynamic_cast<AliMCEvent*>(MCEvent());
      if(fMC)
	fMC->Dump();
      
      fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
      if(!fMCArray){
	Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
      }    
    }
  }
  
  
  // Get trigger decision
  fTriggeredEventMB = 0; //init
  hn1->Fill(0);
  
  if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ftrigBit ){
    fTriggeredEventMB = 1;  //event triggered as minimum bias
  }
  
  // real data that are not triggered we skip
  if(!fAnalysisMC && !fTriggeredEventMB)
    return; 
  hn1->Fill(1);
  

  if (fAnalysisType == "ESD"){
    const AliESDVertex *vtxESD = fESD->GetPrimaryVertexTracks();
    if(vtxESD->GetNContributors()<1) {
      // SPD vertex
      vtxESD = fESD->GetPrimaryVertexSPD();
      /* quality checks on SPD-vertex */
      if (vtxESD->IsFromVertexerZ() && (vtxESD->GetDispersion() > 0.04 || vtxESD->GetZRes() > 0.25))  
	fZvtx  = -1599; //vertex = 0x0; //
      else if (vtxESD->GetNContributors()<1) 
	fZvtx  = -999; //vertex = 0x0; //
      else
	fZvtx = vtxESD->GetZ();
    }  
    else
      fZvtx = vtxESD->GetZ();
    
  }
  else // AOD
    fZvtx = GetVertex(fAOD);
  hVtxBeforeCuts->Fill(fZvtx);
  
  //cut on the z position of vertex
  if (TMath::Abs(fZvtx) > fVtxCut) {	
    return;
  }
  hn1->Fill(2);
  
  Int_t TrackMult03=0;   
  Double_t spherocity=0;
  Double_t sphericity=0;
  if(fTriggeredEventMB) {    // only analyze triggered events
    TrackMult03=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.3); //fine     
    HMultRef->Fill(TrackMult03);    

    spherocity = fESASelection->GetEventShape( event, "SO", kTRUE );
    sphericity = fESASelection->GetEventShape( event, "ST", kTRUE );
    cout<<"not cuts: So="<<spherocity<<"\t St = "<< sphericity <<endl;
    hso->Fill(spherocity);
    hst->Fill(spherocity);
    
    if(0 < TrackMult03 && TrackMult03 <= 20){
      cout<<"  cut 1: So="<< spherocity <<"\t St = "<< sphericity <<endl;
      HSoMultRef0->Fill(spherocity);
      HStMultRef0->Fill(sphericity);      
    }  
  }
  
  hVtxAfterCuts->Fill(fZvtx);
  
  // Post output data.
  PostData(1, fListOfObjects);
}

//_____________________________________________________________________________
Float_t AliAnaTransverseEventShapeTask::GetVertex(const AliVEvent* event) const
{
  Float_t zvtx = -999;
  
  const AliVVertex* primaryVertex = event->GetPrimaryVertex();   
  if(primaryVertex->GetNContributors()>0)
    zvtx = primaryVertex->GetZ();
  
  return zvtx;
}
//____________________________________
Float_t AliAnaTransverseEventShapeTask::GetTest(){
  
  return 10.0;
  
}
//____________________________________________________________________________
void AliAnaTransverseEventShapeTask::Terminate(Option_t *)
{  
  TFile* fout = new TFile("esa_qa.root", "RECREATE"); 
  if (fESASelection)
    {
      fESASelection->SaveHistos();
    }
  
  fout->Write();
  fout->Close();
  
  Printf("Writing result to esa_qa.root");
}



