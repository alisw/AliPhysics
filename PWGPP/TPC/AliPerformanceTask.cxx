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
// Implementation of the AliPerformanceTask class. It checks reconstruction performance 
// for the reconstructed vs MC particle tracks under several conditions. For real data 
// the control QA histograms are filled.
//
// The comparison output objects deriving from AliPerformanceObject 
// (e.g. AliPerformanceRes, AliPerformanceEff, AliPerformanceDEdx, AliPerformanceDCA ...) 
// are stored in the output file (details in description of these classes).
// 
// Author: J.Otwinowski 01/04/2009 
// Changes by M.Knichel 15/10/2010
// Changes by J.Salzwedel 30/9/2014
//------------------------------------------------------------------------------

#include "iostream"

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TFile.h"
#include "TSystem.h"
#include "TBufferFile.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliVfriendEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDVertex.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"

#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"

#include "AliMCInfo.h"
#include "AliESDRecInfo.h"
#include "AliMCInfoCuts.h"
#include "AliRecInfoCuts.h"
#include "AliPerformanceObject.h"
#include "AliTPCPerformanceSummary.h"
#include "AliPerformanceTPC.h"
#include "AliPerformanceDEdx.h"
#include "AliPerformanceMatch.h"
#include "AliPerformanceTask.h"

#include <AliSysInfo.h>

using namespace std;

ClassImp(AliPerformanceTask)

//_____________________________________________________________________________
AliPerformanceTask::AliPerformanceTask() 
  : AliAnalysisTaskSE("TPCqa")
  , fVEvent(0)
  , fVfriendEvent(0)
  , fMC(0)
  , fOutput(0)
  , fOutputSummary(0)
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
  , fUseVfriend(kFALSE)
  , fUseHLT(kFALSE)
  , fUseTerminate(kTRUE)
  , fUseCentrality(0)
  , fUseOCDB(kTRUE)
  , fUseCentralityBin(0)
{
    fEvents = 0;
    fDebug = 0;
    // Dummy Constructor
  // should not be used
}

//_____________________________________________________________________________
AliPerformanceTask::AliPerformanceTask(const char *name, const char */*title*/)
  : AliAnalysisTaskSE(name)
  , fVEvent(0)
  , fVfriendEvent(0)
  , fMC(0)
  , fOutput(0)
  , fOutputSummary(0)
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
  , fUseVfriend(kFALSE)
  , fUseHLT(kFALSE)
  , fUseTerminate(kTRUE)
  , fUseCentrality(0)
  , fUseOCDB(kTRUE)
  , fUseCentralityBin(0)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(0, TTree::Class());
  DefineOutput(1, TList::Class());

    // create the list for comparison objects
  fCompList = new TList;
  fEvents = 0;
  fDebug = 0;
}

//_____________________________________________________________________________
AliPerformanceTask::~AliPerformanceTask()
{
  if (!(AliAnalysisManager::GetAnalysisManager() && AliAnalysisManager::GetAnalysisManager()->IsProofMode())) {
    if (fOutput)     delete fOutput;    fOutput   = 0; 
    if (fOutputSummary) delete fOutputSummary; fOutputSummary = 0;
    if (fCompList)   delete fCompList;  fCompList = 0; 
  }
}

//_____________________________________________________________________________
Bool_t AliPerformanceTask::AddPerformanceObject(AliPerformanceObject *pObj) 
{

    // add comparison object to the list
  if(pObj == 0) {
    Printf("ERROR: Could not add comparison object");
    return kFALSE;
  }

  // add object to the list
  fCompList->AddLast(pObj);
       
return kTRUE;
}

//_____________________________________________________________________________
void AliPerformanceTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

    // create output list
  fOutput = new TList;
  fOutput->SetOwner();
  fPitList = fOutput->MakeIterator();
  
  // create output list
  //fOutputSummary = new TTree;
  
  // add comparison objects to the output
  AliPerformanceObject *pObj=0;
  Int_t count=0;
  TIterator *pitCompList = fCompList->MakeIterator();
  pitCompList->Reset();
  while(( pObj = (AliPerformanceObject *)pitCompList->Next()) != NULL) {
    fOutput->Add(pObj);
    count++;
  }
    Printf("UserCreateOutputObjects(): Number of output comparison objects: %d \n", count);
  
  PostData(1, fOutput);  
  PostData(0, fOutputSummary);  
}

//_____________________________________________________________________________
void AliPerformanceTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
   
  // Decide whether to use HLT or Offline events
  fEvents++;
  cout <<"Event number "<<fEvents<<endl;
  //if(fDebug) AliSysInfo::AddStamp("memleak",fEvents);
  
// Decide whether to use HLT ESD or Offline ESD/AOD
  if(fUseHLT){
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> 
      (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
      return;
    }
    fVEvent = esdH->GetHLTEvent();
    if(!fVEvent) {
      Printf("ERROR: HLTEvent unavailable from ESDInputHandler");
      return;
    }
  }// end if fUseHLT
  else {
    // Get an offline event
      AliVEventHandler *vH = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
      if (!vH) {
          Printf("ERROR: Could not get VEventHandler");
          return;
      }
      fVEvent = vH->GetEvent();
      if(!fVEvent) { Printf("ERROR: Event not available!"); return;}
  }
  
    if(fUseVfriend) {
    if (fUseHLT)
    {
        AliVEventHandler *vH = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
        AliVEvent *offlineVEvent = vH->GetEvent();
        if(!offlineVEvent) {
            Printf("ERROR: Could not get offline event");
            return;
      }
    }
    else
    {
      fVfriendEvent = fVEvent->FindFriend();
    }
    if(!fVfriendEvent) {
      Printf("ERROR: ESD friends not available");
    }
  } // end if fUseVfriend
  
  if(fUseMCInfo) {
      fMC = MCEvent();
  }  

  if (fUseMCInfo && !fMC) {
    Printf("ERROR: MC event not available");
    return;
  }

  // Process analysis
  Bool_t process = kTRUE;

  // Check for centrality
  if (fUseCentrality) {
    if ( CalculateCentralityBin() != fUseCentralityBin ) {
      process = kFALSE;
    }
  }

  // Process comparison
    if(fEvents==1) fVEvent->InitMagneticField();
    if (process) {
    AliPerformanceObject *pObj=0;
    fPitList->Reset();
    while(( pObj = (AliPerformanceObject *)fPitList->Next()) != NULL) {
      //AliInfo(pObj->GetName());
    pObj->Exec(fMC,fVEvent,fVfriendEvent,fUseMCInfo,fUseVfriend);
    }
  }

    if(fDebug) {
        TBufferFile tempMem(TBuffer::kWrite);
        tempMem.WriteObject(fOutput);
        AliSysInfo::AddStamp("memleak",fEvents,tempMem.Length()/1024.);
    }
  // Post output data.
  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AliPerformanceTask::Terminate(Option_t *) 
{
  // Called once at the end 

  if ( !fUseTerminate )
    return;
  
  // check output data
    fOutputSummary = dynamic_cast<TTree*> (GetOutputData(0));
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutput) {
        Printf("ERROR: AliPerformanceTask::Terminate(): fOutput data not available  ..." );
        return;
   }
    if (fOutputSummary) { delete fOutputSummary; fOutputSummary=0; }      
    AliPerformanceObject* pObj=0;
    AliPerformanceTPC*  pTPC = 0;
    AliPerformanceDEdx* pDEdx = 0;
    AliPerformanceMatch* pMatch = 0;
    AliPerformanceMatch* pPull = 0;
    AliPerformanceMatch* pConstrain = 0;
    TIterator* itOut = fOutput->MakeIterator();
    itOut->Reset();
    while(( pObj = dynamic_cast<AliPerformanceObject*>(itOut->Next())) != NULL) { 
      pObj->AnalyseFinal();
      /*      if (!  pTPC)  {    pTPC = dynamic_cast<AliPerformanceTPC*>(pObj); }
        if (! pDEdx)  {   pDEdx = dynamic_cast<AliPerformanceDEdx*>(pObj); }
        if (! pMatch) {  pMatch = dynamic_cast<AliPerformanceMatch*>(pObj); }
        if ((! pPull) && pMatch ) {  pPull = dynamic_cast<AliPerformanceMatch*>(pObj);*/

	if (!strcmp(pObj->GetName(),"AliPerformanceTPC"))  {    pTPC = dynamic_cast<AliPerformanceTPC*>(pObj); }
        if (!strcmp(pObj->GetName(),"AliPerformanceDEdxTPCInner"))  {   pDEdx = dynamic_cast<AliPerformanceDEdx*>(pObj); }
        if (!strcmp(pObj->GetName(),"AliPerformanceMatchTPCITS")) {  pMatch = dynamic_cast<AliPerformanceMatch*>(pObj); }
        if (!strcmp(pObj->GetName(),"AliPerformanceMatchITSTPC")) {  pPull = dynamic_cast<AliPerformanceMatch*>(pObj);}
        if (!strcmp(pObj->GetName(),"AliPerformanceMatchTPCConstrain")) {  pConstrain = dynamic_cast<AliPerformanceMatch*>(pObj);}
    }
  
   
    if(!fUseOCDB)  { 
      printf("DO NOT USE OCDB \n");
      return;
    }
  
/*    if (! AliCDBManager::Instance()->GetDefaultStorage()) { AliCDBManager::Instance()->SetDefaultStorage("raw://"); }
    TUUID uuid;
    TString tmpFile = gSystem->TempDirectory() + TString("/TPCQASummary.") + uuid.AsString() + TString(".root");
    AliTPCPerformanceSummary::WriteToFile(pTPC, pDEdx, pMatch, pPull, pConstrain, tmpFile.Data());
    TChain* chain = new TChain("tpcQA");
    if(!chain) return;
    chain->Add(tmpFile.Data());
    TTree *tree = chain->CopyTree("1");
    if (chain) { delete chain; chain=0; }
    fOutputSummary = tree;*/
      
//      // Post output data.
//      PostData(0, fOutputSummary);

}

//_____________________________________________________________________________
void AliPerformanceTask::FinishTaskOutput()
{
    // called once at the end of each job (on the workernode)
    //
    // projects THnSparse to TH1,2,3
    
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutput) {
        Printf("ERROR: AliPerformanceTask::FinishTaskOutput(): fOutput data not available  ..." );
        return;
   }

      AliPerformanceObject* pObj=0;
      TIterator* itOut = fOutput->MakeIterator();  
      itOut->Reset();
      while(( pObj = dynamic_cast<AliPerformanceObject*>(itOut->Next())) != NULL) {
          pObj->SetRunNumber(fCurrentRunNumber);
          pObj->Analyse();
      }
    
    if(fDebug) {
        TBufferFile tempMem(TBuffer::kWrite);
        tempMem.WriteObject(fOutput);
        AliSysInfo::AddStamp("memleak",fEvents+1,tempMem.Length()/1024.);
    }

    
     // Post output data.
     PostData(1, fOutput);
}

//_____________________________________________________________________________
Bool_t AliPerformanceTask::Notify()
{
  static Int_t count = 0;
  count++;
  Printf("Processing %d. file", count);

  return kTRUE;
}

//________________________________________________________________________
Int_t AliPerformanceTask::CalculateCentralityBin(){
  // Get Centrality bin

  Int_t centrality = -1;
  Float_t centralityF = -1;

  if (fUseCentrality == 0)
    return centrality;
  
  AliCentrality *eventCentrality = NULL;
  if(fVEvent) eventCentrality = fVEvent->GetCentrality();
  if(!eventCentrality) {
    Printf("ERROR: Could not obtain event centrality");
    return centrality;
  }
  // New : 2010-11-18 JMT 
  if ( fUseCentrality == 1 )
    centralityF = eventCentrality->GetCentralityPercentile("V0M");
  else if ( fUseCentrality == 2 )
    centralityF = eventCentrality->GetCentralityPercentile("CL1");
  else if ( fUseCentrality == 3 )
    centralityF = eventCentrality->GetCentralityPercentile("TRK"); 
  if (centralityF == 0.)
    centralityF = 100.;

  if      ( centralityF >=  0. && centralityF <   5.) centrality =  0;
  else if ( centralityF >=  5. && centralityF <  10.) centrality =  5;
  else if ( centralityF >= 10. && centralityF <  20.) centrality = 10;
  else if ( centralityF >= 20. && centralityF <  30.) centrality = 20;
  else if ( centralityF >= 30. && centralityF <  40.) centrality = 30;
  else if ( centralityF >= 40. && centralityF <  50.) centrality = 40;
  else if ( centralityF >= 50. && centralityF <  60.) centrality = 50;
  else if ( centralityF >= 60. && centralityF <  70.) centrality = 60;
  else if ( centralityF >= 70. && centralityF <  80.) centrality = 70;
  else if ( centralityF >= 80. && centralityF <  90.) centrality = 80;
  else if ( centralityF >= 90. && centralityF <=100.) centrality = 90;
  
  return centrality;
}

Bool_t AliPerformanceTask::ResetOutputData(){

    AliPerformanceObject* pObj=0;
    AliPerformanceTPC*  pTPC = 0;
    AliPerformanceDEdx* pDEdx = 0;
    AliPerformanceMatch* pMatch = 0;
    AliPerformanceMatch* pPull = 0;
    AliPerformanceMatch* pConstrain = 0;
    
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    TIterator* itOut = fOutput->MakeIterator();
    itOut->Reset();
    
    while(( pObj = dynamic_cast<AliPerformanceObject*>(itOut->Next())) != NULL) {
        if (!strcmp(pObj->GetName(),"AliPerformanceTPC"))  {
            pTPC = dynamic_cast<AliPerformanceTPC*>(pObj);
            pTPC->GetTPCClustHisto()->Reset("ICE");
            pTPC->GetTPCEventHisto()->Reset("ICE");
            pTPC->GetTPCTrackHisto()->Reset("ICE");
        }
        if (!strcmp(pObj->GetName(),"AliPerformanceDEdxTPCInner"))  {
            pDEdx = dynamic_cast<AliPerformanceDEdx*>(pObj);
            pDEdx->GetDeDxHisto()->Reset("ICE");
        }
        if (!strcmp(pObj->GetName(),"AliPerformanceMatchTPCITS")) {
            pMatch = dynamic_cast<AliPerformanceMatch*>(pObj);
            pMatch->GetResolHisto()->Reset("ICE");
            pMatch->GetPullHisto()->Reset("ICE");
            pMatch->GetTrackEffHisto()->Reset("ICE");
            pMatch->GetTPCConstrain()->Reset("ICE");

        }
        if (!strcmp(pObj->GetName(),"AliPerformanceMatchITSTPC")) {
            pPull = dynamic_cast<AliPerformanceMatch*>(pObj);
            pPull->GetResolHisto()->Reset("ICE");
            pPull->GetPullHisto()->Reset("ICE");
            pPull->GetTrackEffHisto()->Reset("ICE");
            pPull->GetTPCConstrain()->Reset("ICE");
        }
        if (!strcmp(pObj->GetName(),"AliPerformanceMatchTPCConstrain")) {
            pConstrain = dynamic_cast<AliPerformanceMatch*>(pObj);
            pConstrain->GetResolHisto()->Reset("ICE");
            pConstrain->GetPullHisto()->Reset("ICE");
            pConstrain->GetTrackEffHisto()->Reset("ICE");
            pConstrain->GetTPCConstrain()->Reset("ICE");
        }
    }
    return kFALSE;
}
