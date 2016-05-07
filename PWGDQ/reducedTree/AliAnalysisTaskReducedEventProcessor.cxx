/*
 ***********************************************************
 Wrapper class for AliReducedAnalysisTaskSE to be used in the AliAnalysisTask framework*
Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
Instructions in AddTask_EPcorrectionsExample.C
2015/10/01
 *********************************************************
 */

#include <AliAnalysisTaskReducedEventProcessor.h>

#include <iostream>

#include <TROOT.h>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <TChain.h>
#include <THashList.h>
#include <AliInputEventHandler.h>
#include <AliMultiInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliESDEvent.h>
#include "AliReducedBaseEvent.h"
#include "AliReducedEventInfo.h"
#include "AliHistogramManager.h"
#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedEventInputHandler.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskReducedEventProcessor);


//_________________________________________________________________________________
AliAnalysisTaskReducedEventProcessor::AliAnalysisTaskReducedEventProcessor() :
  AliAnalysisTaskSE(),
  fReducedTask(0x0),
  fOutputSlot(),
  fContainerType(),
  fNoutputSlots(),
  fRunningMode(kUseEventsFromTree),
  fReducedEvent()
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskReducedEventProcessor::AliAnalysisTaskReducedEventProcessor(const char* name, Int_t runningMode) :
  AliAnalysisTaskSE(name),
  fReducedTask(0x0),
  fOutputSlot(),
  fContainerType(),
  fNoutputSlots(0),
  fRunningMode(runningMode),
  fReducedEvent()
{
  //
  // Constructor
  //
   if(fRunningMode==kUseOnTheFlyReducedEvents)
      DefineInput(0,AliReducedBaseEvent::Class());
   if(fRunningMode==kUseEventsFromTree)
      DefineInput(0,TChain::Class());
   
  DefineOutput(1,THashList::Class());
}


//______________________________________________________________________________
void AliAnalysisTaskReducedEventProcessor::ConnectInputData(Option_t* /*option*/)
{
   //
   // Special implementation for ConnectInputData() in order to allow using non AliVEvent events
   //

   // Connect input handlers (multi input handler is handled)
   ConnectMultiHandler();
   
   if (fInputHandler && fInputHandler->GetTree()) {
      if(fInputHandler->IsA()==AliESDInputHandler::Class() || fInputHandler->IsA()==AliAODInputHandler::Class()) {
         fInputEvent = fInputHandler->GetEvent();
      }
      if(fInputHandler->IsA()==AliReducedEventInputHandler::Class()) {
        fReducedEvent = ((AliReducedEventInputHandler*)fInputHandler)->GetReducedEvent();
      }
   } else {
      AliError("No Input Event Handler connected") ; 
      return ; 
   }
   // Disconnect multi handler
   DisconnectMultiHandler();
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedEventProcessor::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //
  fReducedTask->GetHistogramManager()->AddHistogramsToOutputList();
  PostData(1, fReducedTask->GetHistogramManager()->GetHistogramOutputList());
                                          
  //for(Int_t i=0; i<fNoutputSlots; i++)   DefineOutput(1, fOutputSlot[i]->Class());
  return;
}


//________________________________________________________________________________________________________
void AliAnalysisTaskReducedEventProcessor::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //   
  AliReducedBaseEvent* event = NULL;
  if(fRunningMode==kUseOnTheFlyReducedEvents) 
     event = dynamic_cast<AliReducedBaseEvent*>(GetInputData(0)); 
  
  if(fRunningMode==kUseEventsFromTree) {
     fInputHandler = (AliInputEventHandler *)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
     fMultiInputHandler = dynamic_cast<AliMultiInputEventHandler *>(fInputHandler);
     if (fMultiInputHandler)
        fInputHandler = dynamic_cast<AliInputEventHandler *>(fMultiInputHandler->GetFirstInputEventHandler());
     
     AliReducedEventInputHandler* handler = dynamic_cast<AliReducedEventInputHandler *>(fInputHandler);
     if(handler)
       event = handler->GetReducedEvent();
  }
  
  if(!event) return;
  
  fReducedTask->SetEvent(event);
  fReducedTask->Process();
  //for(Int_t i=0; i<fNoutputSlots; i++)   if(fContainerType[i]==1) PostData(i, fOutputSlot[i]);
} 


//__________________________________________________________________
void AliAnalysisTaskReducedEventProcessor::FinishTaskOutput()
{
    //
    // Finish Task 
    //
  fReducedTask->Finish();
  //for(Int_t i=0; i<fNoutputSlots; i++)   if(fContainerType[i]==0) PostData(i, fOutputSlot[i]);
  PostData(1, fReducedTask->GetHistogramManager()->GetHistogramOutputList());
  return;
}
