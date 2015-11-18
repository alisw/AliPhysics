/*
 ***********************************************************
 Wrapper class for AliReducedAnalysisTaskSE to be used in the AliAnalysisTask framework*
Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
Instructions in AddTask_EPcorrectionsExample.C
2015/10/01
 *********************************************************
 */

#include <AliAnalysisTaskReducedEventProcessor.h>

//#include "AliSysInfo.h"
#include <iostream>

#include <TROOT.h>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <TChain.h>
#include <THashList.h>
#include <AliInputEventHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliESDEvent.h>
#include "AliReducedEventInfo.h"
#include "AliHistogramManager.h"
#include "AliReducedAnalysisTaskSE.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskReducedEventProcessor);


//_________________________________________________________________________________
AliAnalysisTaskReducedEventProcessor::AliAnalysisTaskReducedEventProcessor() :
  AliAnalysisTaskSE(),
  fReducedTask(0x0),
  fOutputSlot(),
  fContainerType(),
  fNoutputSlots()
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskReducedEventProcessor::AliAnalysisTaskReducedEventProcessor(const char* name) :
  AliAnalysisTaskSE(name),
  fReducedTask(0x0),
  fOutputSlot(),
  fContainerType(),
  fNoutputSlots(0)
{
  //
  // Constructor
  //
  //DefineInput(0,TChain::Class());
  DefineInput(0,AliReducedEventInfo::Class());
  DefineOutput(1,THashList::Class());


}


//_________________________________________________________________________________
void AliAnalysisTaskReducedEventProcessor::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //


  //DefineInput(0,TChain::Class());

  fReducedTask->GetHistogramManager()->AddHistogramsToOutputList();
  PostData(1, fReducedTask->GetHistogramManager()->GetHistogramOutputList());
                                          
  //for(Int_t i=0; i<fNoutputSlots; i++)   DefineOutput(1, fOutputSlot[i]->Class());
  //fReducedTask->Init();
  //PostData(1,*(fReducedTask->GetHistogramManager()->GetOutputHistogramList()));

  return;
}



//________________________________________________________________________________________________________
void AliAnalysisTaskReducedEventProcessor::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //

  //((TH1F*)fReducedTask->GetHistogramManager()->GetHistogram("Event_NoCuts","VtxX"))->Fill(0.);
  AliReducedEventInfo* event = dynamic_cast<AliReducedEventInfo*>(GetInputData(0)); 
  if(!event) return;
  fReducedTask->SetEvent(event);
  fReducedTask->Process();

  //PostData(1, fReducedTask->GetHistogramManager()->GetHistogramOutputList());
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

