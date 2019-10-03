// Task to hold TenderSupply in case of running on AOD.
// 
// Author: D.Peressounko after EMCAL Tender task

#include <TChain.h>
#include <TFile.h>

#include "AliAnalysisManager.h"
#include "AliPHOSTenderSupply.h"
#include "AliAODEvent.h"

#include "AliPHOSTenderTask.h"

ClassImp(AliPHOSTenderTask)

//______________________________________________________________________________
AliPHOSTenderTask::AliPHOSTenderTask():
  AliAnalysisTaskSE(),
  fPHOSTender(NULL)
{
  // Default constructor.
}

//______________________________________________________________________________
AliPHOSTenderTask::AliPHOSTenderTask(const char* name):
  AliAnalysisTaskSE(name),
  fPHOSTender(NULL)
{
  // Constructor.
  DefineOutput(1,  AliAODEvent::Class());
}

//______________________________________________________________________________
AliPHOSTenderTask::~AliPHOSTenderTask()
{
  // Destructor

  if (fPHOSTender)
    fPHOSTender->Delete();
}

//______________________________________________________________________________
void AliPHOSTenderTask::SetPHOSTenderSupply(AliPHOSTenderSupply *supply)
{
  // Set tender supply.

  fPHOSTender = supply;
  supply->SetTask(this);
}
   
//______________________________________________________________________________
void AliPHOSTenderTask::ConnectInputData(Option_t *option)
{
  // Connect input data.

  AliAnalysisTaskSE::ConnectInputData(option);
  fPHOSTender->Init();
}

//______________________________________________________________________________
void AliPHOSTenderTask::UserCreateOutputObjects()
{
  // Nothing to be done.
}

//______________________________________________________________________________
void AliPHOSTenderTask::UserExec(Option_t*)
{
  // Process the event.

  fPHOSTender->ProcessEvent();
}
//______________________________________________________________________________
void AliPHOSTenderTask::NotifyRun(){
 //Change of the run number 
  
 fPHOSTender->InitTender(); 
}
