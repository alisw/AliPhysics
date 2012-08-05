// $Id$
//
// Task to hold TenderSupply in case of running on AOD.
// 
// Author: S.Aiola, C.Loizides

#include <TChain.h>
#include <TFile.h>

#include "AliAnalysisManager.h"
#include "TenderSupplies/AliEMCALTenderSupply.h"
#include "AliAODEvent.h"

#include "AliEmcalTenderTask.h"

ClassImp(AliEmcalTenderTask)

//______________________________________________________________________________
AliEmcalTenderTask::AliEmcalTenderTask():
  AliAnalysisTaskSE(),
  fEMCALTender(NULL)
{
  // Default constructor.
}

//______________________________________________________________________________
AliEmcalTenderTask::AliEmcalTenderTask(const char* name):
  AliAnalysisTaskSE(name),
  fEMCALTender(NULL)
{
  // Constructor.
  DefineOutput(1,  AliAODEvent::Class());
}

//______________________________________________________________________________
AliEmcalTenderTask::~AliEmcalTenderTask()
{
  // Destructor

  if (fEMCALTender)
    fEMCALTender->Delete();
}

//______________________________________________________________________________
void AliEmcalTenderTask::SetEMCALTenderSupply(AliEMCALTenderSupply *supply)
{
  // Set tender supply.

  fEMCALTender = supply;
  supply->SetTask(this);
}
   
//______________________________________________________________________________
void AliEmcalTenderTask::ConnectInputData(Option_t *option)
{
  // Connect input data.

  AliAnalysisTaskSE::ConnectInputData(option);
  fEMCALTender->Init();
}

//______________________________________________________________________________
void AliEmcalTenderTask::UserCreateOutputObjects()
{
  // Nothing to be done.
}

//______________________________________________________________________________
void AliEmcalTenderTask::UserExec(Option_t*)
{
  // Process the event.

  fEMCALTender->ProcessEvent();
}
