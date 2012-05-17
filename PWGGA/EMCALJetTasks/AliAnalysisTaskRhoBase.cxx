// $Id$
//
// Base class for rho calculation
//
// Author: A.Saiola

#include <TF1.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliVCluster.h"

#include "AliAnalysisTaskRhoBase.h"

ClassImp(AliAnalysisTaskRhoBase)

//________________________________________________________________________
AliAnalysisTaskRhoBase::AliAnalysisTaskRhoBase() : 
  AliAnalysisTaskSE(),
  fRhoName("Rho"),
  fRhoFunction(0x0),
  fCent(-1),
  fRho(0)
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskRhoBase::AliAnalysisTaskRhoBase(const char *name) :
  AliAnalysisTaskSE(name),
  fRhoName("Rho"),
  fRhoFunction(0x0),
  fCent(-1),
  fRho(0)
{
  // Constructor.
}

//________________________________________________________________________
void AliAnalysisTaskRhoBase::UserCreateOutputObjects()
{
  // Run at beginning of task.

  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!handler) {
    AliError("Input handler not available!");
    return;
  }

  fRho = new TParameter<Double_t>(fRhoName, 0);
}

//________________________________________________________________________
Double_t AliAnalysisTaskRhoBase::GetRhoFactor(Double_t cent)
{
  // Return rho per centrality.

  Double_t rho = -1;
  if (fRhoFunction)
    rho = fRhoFunction->Eval(cent);
  return rho;
}

//________________________________________________________________________
void AliAnalysisTaskRhoBase::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // add rho to event if not yet there
  if (!(InputEvent()->FindListObject(fRhoName))) {
    new(fRho) TParameter<Double_t>(fRhoName, 0);
    InputEvent()->AddObject(fRho);
  }

  // get centrality 
  AliCentrality *centrality = InputEvent()->GetCentrality() ;
  if (centrality)
    fCent = centrality->GetCentralityPercentile("V0M");
  else
    fCent = 99; // probably pp data
  if (fCent < 0) {
    AliError(Form("Centrality negative: %f", fCent));
    return;
  }

  Double_t rhochem = GetRhoFactor(fCent);
  fRho->SetVal(rhochem);
}      

//________________________________________________________________________
void AliAnalysisTaskRhoBase::Terminate(Option_t *) 
{
  // Run at the end of the task.
}
