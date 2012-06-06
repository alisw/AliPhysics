// $Id$
//
// Base class for rho calculation
//
// Author: S.Aiola

#include <TF1.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"

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

  fRho = new AliRhoParameter(fRhoName, 0);
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

//_____________________________________________________
TString AliAnalysisTaskRhoBase::GetBeamType()
{
  // Get beam type : pp-AA-pA
  // ESDs have it directly, AODs get it from hardcoded run number ranges
  
  AliVEvent *event = InputEvent();

  if (!event) { 
    AliError("Couldn't retrieve event!");
    return "";
  }

  TString beamType;

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  if (esd) {
    const AliESDRun *run = esd->GetESDRun();
    beamType = run->GetBeamType();
  }
  else
  {
    Int_t runNumber = event->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
	(runNumber >= 166529 && runNumber <= 170593))    // LHC11h
    {
      beamType = "A-A";
    }
    else 
    {
      beamType = "p-p";
    }
  }

  return beamType;    
}

//________________________________________________________________________
void AliAnalysisTaskRhoBase::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // add rho to event if not yet there
  if (!(InputEvent()->FindListObject(fRhoName))) {
    InputEvent()->AddObject(fRho);
  }

  // determine centrality 
  fCent = 99; 
  
  if (GetBeamType() == "A-A") {
    AliCentrality *centrality = InputEvent()->GetCentrality();
    
    if (centrality)
      fCent = centrality->GetCentralityPercentile("V0M");
    else
      fCent = 99; // probably pp data
    
    if (fCent < 0) {
      AliWarning(Form("Centrality negative: %f, assuming 99", fCent));
      fCent = 99;
    }
  }

  Double_t rhochem = GetRhoFactor(fCent);
  fRho->SetVal(rhochem);
}      

//________________________________________________________________________
void AliAnalysisTaskRhoBase::Terminate(Option_t *) 
{
  // Run at the end of the task.
}
