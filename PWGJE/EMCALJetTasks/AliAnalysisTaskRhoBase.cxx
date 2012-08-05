// $Id$
//
// Base class for rho calculation.
// Calculates parameterized rho for given centrality independent of input.
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
  fRho(0),
  fDoCent(0),
  fIsInit(0)
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskRhoBase::AliAnalysisTaskRhoBase(const char *name) :
  AliAnalysisTaskSE(name),
  fRhoName("Rho"),
  fRhoFunction(0x0),
  fCent(-1),
  fRho(0),
  fDoCent(0),
  fIsInit(0)
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
void AliAnalysisTaskRhoBase::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (!fIsInit) {
    ExecOnce();
    fIsInit = 1;
  }

  DetermineCent();

  Double_t rho = GetRhoFactor(fCent);
  fRho->SetVal(rho);
}      

//________________________________________________________________________
void AliAnalysisTaskRhoBase::DetermineCent() 
{
  // Determine centrality.

  fCent = 99; 
  
  if (fDoCent) {
    AliCentrality *centrality = InputEvent()->GetCentrality();
    
    if (centrality)
      fCent = centrality->GetCentralityPercentile("V0M");
    else
      fCent = 99; // probably pp data
    
    if (fCent < 0) {
      AliWarning(Form("%s: Centrality negative: %f, assuming 99", GetName(), fCent));
      fCent = 99;
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskRhoBase::ExecOnce() 
{
  // Initialize some settings that need to be determined in UserExec.

  // add rho to event if not yet there
  if (!(InputEvent()->FindListObject(fRhoName))) {
    InputEvent()->AddObject(fRho);
  } else {
    AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fRhoName.Data()));
    return;
  }

  // determine if centrality should be used
  TString bt(GetBeamType());
  if (bt == "A-A")
    fDoCent = 1;
  else fDoCent = 0;
}

//_____________________________________________________
TString AliAnalysisTaskRhoBase::GetBeamType()
{
  // Get beam type : pp-AA-pA
  // ESDs have it directly, AODs get it from hardcoded run number ranges
  
  AliVEvent *event = InputEvent();
  if (!event) { 
    AliError(Form("%s: Couldn't retrieve event!", GetName()));
    return "";
  }

  TString beamType;

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  if (esd) {
    const AliESDRun *run = esd->GetESDRun();
    beamType = run->GetBeamType();
  } else {
    Int_t runNumber = event->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
	(runNumber >= 166529 && runNumber <= 170593)) {  // LHC11h
      beamType = "A-A";
    } else {
      beamType = "p-p";
    }
  }

  return beamType;    
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
