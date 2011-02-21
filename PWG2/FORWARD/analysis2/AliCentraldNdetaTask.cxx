//====================================================================
#include "AliCentraldNdetaTask.h"
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TList.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"


//____________________________________________________________________
TH2D*
AliCentraldNdetaTask::GetHistogram(AliAODEvent* aod, Bool_t mc) 
{
  // Get objects from the event structure 
  TObject* obj = 0;
  if (mc) obj = aod->FindListObject("CentralClustersMC");
  else    obj = aod->FindListObject("CentralClusters");

  // We should have a central object at least 
  if (!obj) { 
    if (!mc) AliWarning("No Central object found AOD");
    return 0;
  }

  // Cast to good types 
  AliAODCentralMult* central   = static_cast<AliAODCentralMult*>(obj);

  return &(central->GetHistogram());
}

//
// EOF
//
