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

ClassImp(AliCentraldNdetaTask)
#ifdef DOXY_INPUT
;
#endif

//____________________________________________________________________
AliCentraldNdetaTask::AliCentraldNdetaTask()
  : AliBasedNdetaTask() 
{ 
  DGUARD(fDebug,3,"Default CTOR of AliCentraldNdetaTask");
}

//____________________________________________________________________
AliCentraldNdetaTask::AliCentraldNdetaTask(const char*)
  : AliBasedNdetaTask("Central") 
{ 
  DGUARD(fDebug,3,"Named CTOR of AliCentraldNdetaTask");
  fCorrEmpty  = false;
  // SetTitle("Central");
}

//____________________________________________________________________
TH2D*
AliCentraldNdetaTask::GetHistogram(const AliAODEvent& aod, Bool_t mc) 
{
  // Get objects from the event structure 
  DGUARD(fDebug,2,"Get our histogram for AliCentraldNdetaTask");
  // Cast to good types 
  AliAODCentralMult* central   = GetCentral(aod, mc, !mc);
  if (!central) return 0;
  return &(central->GetHistogram());
}

//
// EOF
//
