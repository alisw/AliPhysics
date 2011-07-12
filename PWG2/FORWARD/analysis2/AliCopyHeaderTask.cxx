/**
 * @file   AliCopyHeaderTask.cxx
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:59:32 2011
 * 
 * @brief  Task to copy ESD header to AOD 
 * 
 * @ingroup pwg2_forward_tasks 
 */

#include "AliCopyHeaderTask.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"

ClassImp(AliCopyHeaderTask)
#if 0 
; // for emacs - do not remove 
#endif

void
AliCopyHeaderTask::UserExec(Option_t*)
{
  // 
  // Called at every event 
  //
  // Copies information from ESD header to AOD header
  // 
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(AODEvent());

  if (!esd) { 
    AliWarning("Missing ESD event");
    return;
  }
  if (!aod) { 
    AliWarning("Missing AOD event");
    return;
  }
  
  AliAODHeader* aodHeader = aod->GetHeader();
  if (!aodHeader) { 
    AliWarning("Missing AOD header");
    aodHeader = new AliAODHeader(esd->GetRunNumber(),
				 esd->GetBunchCrossNumber(),
				 esd->GetOrbitNumber(),
				 esd->GetPeriodNumber());
    aod->AddHeader(aodHeader);
  }

  aodHeader->SetRunNumber(esd->GetRunNumber());
  aodHeader->SetOfflineTrigger(fInputHandler->IsEventSelected());
  aodHeader->SetBunchCrossNumber(esd->GetBunchCrossNumber());
  aodHeader->SetOrbitNumber(esd->GetOrbitNumber());
  aodHeader->SetPeriodNumber(esd->GetPeriodNumber());
  aodHeader->SetEventType(esd->GetEventType());
  aodHeader->SetEventNumberESDFile(esd->GetHeader()->GetEventNumberInFile());
  if(esd->GetCentrality())
    aodHeader->SetCentrality(new AliCentrality(*(esd->GetCentrality())));
  else
    aodHeader->SetCentrality(0);

  aodHeader->SetFiredTriggerClasses(esd->GetFiredTriggerClasses());
  aodHeader->SetTriggerMask(esd->GetTriggerMask()); 
  aodHeader->SetTriggerCluster(esd->GetTriggerCluster());
  aodHeader->SetL0TriggerInputs(esd->GetHeader()->GetL0TriggerInputs());    
  aodHeader->SetL1TriggerInputs(esd->GetHeader()->GetL1TriggerInputs());    
  aodHeader->SetL2TriggerInputs(esd->GetHeader()->GetL2TriggerInputs());    
  
  aodHeader->SetMagneticField(esd->GetMagneticField());
  aodHeader->SetMuonMagFieldScale(esd->GetCurrentDip()/6000.);
  aodHeader->SetZDCN1Energy(esd->GetZDCN1Energy());
  aodHeader->SetZDCP1Energy(esd->GetZDCP1Energy());
  aodHeader->SetZDCN2Energy(esd->GetZDCN2Energy());
  aodHeader->SetZDCP2Energy(esd->GetZDCP2Energy());
  aodHeader->SetZDCEMEnergy(esd->GetZDCEMEnergy(0),esd->GetZDCEMEnergy(1));
}

void
AliCopyHeaderTask::Terminate(Option_t*)
{
  // Called at the end of the job  - does nothing 
}

//
// EOF
//
