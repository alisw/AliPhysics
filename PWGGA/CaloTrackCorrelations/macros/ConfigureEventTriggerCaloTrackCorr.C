/// \file ConfigureEventTriggerCaloTrackCorr.C
/// \ingroup CaloTrackCorrMacros
/// \brief Configuration of the event trigger
///
/// Configuration macro of the task event trigger, depending on trigger name,
/// specially for EMCal triggers. 
/// The trigger event settings depend on special trigger names:
///  * default: Min bias like triggers, kMB or kINT7
///  * EMCAL_L0, DCAL_L0: EMCAL or DCAL L0 triggers kEMC7, kEMC8, kEMC1 
///  * EMCAL_L1, DCAL_L1: EMCAL or DCAL L1 triggers, kEMCEGA, with name EG1 or EGA 
///  * EMCAL_L2, DCAL_L2: EMCAL or DCAL L2 triggers, kEMCEGA, with name EG2 
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TString.h>

#include "AliAnalysisTaskCaloTrackCorrelation.h"
#include "AliAnaCaloTrackCorrMaker.h"

#endif

///
/// Main method calling all the configuration
///
/// The options that can be passed to the macro are:
/// \param task : analysis task pointer
/// \param trigger : name of the desired trigger. Special names.
/// \param year: year
///
void ConfigureEventTriggerCaloTrackCorr
(AliAnalysisTaskCaloTrackCorrelation* task, TString trigger, Int_t year)
{
  //printf("ConfigureEventTriggerCaloTrackCorr::Set event trigger class for %s in year %d\n",
  //       trigger.Data(),year);
  
  AliAnaCaloTrackCorrMaker * maker = task->GetAnalysisMaker();
  
  if( trigger.Contains("INT") || trigger.Contains("MB") || trigger.Contains("default") )
  {
    task->SelectCollisionCandidates( AliVEvent::kINT7 | AliVEvent::kMB );
  }
  else if(trigger.Contains("EMCAL_L0"))
  {
    task ->SelectCollisionCandidates( AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 );
    maker->GetReader()->SetFiredTriggerClassName("EMC");
  }
  else if(trigger.Contains("DCAL_L0"))
  {
    task ->SelectCollisionCandidates( AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 );
    maker->GetReader()->SetFiredTriggerClassName("DMC");
  }
  else if(trigger.Contains("EMCAL_L1"))
  {
    task ->SelectCollisionCandidates( AliVEvent::kEMCEGA );
    if(year > 2012) maker->GetReader()->SetFiredTriggerClassName("EG1");
    // before 2013 only one kind of L1 trigger
  }
  else if(trigger.Contains("DCAL_L1"))
  {
    task ->SelectCollisionCandidates( AliVEvent::kEMCEGA );
    maker->GetReader()->SetFiredTriggerClassName("DG1");
  }
  else if(trigger.Contains("EMCAL_L1_Run1"))
  {
    task ->SelectCollisionCandidates( AliVEvent::kEMCEGA );
    //maker->GetReader()->SetFiredTriggerClassName("EGA");
  }
  else if(trigger.Contains("EMCAL_L2"))
  {
    task ->SelectCollisionCandidates( AliVEvent::kEMCEGA );
    maker->GetReader()->SetFiredTriggerClassName("EG2");
  }
  else if(trigger.Contains("DCAL_L2"))
  {
    task ->SelectCollisionCandidates( AliVEvent::kEMCEGA );
    maker->GetReader()->SetFiredTriggerClassName("DG2");
  }
}

