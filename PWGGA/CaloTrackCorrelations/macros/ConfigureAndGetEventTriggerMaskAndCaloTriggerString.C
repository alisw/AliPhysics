/// \file ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C
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

#include "AliVEvent.h"

#endif

///
/// Main method 
///
/// The options that can be passed to the macro are:
/// \param trigger       : name of the desired trigger. Special names depending on Calo triggers.
/// \param year          : year
/// \param triggerString : specific trigger name string to select at analysis level
///
UInt_t ConfigureAndGetEventTriggerMaskAndCaloTriggerString
(TString trigger, Int_t year, TString & triggerString)
{
  //printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set event trigger class for %s in year %d\n",
  //       trigger.Data(),year);
  
  triggerString = "";
  
  if( trigger.Contains("INT") || trigger.Contains("MB") || trigger.Contains("default") )
  {
    return ( AliVEvent::kINT7 | AliVEvent::kMB );
  }
  else if(trigger.Contains("EMCAL_L0"))
  {
    triggerString = "EMC";
    return ( AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 );
  }
  else if(trigger.Contains("DCAL_L0"))
  {
    triggerString = "DMC";
    return ( AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 );
  }
  else if(trigger.Contains("EMCAL_L1"))
  {
    if(year > 2012) triggerString = "EG1";
    // before 2013 only one kind of L1 trigger
    
    return ( AliVEvent::kEMCEGA );
  }
  else if(trigger.Contains("DCAL_L1"))
  {
    triggerString = "DG1";
    return ( AliVEvent::kEMCEGA );
  }
  else if(trigger.Contains("EMCAL_L2"))
  {
    triggerString = "EG2";
    return ( AliVEvent::kEMCEGA );
  }
  else if(trigger.Contains("DCAL_L2"))
  {
    triggerString = "DG2";
    return ( AliVEvent::kEMCEGA );
  }
  
  return 0;
}

