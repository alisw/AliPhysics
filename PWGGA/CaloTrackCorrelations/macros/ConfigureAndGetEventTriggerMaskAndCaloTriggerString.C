/// \file ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C
/// \ingroup CaloTrackCorrMacros
/// \brief Configuration of the event trigger
///
/// Configuration macro of the task event trigger, depending on trigger name,
/// specially for EMCal triggers. 
/// The trigger event settings depend on special trigger names:
///  * default: Min bias like triggers, kMB or kINT7
///  * EMCAL_L0, DCAL_L0: EMCal or DCal L0 triggers kEMC7, kEMC8, kEMC1 
///  * EMCAL_L1, DCAL_L1: EMCal or DCal L1 triggers, kEMCEGA, with name EG1 or EGA 
///  * EMCAL_L2, DCAL_L2: EMCal or DCal L1 triggers, kEMCEGA, with name EG2 
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
  printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set event trigger class for %s in year %d\n",
         trigger.Data(),year);
  
  triggerString = "";
  
  UInt_t mask = AliVEvent::kAny;
  
  if( trigger.Contains("MBOrSemiOrCentral") )
  {
    mask = (  AliVEvent::kINT7 | AliVEvent::kMB | AliVEvent::kSemiCentral | AliVEvent::kCentral );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set MB | Semi | Central mask %d\n",mask);
  }
  else if( trigger.Contains("INT") || trigger.Contains("MB") || trigger.Contains("default") )
  {
    mask = ( AliVEvent::kINT7 | AliVEvent::kMB );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set MB/INT7 mask %d\n",mask);
  }  
  else if( trigger.Contains("SemiCentral") )
  {
    mask = ( AliVEvent::kSemiCentral );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set SemiCentral mask %d\n",mask);
  }  
  else if( trigger.Contains("Central") )
  {
    mask = ( AliVEvent::kCentral );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set Central mask %d\n",mask);
  }
  else if(trigger.Contains("EMDCAL_L0"))
  {
    triggerString = "MC";
    mask = ( AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L0 EMC/DMC/1/7/8 mask %d\n",mask);
  }
  else if(trigger.Contains("EMCAL_L0"))
  {
    triggerString = "EMC";
    mask = ( AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L0 EMC1/7/8 mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("DCAL_L0"))
  {
    triggerString = "DMC";
    mask = ( AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L0 EMC1/7/8 mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("EMDCAL_L1"))
  {    
    if(year > 2012) triggerString = "G1";
    else            triggerString = "EGA";

    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("EMCAL_L1"))
  {
    if(year > 2012) triggerString = "EG1";
    else            triggerString = "EGA";
    // before 2013 only one kind of L1 trigger
    
    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("DCAL_L1"))
  {
    triggerString = "DG1";
    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("EMDCAL_L2"))
  {
    triggerString = "G2";
    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("EMCAL_L2"))
  {
    triggerString = "EG2";
    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("DCAL_L2"))
  {
    triggerString = "DG2";
    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("EMDCAL_GA"))
  {
    triggerString = "G";
    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("EMCAL_GA"))
  {
    triggerString = "EG";
    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  else if(trigger.Contains("DCAL_GA"))
  {
    triggerString = "DG";
    mask = ( AliVEvent::kEMCEGA );
    printf("ConfigureAndGetEventTriggerCaloTrackCorr - Set L1 EMCEGA mask %d and string %s\n",mask,triggerString.Data());
  }
  if ( trigger.Contains("CaloOnly") )
  {
    mask = AliVEvent::kCaloOnly;
    printf("\t If periods without TPC, change trigger mask %d!!!!\n",mask);
    
    if ( trigger.Contains("L0") )
      triggerString = "CDMC7PER";
    
    if ( trigger.Contains("MB") )
      triggerString = "CINT7-B-NOPF-CALOPLUS";
  }
  
  return mask;
}

