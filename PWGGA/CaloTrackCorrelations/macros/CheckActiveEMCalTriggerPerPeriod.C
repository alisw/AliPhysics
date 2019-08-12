/// \file CheckActiveEMCalTriggerPerPeriod.C
/// \ingroup CaloTrackCorrMacros
/// \brief Check EMCal trigger availability in period
///
/// Check if the selected EMCal trigger is appropriate
/// to run the analysis, depending on the period
/// certain triggers were not available.
///
/// The trigger names checked are:
///  * default: Min bias like triggers, kMB or kINT7
///  * EMCAL_L0, DCAL_L0: EMCal or DCal L0 triggers kEMC7, kEMC8, kEMC1 
///  * EMCAL_L1, DCAL_L1: EMCal or DCal L1 triggers, kEMCEGA, with name EG1 or EGA 
///  * EMCAL_L2, DCAL_L2: EMCal or DCal L1 triggers, kEMCEGA, with name EG2 
/// Run MC analysis for no trigger, default/MB option.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TString.h>

#endif

///
/// Main method 
///
/// The options that can be passed to the macro are:
///
/// \param simulation: bool with data (0) or MC (1) condition
/// \param trigger: trigger string name (EMCAL_L0, EMCAL_L1, EMCAL_L2, DCAL_L0, DCAL_L1, DCAL_L2), it can be modified for CaloOnly periods.
/// \param period: LHCXXx
/// \param year: 2011, ...
///
/// \return True if analysis can be done.
///
Bool_t CheckActiveEMCalTriggerPerPeriod(Bool_t simulation, TString & trigger, TString period, Int_t year)
{
  // Accept directly all MB kind of events
  //
  if ( trigger.Contains("default") || trigger.Contains("INT") || trigger.Contains("MB") ) 
    return kTRUE;
  
  // MC analysis has no trigger dependence, execute only for the default case
  //
  if ( simulation )
  {
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : Triggered events not checked in simulation, SKIP trigger %s! \n", 
           trigger.Data());
    
    return kFALSE;
  }
  
  // Triggers introduced in 2011
  //
  if ( year < 2011 && ( trigger.Contains("EMCAL") || trigger.Contains("DCAL") ) )
  {
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No triggered events for year < 2011, SKIP trigger %s! \n", 
           trigger.Data());
    
    return kFALSE;
  }
  
  // DCal Triggers introduced in 2015
  //
  if ( year < 2014 && trigger.Contains("DCAL") )
  {
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No triggered events by DCal for year < 2014, SKIP trigger %s! \n", 
           trigger.Data());
    
    return kFALSE;
  }
  
  // EG2 trigger only activated from 2013
  //
  if ( year  < 2013 && trigger.Contains("L2") )
  { 
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : EG2 trigger not available for year < 2012, SKIP trigger %s in %s \n", 
           trigger.Data(),period.Data());
    
    return kFALSE;
  }
  
  // Triggers only activated in 2013 from LHC13d for physics (it might be there are in b and c but not taking data)
  //
  if ( year == 2013 && trigger.Contains("L") && ( period.Contains("b") || period.Contains("c") ) )
  { 
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : Triggers not available for year 2013 in period %s, SKIP trigger %s! \n",
           period.Data(), trigger.Data());
    
    return kFALSE;
  }
  
  // DCal Triggers introduced in 2015
  //
  if ( year < 2014 && ( trigger.Contains("DCAL") ) )
  {
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No triggered events by DCal for year < 2014, SKIP trigger %s! \n", 
           trigger.Data());
    
    return kFALSE;
  }
  
  // L0 trigger used for periods below LHC11e? 
  //
  if ( period == "LHC11h" && trigger.Contains("EMCAL_L0") )
  {
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No EMCAL_L0 triggered events by EMCal for period LHC11h, SKIP trigger %s! \n", 
           trigger.Data());
    
    return kFALSE;
  }
  
  // L1 trigger not used until LHC11e? period, what about LHC11f?
  //
  if ( period.Contains("LHC11") && period != "LHC11h" && trigger.Contains("EMCAL_L1") )
  {
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal for period %s, SKIP \n", 
           trigger.Data(),period.Data());
    
    return kFALSE;
  }
  
  // L1 trigger not used again until LHC12c period
  //
  if ( ( period == "LHC12a" ||  period == "LHC12b" ) && trigger.Contains("EMCAL_L1") )
  { 
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal for period %s, SKIP \n", 
           trigger.Data(),period.Data());
    
    return kFALSE;
  }
  
  // Run2: No trigger used again until LHC15i period
  //
  if ( year == 2015 && ( period == "LHC15h" ||  period == "LHC15g" || period == "LHC15f" || period == "LHC15e" ||  
                         period == "LHC15d" ||  period == "LHC15c" || period == "LHC15b" || period == "LHC15a"    ) )
  { 
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal for period %s, SKIP \n", 
           trigger.Data(),period.Data());
    
    return kFALSE;
  }
  
  // Run2: L1 trigger not used again until LHC15o period
  //
  if ( year == 2015 && period != "LHC15o" && !trigger.Contains("L0") )
  { 
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal for period %s, SKIP \n", 
           trigger.Data(),period.Data());
    
    return kFALSE;
  }
  
  // Run2: L0 and L2 trigger not used in LHC15o period
  //
  if ( year == 2015 && period == "LHC15o" && ( trigger.Contains("L0") || trigger.Contains("L2") ) )
  { 
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal for period %s, SKIP \n", 
           trigger.Data(),period.Data());
    
    return kFALSE;
  }

  // Run2: Triggers not used in first LHC16 periods 
  //
  if ( year == 2016 && trigger.Contains("L") )
  { 
    if ( period == "LHC16b" || period == "LHC16c" || period == "LHC16d" ||
         period == "LHC16e" || period == "LHC16f" || period == "LHC16g" || 
         period == "LHC16h" ) 
    {
      printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal for period %s, SKIP \n", 
             trigger.Data(),period.Data());
      
      return kFALSE;
    }
  }
  
  // Run2: triggers not used in first LHC17 periods and LHC17n period XeXe
  //
  if ( year == 2017 && trigger.Contains("L") ) 
  { 
    if ( period == "LHC17a" || period == "LHC17b" || period == "LHC17c" ||
         period == "LHC17d" || period == "LHC17e" || period == "LHC17f" || 
         period == "LHC17g" || period == "LHC17n" ) 
    {
      printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal triggers for period %s, SKIP \n", 
             trigger.Data(),period.Data());
      
      return kFALSE;
    }
  }
  
  // Run2: L1 trigger not used in LHC17pq period
  //
  if ( year == 2017 && ( period == "LHC17p" || period == "LHC17q" ) && ( trigger.Contains("L1") ) )
  { 
    printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal L1 trigger for period %s, SKIP \n", 
           trigger.Data(),period.Data());
    
    return kFALSE;
  }
  
  // Some periods with triggered events do not have TPC and the trigger mask is kCaloOnly, 
  // indicate this via the trigger string, so that in macro ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C
  // the proper trigger mask AliVEvent::kCaloOnly is applied.
  if ( !trigger.Contains("CaloOnly") && 
      ( period == "LHC15n" || period == "LHC17p" || period == "LHC17q") ) 
  {
    trigger+="_CaloOnly";
    printf("CheckActiveEMCalTriggerPerPeriod() - Add <_CaloOnly> to trigger string: %s!!!\n",trigger.Data());
  }
  
  
  // Run2: triggers not used in first LHC18 periods
  //
  if ( year == 2018 && trigger.Contains("L") ) 
  { 
    if ( period == "LHC18a" || period == "LHC18b" || period == "LHC18c" ) 
    {
      printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal triggers for period %s, SKIP \n", 
             trigger.Data(),period.Data());
      
      return kFALSE;
    }
  }
 
  // Run2: triggers EG2 and L0 not used in LHC18 Pb-Pb periods
  //
  if ( year == 2018 && (trigger.Contains("L2") || trigger.Contains("L0")) ) 
  { 
    if ( period == "LHC18q" || period == "LHC18r" ) 
    {
      printf("CheckActiveEMCalTriggerPerPeriod() - REMARK! : No %s triggered events by EMCal triggers for period %s, SKIP \n", 
             trigger.Data(),period.Data());
      
      return kFALSE;
    }
  }
  
  // Everything is ok accept this configuration
  return kTRUE;
}



