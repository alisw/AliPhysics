#ifndef __CINT__
#include "AliAnalysisTaskEmcalJetEnergyScale.h"
#endif

EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale *AddTaskEmcalJetEnergyScale(AliJetContainer::EJetType_t jettype, double jetradius, const char *trigger){
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale(jettype, jetradius, trigger);
}