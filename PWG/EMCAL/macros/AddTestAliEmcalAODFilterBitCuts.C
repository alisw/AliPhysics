#if !defined(__CLING__) && !defined(__CINT__)
#include <TestAliEmcalAODFilterBitCuts.h>
#endif

PWG::EMCAL::TestAliEmcalAODFilterBitCuts *AddTestAliEmcalAODFilterBitCuts(const char *name){
  return PWG::EMCAL::TestAliEmcalAODFilterBitCuts::AddTestAliEmcalAODFilterBitCuts(name);
}
