#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/AliTRDefficiency.h"
#include "PWG1/TRD/AliTRDefficiencyMC.h"
#include "PWG1/TRD/AliTRDmultiplicity.h"
#endif

#include "PWG1/TRD/macros/helper.C"
void AddTRDefficiency(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Int_t map = ParseOptions(trd);
  if(!(TSTBIT(map, kEfficiency))) return;

  Bool_t mc = mgr->GetMCtruthEventHandler();
  AliTRDefficiency *eff = 0x0;
  mgr->AddTask(eff = new AliTRDefficiency((char*)"TRDefficiency"));
  eff->SetDebugLevel(0);
  mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  
  mgr->ConnectInput(eff,  1, ci[0]);
  mgr->ConnectOutput(eff, 1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));
    

  // TRD combined tracking efficiency
  AliTRDefficiencyMC *effMC = 0x0;
  if(mc && TSTBIT(map, kEfficiencyMC)){
    mgr->AddTask(effMC = new AliTRDefficiencyMC((char*)"TRDefficiencyMC"));
    effMC->SetDebugLevel(0);

    // Create containers for input/output
    mgr->ConnectInput(effMC, 0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(effMC, 1, ci[0]);
    mgr->ConnectOutput(effMC, 1, mgr->CreateContainer(effMC->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", effMC->GetName())));
  }

  // TRD single track selection
  if(!(TSTBIT(map, kMultiplicity))) return;

  AliTRDmultiplicity *mult = 0x0;
  mgr->AddTask(mult = new AliTRDmultiplicity((char*)"TRDmultiplicity"));
  mult->SetDebugLevel(0);
  mgr->ConnectInput(mult, 0, mgr->GetCommonInputContainer());  
  mgr->ConnectInput(mult, 1, ci[0]);
  mgr->ConnectOutput(mult, 1, mgr->CreateContainer(mult->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", mult->GetName())));
}

