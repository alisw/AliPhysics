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
  printf("AddTRDefficiency <- [0]=\"%s\" [1]=\"%s\" [2]=\"%s\"\n", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName());

  AliTRDrecoTask *eff(NULL);
  mgr->AddTask(eff = new AliTRDefficiency((char*)"TRDefficiency"));
  eff->SetDebugLevel(0);
  //AliLog::SetClassDebugLevel("AliTRDefficiency", 5);  
  mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  
  mgr->ConnectInput(eff, 1, ci[0]);
  mgr->ConnectOutput(eff,1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance", mgr->GetCommonFileName())));
    

  // TRD combined tracking efficiency
  if(mgr->GetMCtruthEventHandler() && TSTBIT(map, kEfficiencyMC)) {
    mgr->AddTask(eff = new AliTRDefficiencyMC((char*)"TRDefficiencyMC"));
    eff->SetDebugLevel(0);
    //AliLog::SetClassDebugLevel("AliTRDefficiencyMC", 5);  

    // Create containers for input/output
    mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(eff, 1, ci[0]);
    mgr->ConnectOutput(eff, 1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance", mgr->GetCommonFileName(),eff->GetName())));
  }

  // TRD single track selection
  if(!(TSTBIT(map, kMultiplicity))) return;

  mgr->AddTask(eff = new AliTRDmultiplicity((char*)"TRDmultiplicity"));
  eff->SetDebugLevel(0);
  //AliLog::SetClassDebugLevel("AliTRDmultiplicity", 5);  
  mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  
  mgr->ConnectInput(eff, 1, ci[0]);
  mgr->ConnectOutput(eff, 1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.Calib%s", mgr->GetCommonFileName(),eff->GetName())));
}

