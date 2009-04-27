#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TRD/qaRec/macros/AliTRDperformanceTrain.h"
#include "TRD/qaRec/AliTRDefficiency.h"
#include "TRD/qaRec/AliTRDefficiencyMC.h"
#include "TRD/qaRec/AliTRDmultiplicity.h"
#endif


void AddTRDefficiency(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Int_t map = ParseOptions(trd);
  if(!(TSTBIT(map, kEfficiency))) return;

  Bool_t mc = mgr->GetMCtruthEventHandler();
  AliTRDefficiency *eff = 0x0;
  mgr->AddTask(eff = new AliTRDefficiency());
  eff->SetDebugLevel(0);
  mgr->ConnectInput(eff, 0, ci[0]);
  mgr->ConnectOutput(eff, 0, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));
    

  // TRD combined tracking efficiency
  AliTRDefficiencyMC *effMC = 0x0;
  if(mc && TSTBIT(map, kEfficiencyMC)){
    mgr->AddTask(effMC = new AliTRDefficiencyMC());
    effMC->SetDebugLevel(0);

    // Create containers for input/output
    mgr->ConnectInput(effMC, 0, ci[0]);
    mgr->ConnectOutput(effMC, 0, mgr->CreateContainer(effMC->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", effMC->GetName())));
  }

  // TRD single track selection
  if(!(TSTBIT(map, kMultiplicity))) return;

  AliTRDmultiplicity *mult = 0x0;
  mgr->AddTask(mult = new AliTRDmultiplicity());
  mult->SetDebugLevel(0);
  mgr->ConnectInput(mult, 0, ci[0]);
  mgr->ConnectOutput(mult, 0, mgr->CreateContainer(mult->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", mult->GetName())));
}

