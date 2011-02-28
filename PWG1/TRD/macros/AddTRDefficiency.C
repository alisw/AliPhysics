#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TError.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/AliTRDpwg1Helper.h"
#include "PWG1/TRD/AliTRDefficiency.h"
#include "PWG1/TRD/AliTRDefficiencyMC.h"
#include "PWG1/TRD/AliTRDmultiplicity.h"
#endif

void AddTRDefficiency(AliAnalysisManager *mgr, Int_t map, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Info("AddTRDefficiency", Form("[0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName(), ci[3]->GetName()));

  AliAnalysisDataContainer *evInfoContainer = ci[3];
  AliTRDrecoTask *eff(NULL);
  mgr->AddTask(eff = new AliTRDefficiency((char*)"TRDefficiency"));
  eff->SetDebugLevel(0);
  //AliLog::SetClassDebugLevel("AliTRDefficiency", 5);  
  Int_t trackStatus = 0; // barrel tracks
//                    = 1; // kink tracks
//                    = 2; // SA tracks
  mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  // connect main (ESD) container
  mgr->ConnectInput(eff, 1, ci[trackStatus]);                 // conect track info container
  mgr->ConnectInput(eff, 2, evInfoContainer);                 // conect event info container
  mgr->ConnectOutput(eff,1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance", mgr->GetCommonFileName())));
    

  // TRD combined tracking efficiency
  if(mgr->GetMCtruthEventHandler() && TESTBIT(map, AliTRDpwg1Helper::kEfficiencyMC)) {
    mgr->AddTask(eff = new AliTRDefficiencyMC((char*)"TRDefficiencyMC"));
    eff->SetDebugLevel(0);
    //AliLog::SetClassDebugLevel("AliTRDefficiencyMC", 5);  

    // Create containers for input/output
    mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(eff, 1, ci[trackStatus]);
    mgr->ConnectInput(eff, 2, evInfoContainer);
    mgr->ConnectOutput(eff, 1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance", mgr->GetCommonFileName())));
  }

  // TRD single track selection
  if(!(TESTBIT(map, AliTRDpwg1Helper::kMultiplicity))) return;

  mgr->AddTask(eff = new AliTRDmultiplicity((char*)"TRDmultiplicity"));
  eff->SetDebugLevel(0);
  //AliLog::SetClassDebugLevel("AliTRDmultiplicity", 5);  
  mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  
  mgr->ConnectInput(eff, 1, ci[0]);
  mgr->ConnectOutput(eff, 1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
}

