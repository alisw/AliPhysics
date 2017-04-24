#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TString.h"
#include "TSystem.h"
#include "TGrid.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisAlien.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskSE.h"
#include "TROOT.h"
#include "AliPhysicsSelectionTask.h"
#include "AliOADBPhysicsSelection.h"
#include "AliMultSelectionTask.h"

#endif


Int_t nRuns = 1;



void runGridEsd(
                Int_t ntrk = 1, // cut to define offline trigger. default (1) is INEL>0
                const char * gridMode = "full", // Try with test first
                TString dataDir = "/alice/data/2016/LHC16h",
                TString dataPattern = "/muon_calo_pass1/*/AliESDs.root",
                TString workingDir = "lhc16h",
                TString runListString = "254476",
                Bool_t isMC = 0,
                const char * addTaskString = "AddTaskTrigEff.C(\"outfile.root\",%d,%d)"
                )
{
  
  gSystem->AddIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  if (!TGrid::Connect("alien://")) return;


  TObjArray * runs = runListString.Tokenize(",");
  
  Int_t nruns = runs->GetEntries();
  
  Int_t * runList = new Int_t[nruns];
  for(Int_t iruns = 0; iruns < nruns; iruns++){
    runList[iruns] = ((TObjString*)runs->At(iruns))->GetString().Atoi();
    //    std::cout << " " << iruns << " " << runList[iruns] << std::endl;    
  }
  delete runs;
  
  AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
  AliESDInputHandler* handler = new AliESDInputHandler();
  mgr->SetInputEventHandler(handler);
  



  // ### PHYSICS SELECTIOPN ###
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(isMC);
  if(!physSelTask) { Printf("no physSelTask"); return; }
  AliOADBPhysicsSelection * oadbCustom = new AliOADBPhysicsSelection("oadbCustom");
  Int_t triggerCount = 0;
  // WARNING: INT11 hooked to INT1 (we have no INT1 in run 2)
  oadbCustom->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT11-B-NOPF-CENTNOTRD","B",triggerCount);
  oadbCustom->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C || ADA || ADC");
  oadbCustom->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");
  triggerCount++;
  oadbCustom->AddCollisionTriggerClass(AliVEvent::kINT7,"+CINT7-[I|B|S]-NOPF-[ALL|CENT][NOTRD|]","B",triggerCount);
  oadbCustom->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbCustom->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");


  // ### MULTIPLICITY SELECTIOPN ###
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask * multSelTask =  AddTaskMultSelection();
  multSelTask->SetSelectedTriggerClass(AliVEvent::kAny);
  multSelTask->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/m/mfloris/lhc16h/OADB-LHC16h-Michele.root");
  
  // create task
   AliAnalysisTaskSE * task = 0;
   {
     TString buf1, buf2;
     buf1.Form(".x %s", addTaskString);
     buf2.Form(buf1.Data(), isMC, ntrk);
     std::cout << buf2.Data() << std::endl;
     task = (AliAnalysisTaskSE *)gROOT->ProcessLine( buf2.Data() );
   }

  task->Print();
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(gridMode);
  plugin->SetNtestFiles(2);
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20170409-1");
  plugin->SetGridDataDir(dataDir.Data());
  plugin->SetDataPattern(dataPattern.Data());
  plugin->SetGridWorkingDir(workingDir.Data());
  plugin->SetRunPrefix("000");
  for (Int_t i=0;i<nRuns;i++)  plugin->AddRunNumber(runList[i]);
  plugin->SetGridOutputDir("output");
  //  plugin->SetAnalysisSource("AliAnalysisTaskTest.cxx MyCluster.cxx");
  plugin->SetAnalysisSource("AliAnalysisTaskTrigEff.cxx");
  plugin->SetAdditionalLibs("AliAnalysisTaskTrigEff.h AliAnalysisTaskTrigEff.cxx");
  plugin->SetNrunsPerMaster(1);
  plugin->SetSplitMaxInputFileNumber(100);
  plugin->AddIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  mgr->SetGridHandler(plugin);
  


  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");
 
  //  TChain *chain = new TChain("esdTree");
  //  for (Int_t i=1;i<=1;i++) chain->AddFile(Form("/data/esd/LHC12h/189616/AliESDs.%03i.root",i));
  //  mgr->StartAnalysis("local",chain);
}
