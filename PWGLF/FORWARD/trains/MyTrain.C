#ifndef __CINT__
# include <AliAnalysisManager.h>
#else 
class AliAnalysisManager;
#endif
#include "TrainSetup.C"
class MyTrain : public TrainSetup
{
public:
  MyTrain(const char* name="myTest") : TrainSetup(name) 
  { 
    fOptions.Set("type", "ESD");
  }
  void CreateTasks(AliAnalysisManager* mgr)
  {
    if (!ParUtilities::MakeScriptPAR(fHelper::Mode() != Helper::kLocal, 
				     "MyAnalysis.C", 
				     "STEERBase,ESD,AOD,ANALYSIS,"
				     "OADB,ANALYSISalice"))
      Fatal("CreateTasks", "Failed to create PAR file");
    fHelper->LoadLibrary("MyAnalysis");
    
    Long_t             r = gROOT->ProcessLine("new MyAnalysis(\"test\")");
    AliAnalysisTaskSE* t = reinterpret_cast<AliAnalysisTaskSE*>(r);
    if (!t) Fatal("CreateTasks", "Failed to make task");
    mgr->AddTask(t);
    
    AliAnalysisDataContainer* sums = 
      mgr->CreateContainer("Sums", TList::Class(), 
                           AliAnalysisManager::kOutputContainer,
                           AliAnalysisManager::GetCommonFileName());
    AliAnalysisDataContainer* results = // Needed for output from Terminate
      mgr->CreateContainer("Results", TList::Class(), 
			   AliAnalysisManager::kParamContainer, // Important!
			   AliAnalysisManager::GetCommonFileName());
    
    mgr->ConnectOutput(t, 1, sums);
    mgr->ConnectOutput(t, 2, results);
    mgr->ConnectInput(t, 0, mgr->GetCommonInputContainer());
  }
  void CreateCentralitySelection(Bool_t, AliAnalysisManager*) {}
  AliVEventHandler* CreateOutputHandler(UShort_t type) { return 0; }
  const char* ClassName() const { return "MyTrain"; }
};
//
// EOF
//
