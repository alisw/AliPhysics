void AddPi0FlowAndDependencies()
{
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
  AddTaskEventplane();

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
  AddTaskVZEROEPSelection();
  
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");
  AddAODPHOSTender();

  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb/AddTaskPHOSPi0Flow.C") ;
  AliAnalysisTaskPi0Flow * task = AddTaskPHOSPi0Flow("PHOSPi0Flow", "", AliVEvent::kAny, AliAnalysisTaskPi0Flow::kSemiCentralInclusive);
}
