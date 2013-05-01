void AddPi0FlowAndDependenciesMC()
{
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C");
  AddMCHandler();

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
  AddTaskEventplane();

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
  AddTaskVZEROEPSelection();
  
  // gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");
  // AddAODPHOSTender();

  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb/AddTaskPHOSPi0FlowMC.C") ;
  AliAnalysisTaskPi0FlowMC * task = AddTaskPHOSPi0FlowMC();
}
