AliAnalysisTask* AddTaskQAtrainPIDqa(){
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  return AddTaskPIDqa();
}
