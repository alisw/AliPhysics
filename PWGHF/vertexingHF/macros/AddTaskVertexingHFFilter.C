AliAnalysisTaskSEVertexingHF *AddTaskVertexingHFFilter(TString configPWG3d2h="$ALICE_ROOT/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc.C", Bool_t registerFile=kTRUE)
{
 
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskVertexingHFFilter", "No analysis manager to connect to.");
    return NULL;
  }   


  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/AddTaskVertexingHF.C");
  //  TFile::Cp(gSystem->ExpandPathName(configPWG3d2h.Data()), Form("%s/ConfigVertexingHF.C", train_name.Data()));
  TFile::Cp(gSystem->ExpandPathName(configPWG3d2h.Data()), Form("ConfigVertexingHF.C"));
  AliAnalysisTaskSEVertexingHF *taskvertexingHF = AddTaskVertexingHF();
  // Now we need to keep in sync with the ESD filter
  if (!taskvertexingHF) ::Warning("AddTaskVertexingHFFilter", "AliAnalysisTaskSEVertexingHF cannot run for this train conditions - EXCLUDED");
  
  if(registerFile) mgr->RegisterExtraFile("AliAOD.VertexingHF.root");
  taskvertexingHF->SelectCollisionCandidates(0);

  mgr->AddTask(taskvertexingHF);

  return taskvertexingHF;
}
