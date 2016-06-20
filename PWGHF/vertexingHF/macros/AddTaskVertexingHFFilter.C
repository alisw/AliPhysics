AliAnalysisTaskSEVertexingHF *AddTaskVertexingHFFilter(Int_t collisionSystem, TString configPWGHFd2h="$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_pPbRHF.C", TString localdir="", Int_t runnumber=-1, TString strPeriod="", const char* fname="AliAOD.VertexingHF.root", Bool_t registerFile=kTRUE)
{
 
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskVertexingHFFilter", "No analysis manager to connect to.");
    return NULL;
  }   


  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskVertexingHF.C");
  //  TFile::Cp(gSystem->ExpandPathName(configPWGHFd2h.Data()), Form("%s/ConfigVertexingHF.C", train_name.Data()));
  //  TFile::Cp(gSystem->ExpandPathName(configPWGHFd2h.Data()), Form("ConfigVertexingHF.C"));
  //  TFile::Cp(gSystem->ExpandPathName(configPWGHFd2h.Data()), Form("%s/ConfigVertexingHF.C", localdir.Data()));

  AliAnalysisTaskSEVertexingHF *taskvertexingHF = AddTaskVertexingHF(collisionSystem,localdir,configPWGHFd2h,runnumber,strPeriod,fname);
  // Now we need to keep in sync with the ESD filter
  if (!taskvertexingHF) ::Warning("AddTaskVertexingHFFilter", "AliAnalysisTaskSEVertexingHF cannot run for this train conditions - EXCLUDED");
  
  if(registerFile) mgr->RegisterExtraFile("AliAOD.VertexingHF.root");
  taskvertexingHF->SelectCollisionCandidates(0);

  mgr->AddTask(taskvertexingHF);

  return taskvertexingHF;
}
