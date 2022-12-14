AliAnalysisTaskPHOSPbPbQARun2* AddTaskPHOSPbPbQARun2( Bool_t isMC = kFALSE,
                                                      TString tenderOption = "Run2Tune",
   						      Int_t tenderPass = 1,
   						      TString nonlinearity = "Run2Tune",
						      Int_t recoPass = 1,
                                                      TString name = "PbPbQA")
{
  //Add PHOS PbPb QA task to the PWGPP QA train.
  //See PHOSPbPb.C how to run it locally or standalone.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSPbPbQA", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSPbPbQA", "This task requires an input event handler");
    return NULL;
  }

  //Add multitplicity, tender etc.
  TMacro multadd=(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
  AliMultSelectionTask* multSelectionTask = reinterpret_cast<AliMultSelectionTask*>(multadd.Exec(Form("%d", 0)));
 
  AliPHOSTenderTask *tenderPHOS = reinterpret_cast<AliPHOSTenderTask *>
      (gInterpreter->ExecuteMacro(
        Form("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C(\"%s\", \"%s\", \"%s\", %d, %d, \"%s\", %f)", 
                                                                    "PHOSTenderTask","PHOStender",
								    tenderOption.Data(), tenderPass, 
								    isMC, nonlinearity.Data(), 0.020)));

  TMacro addresp(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
  addresp.Exec(Form("%d", isMC));
  addresp.Exec(Form("%d", recoPass));

  AliAnalysisTaskPHOSPbPbQARun2* task = new AliAnalysisTaskPHOSPbPbQARun2("PbPbQA");
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  const char* filename="AnalysisResults_PHOSPbPbQA.root";
  const char* containername=NULL;

  if(!containername) containername = "PHOSPbPbQAResults";
  if(!filename)     filename =  mgr->GetCommonFileName();

  mgr->ConnectOutput(task, 1, mgr->CreateContainer(containername,TList::Class(), AliAnalysisManager::kOutputContainer, filename));
  return task;
}
