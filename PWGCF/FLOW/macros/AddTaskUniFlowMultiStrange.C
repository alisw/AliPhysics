///////////////////////////////////////////////////////////////////
//                                                               //
// AddUniFlow                                                     //
// Author: Ya Zhu (ya.zhu@cern.ch),CCNU & NBI, 2016       //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class TList;
class AliAnalysisTaskUniFlowMultiStrange;
AliAnalysisTaskUniFlowMultiStrange* AddTaskUniFlowMultiStrange(Bool_t IsGrid, TString name = "UniFlow")
{
  if (IsGrid){
  gGrid->Connect("alien://");
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;

  TString fileName = AliAnalysisManager::GetCommonFileName();   // by default, a file is open for writing. here, we get the filename
  fileName += Form(":%s",name.Data());      // create a subfolder in the file

  AliAnalysisTaskUniFlowMultiStrange* task = new AliAnalysisTaskUniFlowMultiStrange(name.Data()); // now we create an instance of your task
  if(!task) return 0x0;

  mgr->AddTask(task); // add your task to the manager

  // Creating containers
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* cOutput1 = mgr->CreateContainer(Form("Flow_Refs_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput2 = mgr->CreateContainer(Form("Flow_Charged_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput3 = mgr->CreateContainer(Form("Flow_Pion_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput4 = mgr->CreateContainer(Form("Flow_Kaon_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput5 = mgr->CreateContainer(Form("Flow_Proton_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput6 = mgr->CreateContainer(Form("Flow_K0s_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput7 = mgr->CreateContainer(Form("Flow_Lambda_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput8 = mgr->CreateContainer(Form("Flow_Phi_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data())); 
  AliAnalysisDataContainer* cOutput9 = mgr->CreateContainer(Form("Flow_Xi_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput10 = mgr->CreateContainer(Form("Flow_Omega_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput11 = mgr->CreateContainer(Form("QA_Events_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput12 = mgr->CreateContainer(Form("QA_Charged_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput13 = mgr->CreateContainer(Form("QA_PID_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput14 = mgr->CreateContainer(Form("QA_V0s_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput15 = mgr->CreateContainer(Form("QA_Phi_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput16 = mgr->CreateContainer(Form("Flow_Weights_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));

  // Connecting containers to task
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,cOutput1);
  mgr->ConnectOutput(task,2,cOutput2);
  mgr->ConnectOutput(task,3,cOutput3);
  mgr->ConnectOutput(task,4,cOutput4);
  mgr->ConnectOutput(task,5,cOutput5);
  mgr->ConnectOutput(task,6,cOutput6);
  mgr->ConnectOutput(task,7,cOutput7);
  mgr->ConnectOutput(task,8,cOutput8);
  mgr->ConnectOutput(task,9,cOutput9);
  mgr->ConnectOutput(task,10,cOutput10);
  mgr->ConnectOutput(task,11,cOutput11);
  mgr->ConnectOutput(task,12,cOutput12);
  mgr->ConnectOutput(task,13,cOutput13);
  mgr->ConnectOutput(task,14,cOutput14);
  mgr->ConnectOutput(task,15,cOutput15);
  mgr->ConnectOutput(task,16,cOutput16);

  return task;
}
