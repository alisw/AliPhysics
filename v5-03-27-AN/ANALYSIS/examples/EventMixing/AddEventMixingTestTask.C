#ifndef __CINT__
#endif
AliAnalysisTask *AddEventMixingTestTask(TString format = "esd", Bool_t useMC = kFALSE,TString postfix="")
{
   // create manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) mgr = new AliAnalysisManager("MIX test");

   // create our task
   AliAnalysisTaskEx02 *task = new AliAnalysisTaskEx02("AliAnalysisTaskEx02");

   // create output container
   AliAnalysisDataContainer *output1 = mgr->CreateContainer("cEx2", TList::Class(), AliAnalysisManager::kOutputContainer, "MixTestOutput.root");

   // add our task to the manager
   mgr->AddTask(task);

   // finaly connect input and output
   mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output1);


   return task;
}

