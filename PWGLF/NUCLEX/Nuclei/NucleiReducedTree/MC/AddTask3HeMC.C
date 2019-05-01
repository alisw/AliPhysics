AliAnalysisTask *AddTask3HeMC(Bool_t useMC, Bool_t isAOD=kTRUE)
{
   // specify appendix: track cuts or differences in the tasks
   TString appendix("Test");
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
   
   //create output containers
   TString containerName = mgr->GetCommonFileName();
   containerName += ":He3MC"; // create a subfolder in the file
   // Appendix to add the tack cuts to the folder name
   containerName += appendix.Data();
   printf("container name: %s\n", containerName.Data());
   
   AliAnalysisTask3HeMC *task = new AliAnalysisTask3HeMC(Form("He3MCTask%s",appendix.Data()));
   
   // Set data type
   if(isAOD) task->SetAODAnalysis();
   else task->SetESDAnalysis();
   
   // Trigger selection
   task->SelectCollisionCandidates(AliVEvent::kINT7);
   
   if (useMC) task->SetHasMCData(kTRUE);
   else task->SetHasMCData(kFALSE);
   
   // Pseudo-rapidity cut
   task->SetEtaCut(-0.9,0.9);
   // stronger DCA cut
   task->SetDCACut(1,2);

   mgr->AddTask(task);
   
   mgr->ConnectInput(task,  0, cinput);
   // create output containers for Results and QA
   mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("Results_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
   mgr->ConnectOutput(task,2, mgr->CreateContainer(Form("QA_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
   
   return NULL;
}
