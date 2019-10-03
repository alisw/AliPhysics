AliAnalysisTask* AddTaskCompactTreeMaker(const char* outputFileName="", const char* ocdbPath="raw://")
{
   AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr)
   {
       ::Error("AddTaskCompactTreeMaker","No analysis manager to connect to...");
       return 0x0;
   }

   if (!mgr->GetInputEventHandler())
   {
       ::Error("AddTaskCompactTreeMaker","This task requires an input event handler !");
       return 0x0;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD"))
   {
       ::Error("AddTaskCompactTreeMaker","This task needs an ESD input handler !");
       return 0x0;
   }
   
   AliMCEventHandler* mcH = static_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
   if (!mcH)
   {
       ::Error("AddTaskCompactTreeMaker","This task needs a MC event handler to be connected first !");
       return 0x0;
   }

   AliAnalysisTask* task = new AliMuonCompactTreeMaker(ocdbPath);

   mgr->AddTask(task);

   TString outputFile = mgr->GetCommonFileName();
   outputFile += ":";
   outputFile += "MUON_COMPACT";

   if (strlen(outputFileName)) outputFile = outputFileName;

   AliAnalysisDataContainer* output = mgr->CreateContainer("compactevents",TTree::Class(),AliAnalysisManager::kOutputContainer,outputFile.Data());


   mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task,1,output);

   return task;
}
