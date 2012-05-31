#ifndef __CINT__
#endif
void AddAnalysisTaskMixInfo(TString opts = "")
{
   // create manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) return;

   // create our task
   AliAnalysisTaskMixInfo *task = new AliAnalysisTaskMixInfo("AliAnalysisTaskMixInfo");
   Int_t debugLevel = 1;
   TString myclasses = "";
//    myclasses += ":AliAnalysisTaskMixInfo";
//    myclasses += ":AliAnalysisTaskEx02";

   if (!myclasses.IsNull()) task->SetLogType(AliLog::kDebug + debugLevel, myclasses);

   // create mix output container
   AliAnalysisDataContainer *outputMix = mgr->CreateContainer("cMixInfoList", TList::Class(), AliAnalysisManager::kOutputContainer, Form("MixInfo%s.root", opts.Data()));

   // add our task to the manager
   mgr->AddTask(task);

   // finaly connect input and output
   mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, outputMix);

}
