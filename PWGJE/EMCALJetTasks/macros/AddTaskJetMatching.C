/* AddTask macro for class AddTaskJetMatching.C
 * Redmer Alexander Bertens, rbertens@cern.ch
 * Utrecht University, Utrecht, Netherlands             */

AliAnalysisTaskJetMatching* AddTaskJetMatching(
        const char* sourceJets  = "SourceJets", // source jets
        const char* targetJets  = "TargetJets", // target jets
        const char* matchedJets = "MatchedJets",// matched jets
        UInt_t matchingScheme   = AliAnalysisTaskJetMatching::kDeepMatching,
        UInt_t duplicateRecovery= AliAnalysisTaskJetMatching::kDoNothing,
        const char *name        = "AliAnalysisTaskJetMatching",
  )
{ 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)                             return 0x0;
  if (!mgr->GetInputEventHandler())     return 0x0;

  AliAnalysisTaskJetMatching* jetTask = new AliAnalysisTaskJetMatching(name);
  jetTask->SetDebugMode(-1);
  jetTask->SetSourceJetsName(sourceJets);
  jetTask->SetTargetJetsName(targetJets);
  jetTask->SetMatchedJetsName(matchedJets);
  jetTask->SetMatchingScheme(matchingScheme);
  jetTask->SetDuplicateRecoveryScheme(duplicateRecovery);
  mgr->AddTask(jetTask);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  return jetTask;
}
