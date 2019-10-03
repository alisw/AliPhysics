/* AddTask macro for class AliAnalysisTaskJetMatching
 * Redmer Alexander Bertens, rbertens@cern.ch
 * Utrecht University, Utrecht, Netherlands             */

AliAnalysisTaskJetMatching* AddTaskJetMatching(
        const char* sourceJets  = "SourceJets", // source jets
        const char* targetJets  = "TargetJets", // target jets
        const char* matchedJets = "MatchedJets",// matched jets
        UInt_t matchingScheme   = AliAnalysisTaskJetMatching::kGeoR,
        Bool_t matchConstituents= kTRUE,
        Float_t minFrReCon      = .3,
        Float_t minFrReConPt    = .5,
        const char *name        = "AliAnalysisTaskJetMatching",
        Bool_t cut              = kTRUE,
        const char*  sourceType = "TPC",
        Float_t sourceRadius    = 0.3,
        Float_t sourceAreaCut   = .557,
        Float_t sourcePtBias    = 10.,
        const char* targetType  = "TPC",
        Float_t targetRadius    = 0.3,
        Float_t targetAreaCut   = .557,
        Float_t targetPtBias    = 10.
  )
{ 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)                             return 0x0;
  if (!mgr->GetInputEventHandler())     return 0x0;

  AliAnalysisTaskJetMatching* jetTask = new AliAnalysisTaskJetMatching(name);
  jetTask->SetMatchConstituents(matchConstituents);
  jetTask->SetMinFracRecoveredConstituents(minFrReCon);
  jetTask->SetMinFracRecoveredConstituentPt(minFrReConPt);
  jetTask->SetSourceJetsName(sourceJets);
  jetTask->SetTargetJetsName(targetJets);
  jetTask->SetMatchedJetsName(matchedJets);
  jetTask->SetMatchingScheme(matchingScheme);
  // if we want the jet package to cut on the source and target jets
  jetTask->SetUseEmcalBaseJetCuts(cut);
  if(cut) {
      AliJetContainer* sourceContainer = jetTask->AddJetContainer(sourceJets, sourceType, sourceRadius);
      if(sourceContainer) {
          sourceContainer->SetName("sourceJetContainer");
          sourceContainer->SetPercAreaCut(sourceAreaCut);
          sourceContainer->SetPtBiasJetTrack(sourcePtBias);
      }
      AliJetContainer* targetContainer = jetTask->AddJetContainer(targetJets, targetType, targetRadius);
      if(targetContainer) {
          targetContainer->SetName("targetJetContainer");
          targetContainer->SetPercAreaCut(targetAreaCut);
          targetContainer->SetPtBiasJetTrack(targetPtBias);
      }
  }
  
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
