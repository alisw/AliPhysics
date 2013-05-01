AliAnalysisTaskHJetEmbed *AddTaskHJetEmbedPYTHIA(const TString period      = "lhc11h",
						 const Bool_t runQA        = kTRUE, 
						 const Bool_t runHJet      = kTRUE,
						 const Bool_t runMatch     = kTRUE,
						 const TString trackName   = "PicoTracks",
						 const TString clusterName = "",
						 const TString mcName      = "MCParticlesSelected",
						 const Double_t jetRadius  = 0.4,
						 const TString jetArrayPL  = "",
						 const TString jetArrayDL  = "",
						 const TString jetArray    = "",
						 const TString rhoName     = "",
						 const Double_t minTrkPt   = 0.15, const Double_t maxTrkPt = 1e4,
						 const Double_t minTrkEta  = -0.9, const Double_t maxTrkEta = 0.9,
						 const Double_t minTrkPhi  = 0,    const Double_t maxTrkPhi = 2*TMath::Pi()
						 )
{

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHJetEmbedPYTHIA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskHJetEmbedPYTHIA", "This task requires an input event handler");
    return NULL;
  }
  
  //Analyze
  AliAnalysisTaskHJetEmbed *taskEmbed = new AliAnalysisTaskHJetEmbed("AnaEmbedPYTHIA");
  taskEmbed->SetRunPeriod(period.Data());
  taskEmbed->SetTrkPtRange(minTrkPt, maxTrkPt);
  taskEmbed->SetTrkPhiRange(minTrkPhi, maxTrkPhi);
  taskEmbed->SetTrkEtaRange(minTrkEta, maxTrkEta);
  taskEmbed->SetRadius(jetRadius);
  taskEmbed->SetTrackArrName(trackName.Data());
  taskEmbed->SetMCParticleArrName(mcName.Data());
  taskEmbed->SetRhoName(rhoName.Data());
  taskEmbed->SetJetArrName(jetArray.Data());
  taskEmbed->SetPLJetArrName(jetArrayPL.Data());
  taskEmbed->SetDLJetArrName(jetArrayDL.Data());
  taskEmbed->SetRunQA(runQA);
  taskEmbed->SetRunHJet(runHJet);
  taskEmbed->SetRunMatch(runMatch);
  if(period.Contains("lhc10h",TString::kIgnoreCase) || period.Contains("lhc11h",TString::kIgnoreCase))
    taskEmbed->SetCollisionSystem("PbPb");
  else
    taskEmbed->SetCollisionSystem("pp");

  if(period.Contains("lhc11h",TString::kIgnoreCase))
    taskEmbed->SelectCollisionCandidates(AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral);

  AliAnalysisDataContainer *coutput = mgr->CreateContainer("EmbedPYTHIA",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s",  AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput(taskEmbed,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEmbed,1,coutput);

  return taskEmbed;
}
