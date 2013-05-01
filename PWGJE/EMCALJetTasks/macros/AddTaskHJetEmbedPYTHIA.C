AliAnalysisTaskHJetEmbed *AddTaskHJetEmbedPYTHIA(const TString period      = "lhc11h",
						 const Bool_t runQA        = kTRUE, 
						 const Bool_t runHJet      = kTRUE,
						 const Bool_t runMatch     = kTRUE,
						 const char *trackName     = "PicoTracks",
						 const char *clusterName   = "",
						 const char *mcName        = "MCParticlesSelected",
						 const UInt_t jetType      = AliEmcalJetTask::kAKT | AliEmcalJetTask::kChargedJet | AliEmcalJetTask::kR040Jet,
						 const UInt_t jetTypeBkg   = AliEmcalJetTask::kKT | AliEmcalJetTask::kChargedJet | AliEmcalJetTask::kR040Jet,
						 const Double_t jetRadius  = 0.4,
						 const Int_t pSel           = AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral, 
						 const Double_t minTrkPt   = 0.15, const Double_t maxTrkPt = 1e4,
						 const Double_t minTrkEta  = -0.9, const Double_t maxTrkEta = 0.9,
						 const Double_t minTrkPhi  = 0,    const Double_t maxTrkPhi = 2*TMath::Pi(),
						 const Bool_t   hybridRho  = kTRUE,
						 TString outfileName       = "")
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

  // Find particle level jets
  AliEmcalJetTask* jetPLFinderTask =  AddTaskEmcalJet(jetType,mcName,clusterName,minTrkPt,0.30,0.005,jetRadius,"PLJet");
  jetPLFinderTask->SetEtaRange(minTrkEta,maxTrkEta);
  jetPLFinderTask->SetPhiRange(minTrkPhi,maxTrkPhi);
  jetPLFinderTask->SelectCollisionCandidates(pSel);

  // Find detector level jets
  AliEmcalJetTask* jetDLFinderTask =  AddTaskEmcalJet(jetType,trackName,clusterName,minTrkPt,0.30,0.005,jetRadius,"DLJet");
  jetDLFinderTask->SetEtaRange(minTrkEta,maxTrkEta);
  jetDLFinderTask->SetPhiRange(minTrkPhi,maxTrkPhi);
  jetDLFinderTask->SelectConstituents(TObject::kBitMask, 0); // select tracks with label!=0TObject::kBitMap, 0
  jetDLFinderTask->SelectCollisionCandidates(pSel);

  // Calculate rho
  AliEmcalJetTask* jetFinderTaskBkg = AddTaskEmcalJet(jetType, trackName, clusterName,minTrkPt,0.30,0.005,jetRadius,"BkgJet");
  jetFinderTaskBkg->SetEtaRange(minTrkEta,maxTrkEta);
  jetFinderTaskBkg->SetPhiRange(minTrkPhi,maxTrkPhi);
  if(hybridRho) jetFinderTaskBkg->SelectConstituents(0,0);
  else          jetFinderTaskBkg->SelectConstituents(0,TObject::kBitMask);
  jetFinderTaskBkg->SelectCollisionCandidates(pSel);

  if(hybridRho) char *rhoName = Form("Hybrid_Rho");
  else          char *rhoName = Form("Rho");
  AliAnalysisTaskRho *rhochtask = AddTaskRho(jetFinderTaskBkg->GetName(),trackName, clusterName, rhoName, jetRadius, AliAnalysisTaskEmcal::kTPC, 0.001, 0, 0, 2, kFALSE, rhoName);
  rhochtask->SetJetEtaLimits(minTrkEta+jetRadius,maxTrkEta-jetRadius);
  rhochtask->SetJetPhiLimits(minTrkPhi,maxTrkPhi);
  rhochtask->SelectCollisionCandidates(pSel);

  // Find embedded jets
  AliEmcalJetTask* jetFinderTask =  AddTaskEmcalJet(jetType,trackName,clusterName,minTrkPt,0.30,0.005,jetRadius,"Jet");
  jetFinderTask->SetEtaRange(minTrkEta,maxTrkEta);
  jetFinderTask->SetPhiRange(minTrkPhi,maxTrkPhi);
  jetFinderTask->SelectConstituents(0,0); 
  jetFinderTask->SelectCollisionCandidates(pSel);

  
  //Analyze
  AliAnalysisTaskHJetEmbed *taskEmbed = new AliAnalysisTaskHJetEmbed("AnaEmbedPYTHIA");
  taskEmbed->SetRunPeriod(period.Data());
  taskEmbed->SetTrackArrName(trackName);
  taskEmbed->SetRhoName(rhoName);
  taskEmbed->SetJetArrName(jetFinderTask->GetName());
  taskEmbed->SetPLJetArrName(jetPLFinderTask->GetName());
  taskEmbed->SetDLJetArrName(jetDLFinderTask->GetName());
  taskEmbed->SetMCParticleArrName(mcName);
  taskEmbed->SetRunQA(runQA);
  taskEmbed->SetRunHJet(runHJet);
  taskEmbed->SetRunMatch(runMatch);
  if(period.Contains("lhc10h",TString::kIgnoreCase) || period.Contains("lhc11h",TString::kIgnoreCase))
    taskEmbed->SetCollisionSystem("PbPb");
  else
    taskEmbed->SetCollisionSystem("pp");
  taskEmbed->SelectCollisionCandidates(pSel);

  mgr->AddTask(taskEmbed);
  if(outfileName.IsNull()) outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("EmbedPYTHIA",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s", outfileName.Data()));
  mgr->ConnectInput(taskEmbed,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEmbed,1,coutput);

  return taskEmbed;
}
