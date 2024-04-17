AliAnalysisTaskEmcalPIDinJet* AddTaskEmcalPIDinJet(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char *suffix             = ""
)
{

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {  
    ::Error("AddTaskHFjetTagHFE", "No analysis manager to connect to.");
    return 0;
  }   

  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskHFjetTagHFE", "This task requires an input event handler");
    return 0;
  }

  AliAnalysisTaskEmcalPIDinJet *jetTask = 0;
  jetTask = AddTaskEmcalJet2023("usedefault", "usedefault", "usedefault", "");
  jetTask->GetClusterContainer(0)->SetClusECut(0.);
  jetTask->GetClusterContainer(0)->SetClusPtCut(0.);
  jetTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
  jetTask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
  jetTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  jetTask->GetParticleContainer(0)->SetParticlePtCut(0.15);;
  jetTask->SetHistoBins(300, 0, 300);
  jetTask->SelectCollisionCandidates(kPhysSel);

  AliJetContainer* jetCont02 = sampleTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.2, AliEmcalJet::kTPCfid, "Jet");
  AliJetContainer* jetCont03 = sampleTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.3, AliEmcalJet::kTPCfid, "Jet");
  AliJetContainer* jetCont04 = sampleTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.4, AliEmcalJet::kTPCfid, "Jet");  
 
  jetCont02->SetPercAreaCut(0.6);
  jetCont03->SetPercAreaCut(0.6);
  jetCont04->SetPercAreaCut(0.6);
  if (iBeamType != AliAnalysisTaskEmcal::kpp) {
       jetCont02->SetRhoName("Rho");
       jetCont03->SetRhoName("Rho");
       jetCont04->SetRhoName("Rho");
    }


  if(!jetTask) return 0x0;
  
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


