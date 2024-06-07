AliAnalysisTaskEmcalPIDinJet* AddTaskEmcalPIDinJet(
  TString name = "name",
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  Double_t    radius             = 0.2,
  Bool_t iBeamType_PbPb          = kTRUE,          
  const char *suffix             = ""
)
{

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {  
    return 0;
  }   

  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    return 0;
  }

  AliAnalysisTaskEmcalPIDinJet *jetTask = new AliAnalysisTaskEmcalPIDinJet(name.Data());
  //jetTask = AddTaskEmcalJetPIDinJet("usedefault", "usedefault", "usedefault", "");
  jetTask->GetClusterContainer(0)->SetClusECut(0.);
  jetTask->GetClusterContainer(0)->SetClusPtCut(0.);
  jetTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
  jetTask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
  jetTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  jetTask->GetParticleContainer(0)->SetParticlePtCut(0.15);;
  jetTask->SetHistoBins(300, 0, 300);
  jetTask->SelectCollisionCandidates(AliVEvent::kCentral);

  AliJetContainer* jetCont = jetTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, radius, AliEmcalJet::kTPCfid, "Jet");
  if(jetCont)
    {
     jetCont->SetPercAreaCut(0.6);
     if (iBeamType_PbPb) {
        jetCont->SetRhoName("Rho");
     }
   }

  if(!jetTask) return 0x0;
  
  mgr->AddTask(jetTask);
    
  // Create containers for input/output
 /*
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                            TList::Class(),AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  */

    TString containerName = mgr->GetCommonFileName();
    containerName += ":PWGJE_PIDinJet";
    TString SubcontainerName = Form("PIDinJet");
    SubcontainerName += name;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
    mgr->ConnectInput(jetTask, 0, cinput);
    mgr->ConnectOutput(jetTask, 1, coutput1); 
  

  return jetTask;
} 


