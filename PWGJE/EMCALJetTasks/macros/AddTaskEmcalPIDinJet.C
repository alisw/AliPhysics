AliAnalysisTaskEmcalPIDinJet* AddTaskEmcalPIDinJet(
  TString name = "name",
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  Bool_t iBeamType_PbPb          = kTRUE,          
  const char *suffix             = ""
)
{
   cout << "loading PIDinJET task" << endl;
 
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
  cout << "Setting PIDinJET task parameters" << endl;
  jetTask->SetHistoBins(300, 0, 300);
  jetTask->SelectCollisionCandidates(AliVEvent::kCentral);
 
  AliTrackContainer* trackCont = 0;
  trackCont = jetTask->AddTrackContainer("tracks");
  trackCont->SetFilterHybridTracks(kTRUE);

   cout << "++++++++++++++++ Set PIDinJET task parameters" << endl;

  AliJetContainer* jetCont02 = jetTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.2, AliEmcalJet::kTPCfid, "Jet");
  AliJetContainer* jetCont03 = jetTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.3, AliEmcalJet::kTPCfid, "Jet");
  AliJetContainer* jetCont04 = jetTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.4, AliEmcalJet::kTPCfid, "Jet");  
     
     jetCont02->SetPercAreaCut(0.6);
     jetCont03->SetPercAreaCut(0.6);
     jetCont04->SetPercAreaCut(0.6);
     if (iBeamType_PbPb) {
        jetCont02->SetRhoName("Rho");
        jetCont03->SetRhoName("Rho");
        jetCont04->SetRhoName("Rho");
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
    //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectInput(jetTask, 0, cinput);
    mgr->ConnectOutput(jetTask, 1, coutput1); 
  

  return jetTask;
} 


