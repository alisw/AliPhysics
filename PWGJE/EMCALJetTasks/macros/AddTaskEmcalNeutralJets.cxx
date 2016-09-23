AliAnalysisTaskEmcalNeutralJets * AddTaskEmcalNeutralJets(
    const char *clusters = "usedefault",
    const char *suffix = ""
){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname = "NeutralJetTask", suffixstring(suffix);
  if(suffixstring.Length()) taskname += suffixstring;

  AliAnalysisTaskEmcalNeutralJets *jettask = new AliAnalysisTaskEmcalNeutralJets(taskname.Data());
  mgr->AddTask(jettask);

  TString fNClusters = clusters;
  if(fNClusters.Contains("usedefault")){
    AliInputEventHandler *inputhandler = mgr->GetInputEventHandler();
    if(inputhandler->IsA() == AliESDInputHandler::Class()){
      // ESDs
      fNClusters = "CaloClusters";
    } else {
      // AODs
      fNClusters = "caloClusters";
    }
  }

  AliClusterContainer *clustcont = jettask->AddClusterContainer(fNClusters.Data());

  AliJetContainer *r02cont = jettask->AddJetContainer(
      AliJetContainer::kNeutralJet,
      AliJetContainer::antikt_algorithm,
      AliJetContainer::pt_scheme,
      0.2,
      AliJetContainer::kEMCALfid,
      0x0,
      clustcont,
      "Jet"
      );
  r02cont->SetMinPt(10.);

  AliJetContainer *r04cont = jettask->AddJetContainer(
      AliJetContainer::kNeutralJet,
      AliJetContainer::antikt_algorithm,
      AliJetContainer::pt_scheme,
      0.2,
      AliJetContainer::kEMCALfid,
      0x0,
      clustcont,
      "Jet"
      );
  r04cont->SetMinPt(10.);

  jettask->SetContainers(r02cont->GetName(), r04cont->GetName(), clustcont->GetName());

  return jettask;
}
