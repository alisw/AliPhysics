// AddTaskChargedJetsHadronToy.C

AliAnalysisTaskChargedJetsHadronToy* AddTaskChargedJetsHadronToy(
  const char *trackInput         = "tracks",
  const char *trackOutput        = "tracks_toy",
  const char *distributionFile   = "alien:///alice/cern.ch/user/r/rhaake/Distributions_PbPb.root"
)
{
  cout << " ############ MACRO EXECUTION STARTED: AddTaskChargedJetsHadronToy.C ############\n";
  //==============================================================================
  // Prepare analysis manager, containers, etc.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    ::Error("AddTaskChargedJetsHadronToy", "No analysis manager to connect to.");
    return NULL;
  }  
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskChargedJetsHadronToy", "This task requires an input event handler");
    return NULL;
  }
  
  //==============================================================================
  // Adding and configuring tasks

  AliAnalysisDataContainer* contHistos = mgr->CreateContainer("Toy_histos", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));

  AliAnalysisTaskChargedJetsHadronToy* toyModel = new AliAnalysisTaskChargedJetsHadronToy();
  toyModel->SetInputArrayName(trackInput);
  toyModel->SetOutputArrayName(trackOutput);

  if(distributionFile)
  {
    TFile* fileInput = TFile::Open(distributionFile);
    if(!fileInput) std::cout << "Distributions file not found!" << std::endl;
    TH1* distPt = static_cast<TH1*>(fileInput->Get("Pt"));
    TH1* distMult = static_cast<TH1*>(fileInput->Get("Multiplicity"));
    TH1* distPhiEta = static_cast<TH1*>(fileInput->Get("PhiEta"));

    TH2* distV2 = static_cast<TH2*>(fileInput->Get("v2_EP_PbPb"));
    TH2* distV3 = static_cast<TH2*>(fileInput->Get("v3_EP_PbPb"));
    TH2* distV4 = static_cast<TH2*>(fileInput->Get("v4_EP_PbPb"));

    toyModel->SetDistributionMultiplicity(distMult);
    toyModel->SetDistributionPt(distPt);
    toyModel->SetDistributionEtaPhi(distPhiEta);
    toyModel->SetDistributionV2(distV2);
    toyModel->SetDistributionV3(distV3);
    toyModel->SetDistributionV4(distV4);
  }
  mgr->AddTask(toyModel);

  //==============================================================================
  // Finalization

  mgr->ConnectInput  (toyModel, 0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (toyModel, 1, contHistos );

  cout << " ############ MACRO EXECUTION DONE: AddTaskChargedJetsHadronToy.C ############\n";
 
  return toyModel;
}
