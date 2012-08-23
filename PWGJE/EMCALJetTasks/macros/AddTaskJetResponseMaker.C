// $Id$

AliJetResponseMaker* AddTaskJetResponseMaker(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nmcjets            = "MCJets",
  const char *nmctracks          = "MCParticles",
  Double_t    jetradius          = 0.4,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.8,
  Double_t    ptcut              = 0.15,
  Double_t    jetBiasTrack       = 10,
  Double_t    jetBiasClus        = 10,
  Double_t    maxDistance        = 0.25,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
  Int_t       ptHardBin          = -999,
  const char *taskname           = "AliJetResponseMaker",
  AliJetResponseMaker* address   = 0
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetResponseMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetResponseMaker", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(taskname);
  name += "_";
  name += njets;
  name += "_Track";
  name += jetBiasTrack;
  name += "_Clus";
  name += jetBiasClus;
  name += "_R0";
  name += floor(jetradius*10+0.5);
  name += "_PtCut";
  name += floor(ptcut*1000+0.5);
  name += "_";
  if (type == AliAnalysisTaskEmcal::kTPC) 
    name += "TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "EMCAL";
  else if (type == AliAnalysisTaskEmcal::kTPCSmall) 
    name += "TPCSmall";
  else if (type == AliAnalysisTaskEmcal::kEMCALOnly) 
    name += "EMCALOnly";
  if (ptHardBin != -999) {
    name += "_PtHard";
    name += ptHardBin;
  }

  AliJetResponseMaker* jetTask = address;
  if (jetTask)
    new (jetTask) AliJetResponseMaker(name);
  else
    jetTask = new AliJetResponseMaker(name);

  jetTask->SetAnaType(type);
  jetTask->SetTracksName(ntracks);
  jetTask->SetClusName(nclusters);
  jetTask->SetJetsName(njets);
  jetTask->SetMCJetsName(nmcjets);
  jetTask->SetMCTracksName(nmctracks);
  jetTask->SetPtCut(ptcut);
  jetTask->SetJetRadius(jetradius);
  jetTask->SetJetPtCut(jetptcut);
  jetTask->SetPercAreaCut(jetareacut);
  jetTask->SetPtBiasJetTrack(jetBiasTrack);
  jetTask->SetPtBiasJetClus(jetBiasClus);
  jetTask->SetMaxDistance(maxDistance);
  jetTask->SetVzRange(-10,10);
  jetTask->SetPtHardBin(ptHardBin);
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
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
