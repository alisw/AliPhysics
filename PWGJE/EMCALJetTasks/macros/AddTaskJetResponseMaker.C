// $Id$

AliJetResponseMaker* AddTaskJetResponseMaker(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  const char *nmctracks          = "MCParticles",
  const char *nmcjets            = "MCJets",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  Double_t    jetBiasTrack       = 5,
  Double_t    jetBiasClus        = 5,
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

  TString name(Form("%s_%s_Track%d_Clus%d_R0%d_",taskname,njets,(Int_t)floor(jetBiasTrack),(Int_t)floor(jetBiasClus),(Int_t)floor(jetradius*100+0.5)));
  if (type == AliAnalysisTaskEmcal::kTPC)
    name += "TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "EMCAL";
  else if (type == AliAnalysisTaskEmcal::kUser) 
    name += "USER";
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
  jetTask->SetRhoName(nrho);
  jetTask->SetMCJetsName(nmcjets);
  jetTask->SetMCTracksName(nmctracks);
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
