// $Id$

AliJetResponseMaker* AddTaskJetResponseMaker(
  const char *ntracks1           = "Tracks",
  const char *nclusters1         = "CaloClusters",
  const char *njets1             = "Jets",
  const char *nrho1              = "Rho",
  Double_t    jetradius1         = 0.2,
  const char *ntracks2           = "MCParticles",
  const char *nclusters2         = "",
  const char *njets2             = "MCJets",
  const char *nrho2              = "",
  Double_t    jetradius2         = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  Double_t    jetBiasTrack       = 5,
  Double_t    jetBiasClus        = 1000,
  UInt_t      matching           = AliJetResponseMaker::kGeometrical,
  Double_t    maxDistance1       = 0.25,
  Double_t    maxDistance2       = 0.25,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
  Int_t       ptHardBin          = -999,
  Double_t    minCent            = -999,
  Double_t    maxCent            = -999,
  const char *taskname           = "AliJetResponseMaker",
  Bool_t      biggerMatrix       = kFALSE,
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

  TString name(Form("%s_%s_%s_Track%d_Clus%d_",taskname,njets1,njets2,(Int_t)floor(jetBiasTrack),(Int_t)floor(jetBiasClus)));

  if (type == AliAnalysisTaskEmcal::kTPC)
    name += "TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "EMCAL";
  else if (type == AliAnalysisTaskEmcal::kUser) 
    name += "USER";

  if (minCent != -999 && maxCent != -999) 
    name += Form("_Cent%d_%d", (Int_t)floor(minCent), (Int_t)floor(maxCent));

  if (ptHardBin != -999) 
    name += Form("_PtHard%d", ptHardBin);

  AliJetResponseMaker* jetTask = address;
  if (jetTask)
    new (jetTask) AliJetResponseMaker(name);
  else
    jetTask = new AliJetResponseMaker(name);

  jetTask->SetAnaType(type);
  jetTask->SetTracksName(ntracks1);
  jetTask->SetClusName(nclusters1);
  jetTask->SetJetsName(njets1);
  jetTask->SetRhoName(nrho1);
  jetTask->SetJetRadius(jetradius1);
  jetTask->SetTracks2Name(ntracks2);
  jetTask->SetClus2Name(nclusters2);
  jetTask->SetJets2Name(njets2);
  jetTask->SetRho2Name(nrho2);
  jetTask->SetJet2Radius(jetradius2);
  jetTask->SetJetPtCut(jetptcut);
  jetTask->SetPercAreaCut(jetareacut);
  jetTask->SetPtBiasJetTrack(jetBiasTrack);
  jetTask->SetPtBiasJetClus(jetBiasClus);
  jetTask->SetMatching(matching, maxDistance1, maxDistance2);
  jetTask->SetVzRange(-10,10);
  jetTask->SetPtHardBin(ptHardBin);
  jetTask->SetCentRange(minCent,maxCent);

  if (biggerMatrix) 
    jetTask->SetHistoBins(1000,0,500);
  
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
