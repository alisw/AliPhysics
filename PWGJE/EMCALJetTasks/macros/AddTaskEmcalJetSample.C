// $Id$

AliAnalysisTaskEmcalJetSample* AddTaskEmcalJetSample(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  UInt_t      type               = AliAnalysisTaskEmcal::kEMCAL,
  Int_t       leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskEmcalJetSample"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(taskname);
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (type == AliAnalysisTaskEmcal::kTPC) 
    name += "_TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "_EMCAL";
  else if (type == AliAnalysisTaskEmcal::kUser) 
    name += "_USER";

  AliAnalysisTaskEmcalJetSample* jetTask = new AliAnalysisTaskEmcalJetSample(name);
  jetTask->SetAnaType(type);
  jetTask->SetTracksName(ntracks);
  jetTask->SetClusName(nclusters);
  jetTask->SetJetsName(njets);
  jetTask->SetRhoName(nrho);
  jetTask->SetJetRadius(jetradius);
  jetTask->SetJetPtCut(jetptcut);
  jetTask->SetPercAreaCut(jetareacut);
  jetTask->SetLeadingHadronType(leadhadtype);
  
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
