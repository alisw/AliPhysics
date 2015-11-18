// $Id$

AliAnalysisTaskEmcalJetSample* AddTaskEmcalJetSample(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Int_t       nCentBins          = 1,
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.6,
  const char *type               = "EMCAL",
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
  if (strcmp(type,"")) {
    name += "_";
    name += type;
  }

  Printf("name: %s",name.Data());

  AliAnalysisTaskEmcalJetSample* jetTask = new AliAnalysisTaskEmcalJetSample(name);
  jetTask->SetCentRange(0.,100.);
  jetTask->SetNCentBins(nCentBins);

  AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
  if(trackCont) trackCont->SetClassName("AliVTrack");
  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

  TString strType(type);
  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
  if(jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    //jetCont->SetZLeadingCut(0.98,0.98);
    jetCont->SetPercAreaCut(jetareacut);
    jetCont->SetJetPtCut(jetptcut);
    jetCont->SetLeadingHadronType(leadhadtype);
  }

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

AliAnalysisTaskEmcalJetSample* AddTaskEmcalJetSample( AliEmcalJetTask* jetFinderTask,
  Int_t       nCentBins          = 1,
  Double_t    jetareacut         = 0.6,
  const char *type               = "EMCAL",
  Int_t       leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskEmcalJetSample"
)
    {
    const char* ntracks            = jetFinderTask->GetTracksName();
    const char* nclusters          = jetFinderTask->GetClusName();
    const char* njets              = jetFinderTask->GetJetsName();
    const char* nrho               = jetFinderTask->GetRhoName();
    Double_t    jetradius          = jetFinderTask->GetRadius();
    Double_t    jetptcut           = jetFinderTask->GetMinJetPt();

    AliAnalysisTaskEmcalJetSample* jetTask = AddTaskEmcalJetSample(ntracks , nclusters, njets, nrho, nCentBins, jetradius, jetptcut, jetareacut, type, leadhadtype, taskname);

    return jetTask;
    }
