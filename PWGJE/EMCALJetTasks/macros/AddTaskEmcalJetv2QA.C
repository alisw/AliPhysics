// Jet v2 task using QA method, based on jet sample task (S.Aiola).
//
// Authors: Jason Mueller (CERN summer student 2014) & Alice Ohlson

AliAnalysisTaskEmcalJetv2QA* AddTaskEmcalJetv2QA(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Double_t    jetv2              = 0.0,
  Bool_t      ptweight           = kFALSE,
  Int_t       nCentBins          = 1,
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.6,
  const char *type               = "TPC",
  Int_t       leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskEmcalJetv2QA"
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetv2QA", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetv2QA", "This task requires an input event handler");
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
  name += Form("_v%.2i",(Int_t)(10*jetv2));
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (ptweight)
    name += "_ptweight";
  if (strcmp(type,"")) {
    name += "_";
    name += type;
  }

  Printf("name: %s",name.Data());

  AliAnalysisTaskEmcalJetv2QA* jetTask = new AliAnalysisTaskEmcalJetv2QA(name);

  jetTask->SetJetv2(jetv2); // set your jet v2 here
  jetTask->SetDoPtWeight(ptweight); // set doPtWeight

  jetTask->SetCentRange(0.,100.);
  jetTask->SetNCentBins(nCentBins);

  AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
  trackCont->SetClassName("AliVTrack");
  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

  TString strType(type);
  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
  if(jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    jetCont->SetZLeadingCut(0.98,0.98);
    jetCont->SetPercAreaCut(0.6);
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
