// $Id$

AliAnalysisTaskDcalDijetPerf* AddTaskDcalDijetPerf(
						   const char *ntracks            = "Tracks",
						   const char *nclusters          = "CaloClusters",
						   const char *njets              = "Jets",
						   const char *njets2             = "Jets2",
                           const char *njets3             = "Jets3",
						   const char *nrho               = "Rho",
						   Int_t       nCentBins          = 1,
						   Double_t    jetradius          = 0.2,
						   Double_t    jetradius2         = 0.2,
                           Double_t    jetradius3         = 0.2,
						   Double_t    jetptcut           = 1,
						   Double_t    jetareacut         = 0.6,
						   const char *type               = "TPC",
						   Int_t       leadhadtype        = 0,
						   const char *taskname           = "AliAnalysisTaskDcalDijetPerf"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskDcalDijetPerf", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskDcalDijetPerf", "This task requires an input event handler");
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
  if (strcmp(njets2,"")) {
    name += "_";
    name += njets2;
  }
  if (strcmp(njets3,"")) {
    name += "_";
    name += njets3;
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

  AliAnalysisTaskDcalDijetPerf* jetTask = new AliAnalysisTaskDcalDijetPerf(name);
  jetTask->SetCentRange(0.,100.);
  jetTask->SetNCentBins(nCentBins);

  AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
  trackCont->SetClassName("AliVTrack");
  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

  TString strType(type);
  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
  AliJetContainer *jetCont2 = jetTask->AddJetContainer(njets2,strType,jetradius2);
  AliJetContainer *jetCont3 = jetTask->AddJetContainer(njets3,strType,jetradius3);
  if(jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    //jetCont->SetZLeadingCut(0.98,0.98);
    //jetCont->SetPercAreaCut(0.6);
    jetCont->SetJetPtCut(jetptcut);    
    jetCont->SetLeadingHadronType(leadhadtype);
  }
  if(jetCont2) {
    jetCont2->SetRhoName(nrho);
    jetCont2->ConnectParticleContainer(trackCont);
    jetCont2->ConnectClusterContainer(clusterCont);
    //jetCont->SetZLeadingCut(0.98,0.98);
    //jetCont->SetPercAreaCut(0.6);
    jetCont2->SetJetPtCut(jetptcut);
    jetCont2->SetLeadingHadronType(leadhadtype);
  }
    
  if(jetCont3) {
        jetCont3->SetRhoName(nrho);
        jetCont3->ConnectParticleContainer(trackCont);
        jetCont3->ConnectClusterContainer(clusterCont);
        //jetCont->SetZLeadingCut(0.98,0.98);
        //jetCont->SetPercAreaCut(0.6);
        jetCont3->SetJetPtCut(jetptcut);
        jetCont3->SetLeadingHadronType(leadhadtype);
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

AliAnalysisTaskDcalDijetPerf* AddTaskDcalDijetPerf( AliEmcalJetTask* jetFinderTask,
						    Int_t       nCentBins          = 1,
						    Double_t    jetareacut         = 0.6,
						    const char *type               = "EMCAL",
						    Int_t       leadhadtype        = 0,
						    const char *taskname           = "AliAnalysisTaskDcalDijetPerf"
						    )
{
  const char* ntracks            = jetFinderTask->GetTracksName();
  const char* nclusters          = jetFinderTask->GetClusName();
  const char* njets              = jetFinderTask->GetJetsName();
  const char* nrho               = jetFinderTask->GetRhoName();
  Double_t    jetradius          = jetFinderTask->GetRadius();
  Double_t    jetptcut           = jetFinderTask->GetMinJetPt();
  
  AliAnalysisTaskDcalDijetPerf* jetTask = AddTaskDcalDijetPerf(ntracks , nclusters, njets, nrho, nCentBins, jetradius, jetptcut, jetareacut, type, leadhadtype, taskname);
  
  return jetTask;
}
