// AddTaskDmesonJetCorr.C

AliAnalysisTaskDmesonJetCorrelations* AddTaskDmesonJetCorr(AliAnalysisTaskDmesonJetCorrelations::ECandidateType cand = AliAnalysisTaskDmesonJetCorrelations::kDstartoKpipi,
                                                           const char *cutfname           = "",
                                                           const char *ntracks            = "Tracks",
                                                           const char *nclusters          = "CaloClusters",
                                                           const char *njets              = "Jets",
                                                           const char *nrho               = "Rho",
                                                           Double_t    jetradius          = 0.2,
                                                           Double_t    jetptcut           = 1,
                                                           Double_t    jetareacut         = 0.557,
                                                           const char *cutType            = "TPC",
                                                           Int_t       leadhadtype        = 0,
							   const char *DmesonCandName     = "Dcandidates",
							   Bool_t      useExchCont        = kTRUE,
                                                           const char *taskname           = "AliAnalysisTaskDmesonJetCorrelations",
                                                           const char *suffix             = ""
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskDmesonJetCorr", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskDmesonJetCorr", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString cutsname;
  TString candname;
  switch (cand) {
  case AliAnalysisTaskDmesonJetCorrelations::kD0toKpi :
    candname = "D0";
    cutsname = "D0toKpiCuts";
    break;
  case AliAnalysisTaskDmesonJetCorrelations::kDstartoKpipi :
    candname = "DStar";
    cutsname = "DStartoKpipiCuts";
    break;
  default :
    ::Error("AddTaskDmesonJetCorr", "Candidate type %d not recognized!",cand);
    return NULL;
  }

  TString name(Form("%s_%s", taskname, candname.Data()));
  if (strcmp(suffix,"")) {
    name += "_";
    name += suffix;
  }
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  name += "_";
  name += cutType;

  AliRDHFCuts* analysiscuts = 0;
  TFile* filecuts = 0;
  if (strcmp(cutfname, "")) {
    filecuts = TFile::Open(cutfname);
    if (!filecuts || filecuts->IsZombie()) {
      ::Warning("AddTaskDmesonJetCorr", "Input file not found: use std cuts.");
      filecuts = 0;
    }
  }
  else {
    ::Info("AddTaskDmesonJetCorr", "No input file provided: use std cuts.");
  }

  if (filecuts) {
    analysiscuts = dynamic_cast<AliRDHFCuts*>(filecuts->Get(cutsname));
    if (!analysiscuts) {
      ::Warning("AddTaskDmesonJetCorr", "Could not find analysis cuts '%s' in '%s'. Using std cuts.", cutsname.Data(), cutfname);
    }
  }

  if (!analysiscuts) {
    switch (cand) {
    case AliAnalysisTaskDmesonJetCorrelations::kD0toKpi :
      analysiscuts = new AliRDHFCutsD0toKpi();
      break;
      
    case AliAnalysisTaskDmesonJetCorrelations::kDstartoKpipi :
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      break;
      
    default :
      break;
    }

    analysiscuts->SetStandardCutsPP2010();
    analysiscuts->SetName(cutsname);
  }

  AliAnalysisTaskDmesonJetCorrelations* jetTask = new AliAnalysisTaskDmesonJetCorrelations(name, analysiscuts, cand, useExchCont);
  jetTask->SetVzRange(-10,10);
  jetTask->SetNeedEmcalGeom(kFALSE);

  AliParticleContainer* trackCont = jetTask->AddParticleContainer(ntracks);
  AliClusterContainer* clusterCont = jetTask->AddClusterContainer(nclusters);

  AliJetContainer* jetCont = jetTask->AddJetContainer(njets, cutType, jetradius);
  if (jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->SetPercAreaCut(jetareacut);
    jetCont->SetJetPtCut(jetptcut);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    jetCont->SetLeadingHadronType(leadhadtype);
    jetCont->SetMaxTrackPt(1000);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer* cinput1 = mgr->GetCommonInputContainer();
  TString contname1(name);
  contname1 += "_histos";
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(contname1.Data(), 
							    TList::Class(), AliAnalysisManager::kOutputContainer,
							    Form("%s:SA_DmesonJetCorr", AliAnalysisManager::GetCommonFileName()));

  TString contname2(name);
  contname2 += "_cuts";
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contname2.Data(),
                                                            AliRDHFCuts::Class(), AliAnalysisManager::kOutputContainer,
                                                            Form("%s:SA_DmesonJetCorr", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput(jetTask, 0, cinput1);

  if (useExchCont) {
    TString nameContainerFC2(DmesonCandName);
    if (!nameContainerFC2.IsNull()) {
      nameContainerFC2 += candname;
      nameContainerFC2 += suffix;
  
      TObjArray* cnt = mgr->GetContainers();
      AliAnalysisDataContainer* coutputFC2 = static_cast<AliAnalysisDataContainer*>(cnt->FindObject(nameContainerFC2));
      if (!coutputFC2) {
	::Error("AddTaskDmesonJetCorr", "Could not find input container '%s'!", nameContainerFC2.Data());
      }
      mgr->ConnectInput(jetTask, 1, coutputFC2);
    }
  }

  mgr->ConnectOutput(jetTask, 1, coutput1);
  mgr->ConnectOutput(jetTask, 2, coutput2);
  
  return jetTask;
}
