/*
 * AddTask macro for class 
 * Redmer Alexander Bertens, rbertens@cern.ch
 * Utrecht University, Utrecht, Netherlands
 *
 * Note: this macro is pretty much a copy of AddTaskEmcalJetSample.C
 *
 */

AliAnalysisTaskLocalRho* AddTaskLocalRho(
  const char *ntracks           = "Tracks",     // track selection used for vn estimate
  const char *nclusters         = "",           // clusters (not used) 
  const char *njets             = "Jets",       // jet selection for finding leading jet
  const char *nrho              = "Rho",        // name of nominal rho
  const char *lrho              = "LocalRho",   // name of local rho object
  Double_t   jetradius          = 0.3,          // jet radius (to remove leading jet)
  Double_t   jetptcut           = 1,            
  Double_t   jetareacut         = 0.557,
  UInt_t     type               = AliEmcalJet::kTPC,
  Int_t      leadhadtype        = 0,
  const char *name              = "AliAnalysisTaskLocalRho",    // task name
  TString    fitOpts            = "WLQI",                       // options for tfitter
  UInt_t     fitType            = AliAnalysisTaskLocalRho::kCombined,   // fitting strategy
  TArrayI    *centralities      = 0x0,                                  // centrality binning for qa
  UInt_t     runMode            = AliAnalysisTaskLocalRho::kGrid        // run mode 
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

  AliAnalysisTaskLocalRho* jetTask = new AliAnalysisTaskLocalRho(name, runMode);
  // inherited setters
  jetTask->AddJetContainer(njets);
  jetTask->SetAnaType(type);
  jetTask->AddParticleContainer(ntracks);
  jetTask->AddClusterContainer(nclusters);
  jetTask->SetRhoName(nrho);
  jetTask->SetLocalRhoName(lrho);
  jetTask->SetJetRadius(jetradius);
  jetTask->SetJetPtCut(jetptcut);
  jetTask->SetPercAreaCut(jetareacut);
  jetTask->SetLeadingHadronType(leadhadtype);
  // task specific setters
  jetTask->SetModulationFitType(fitType);
  jetTask->SetModulationFitOptions(fitOpts);
  jetTask->SetModulationFitMinMaxP(.001, 1);
  jetTask->SetCentralityClasses(centralities);

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

  switch (runMode) {
      case AliAnalysisTaskLocalRho::kLocal : {
          gStyle->SetOptFit(1);
          AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("good_fits_%s", name), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
          AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("bad_fits_%s", name),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							     Form("%s", AliAnalysisManager::GetCommonFileName()));
          mgr->ConnectOutput (jetTask, 2, coutput2);
          mgr->ConnectOutput (jetTask, 3, coutput3);
      } break;
      default: break;
  }
  return jetTask;
}
