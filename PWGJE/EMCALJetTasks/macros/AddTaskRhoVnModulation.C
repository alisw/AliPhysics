/*
 * AddTask macro for class 
 * Redmer Alexander Bertens, rbertens@cern.ch
 * Utrecht University, Utrecht, Netherlands
 *
 * Note: this macro is pretty much a copy of AddTaskEmcalJetSample.C
 *
 */

AliAnalysisTaskRhoVnModulation* AddTaskRhoVnModulation(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Double_t   jetradius          = 0.2,
  Double_t   jetptcut           = 1,
  Double_t   jetareacut         = 0.557,
  UInt_t     type               = AliAnalysisTaskEmcal::kTPC,
  Int_t      leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskRhoVnModulation",
  UInt_t     runMode            = AliAnalysisTaskRhoVnModulation::kGrid,
  Bool_t     fillQA             = kTRUE,
  TString    fitOpts            = "LWQIM",
  UInt_t     fitType            = AliAnalysisTaskRhoVnModulation::kFourierSeries,
  TArrayI    *centralities      = 0x0,
  TRandom3   *randomizer        = 0x0
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

  AliAnalysisTaskRhoVnModulation* jetTask = new AliAnalysisTaskRhoVnModulation(name, runMode);
  // inherited setters
  jetTask->SetAnaType(type);
  jetTask->SetTracksName(ntracks);
  jetTask->SetClusName(nclusters);
  jetTask->SetJetsName(njets);
  jetTask->SetRhoName(nrho);
  jetTask->SetJetRadius(jetradius);
  jetTask->SetJetPtCut(jetptcut);
  jetTask->SetPercAreaCut(jetareacut);
  jetTask->SetLeadingHadronType(leadhadtype);
  // task specific setters
  jetTask->SetFillQAHistograms(fillQA);
  jetTask->SetDebugMode(-1);
  jetTask->SetModulationFitType(fitType);
  jetTask->SetModulationFitOptions(fitOpts);
  jetTask->SetModulationFitMinMaxP(.05, 1);
  // if centralities haven't been specified use defaults
  if(!centralities) {
     Int_t c[] = {0, 10, 30, 50, 70, 90};
     jetTask->SetCentralityClasses(new TArrayI(sizeof(c)/sizeof(c[0]), c));
  }
  // if a randomized hasn't specified use a safe default 
  if(!randomizer) jetTask->SetRandomSeed(new TRandom3(0));




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
      case AliAnalysisTaskRhoVnModulation::kLocal : {
          gStyle->SetOptFit(1);
          AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("good_fits_%s", name.Data()), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
          AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("bad_fits_%s", name.Data()),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							     Form("%s", AliAnalysisManager::GetCommonFileName()));
          mgr->ConnectOutput (jetTask, 2, coutput2);
          mgr->ConnectOutput (jetTask, 3, coutput3);
      } break;
      default: break;
  }
  return jetTask;
}
