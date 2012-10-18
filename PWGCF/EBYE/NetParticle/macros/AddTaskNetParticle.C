/* *********************************************************************************
 * File    : AddTaskNetParticle.C  
 * Author  : Jochen Thaeder <jochen@thaeder.de>
 * *********************************************************************************
 * Configuring NetParticle Task:
 * - ARGUMENTS : 
 *     name           -> Name of the task, containing partcile type :
 *                       Currently : Proton, Pion, Kaon      
 *     isModeDist     -> Fill Distributions
 *     isModeEff      -> Fill Efficiency/Contamination ThnSparse
 *     isModeDCA      -> Fill DCA ThnSparse
 *     useQAThnSparse -> Fill QA ThnSparse
 *     isCreateCSC    -> Prepare for CrossSectionCorrection 
 *                       - requires isModeEff to be set
 *                       - Proton only
 * 
 * - OUTPUT CONTAINER : #N = 5
 *   (1) - Standard Output, Distributions
 *   (2) - Efficiency ThnSparse
 *   (3) - Contamination ThnSparse
 *   (4) - DCA ThnSparse
 *   (5) - QA ThnSparse
 *
 ********************************************************************************* */

AliAnalysisTask *AddTaskNetParticle(const Char_t * name = "jthaeder_NetProton", 
				    Bool_t isModeDist, Bool_t isModeEff, Bool_t isModeDCA, Bool_t useQAThnSparse = kFALSE,
				    Bool_t isCreateCSC = kFALSE, Bool_t isModeAOD = kFALSE) {

  TString sName(name);

  if (isCreateCSC && !isModeEff) {
    Error("AddTaskNetParticle", "Creating CrossSectionCorrection needs 'isModeEff' to be set.");
    return NULL;
  }

  // ----------------------------------------------
  // -- Get the current analysis manager
  // ----------------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNetParticle", "No analysis manager found.");
    return NULL;
  }

  // ----------------------------------------------
  // -- Check for MC
  // ----------------------------------------------
  Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);
  if (isMC)
    Info("AddTaskNetParticle", "This task has MC.");
  
  // ----------------------------------------------
  // -- Create task configure it
  // ----------------------------------------------
  AliAnalysisTaskNetParticle *task = new AliAnalysisTaskNetParticle("AliAnalysisTaskNetParticle");
  if (!task) {
    Error("AddTaskNetParticle", "Task could not be created.");
    return NULL;
  }

  if (isMC) 
    task->SetIsMC();

  // -- Set particle type
  // ----------------------------------------------
  Float_t minPt, maxPt, minPtEff, maxPtEff, minPtForTOF; 

  if (sName.Contains("Proton")) {
    task->SetParticleSpecies(AliPID::kProton);
    task->SetControlParticleSpecies(3122, kTRUE, "Lambda");
    minPt    = 0.4;    maxPt    = 2.2;
    minPtEff = 0.2;    maxPtEff = 2.6;
    minPtForTOF = 0.8;
    if (isCreateCSC) {
      minPtForTOF = maxPtEff;
    }
  }
  else if (sName.Contains("Pion")) {
    task->SetParticleSpecies(AliPID::kPion);
    task->SetControlParticleSpecies(3122, kTRUE, "Lambda");  /// maybe something else ...
    minPt    = 0.3;    maxPt    = 0.9;
    minPtEff = 0.2;    maxPtEff = 1.2;
    minPtForTOF = 0.8;
  }
  else if (sName.Contains("Kaon")) {
    task->SetParticleSpecies(AliPID::kKaon);
    task->SetControlParticleSpecies(3122, kTRUE, "Lambda");  /// maybe something else ...
    minPt    = 0.2;    maxPt    = 0.4;
    minPtEff = 0.1;    maxPtEff = 2.5;
    minPtForTOF = 0.8;
  }
  else {
    Error("AddTaskNetParticle", "Unknown Particle type.");
    return NULL;
  }

  // -- Configure flags
  // ----------------------------------------------
  if (isModeEff) 
    task->SetModeEffCreation(1);     // => 1 = on    | 0 = off (default)
  if (isModeDCA)
    task->SetModeDCACreation(1);     // => 1 = on    | 0 = off (default)
  if (isModeDist)
    task->SetModeDistCreation(1);    // => 1 = on    | 0 = off (default)
  if(isModeAOD){
    task->SetIsAOD(1);               // => 1 = AODs  | 0 = ESDs
    task->SetTrackFilterBit(1024);   // 1024 = RAA cuts
  }
  // -- Enable QA plots
  if (useQAThnSparse)
    task->SetUseQATHnSparse(kTRUE);

  // -- Configure cuts 
  // ----------------------------------------------

  // -- Set cut flags ...
  task->SetESDTrackCutMode(0);     // => 0 = clean | 1 = dirty

  // -- Set standard event cuts
  task->SetVertexZMax(10.);   
  task->SetCentralityBinMax(7);

  // -- Set track event cuts
  task->SetRapidityMax(0.5); 
  task->SetMinTrackLengthMC(70.);  
  task->SetNSigmaMaxCdd(3.); 
  task->SetNSigmaMaxCzz(3.); 

  task->SetNSigmaMaxTPC(2.5);
  task->SetNSigmaMaxTOF(2.5);
  task->SetMinPtForTOFRequired(minPtForTOF);

  // -- Set analysis ranges
  task->SetEtaMax(0.9);        
  task->SetPtRange(minPt, maxPt);           // pt cut range for the analysis
  task->SetPtRangeEff(minPtEff, maxPtEff);  // pt cut range for the correction / efficiency / contamination creation
  
  // ----------------------------------------------
  // -- Initialize and add task to the ANALYSIS manager
  // ----------------------------------------------
  if (!task->Initialize())
    mgr->AddTask(task);
  else {
    Error("Initialize failed, not adding AliAnalysistaskNetParticle !!!");
    delete task;
    return NULL;
  }

  // ----------------------------------------------
  // -- data containers - input
  // ----------------------------------------------
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  // ----------------------------------------------
  // -- data containers - output
  // ----------------------------------------------
  AliAnalysisDataContainer *coutput     = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s.root", name));
  AliAnalysisDataContainer *coutputEff  = mgr->CreateContainer(Form("%s_eff", name), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s.root", name));
  AliAnalysisDataContainer *coutputCont = mgr->CreateContainer(Form("%s_cont", name), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s.root", name));
  AliAnalysisDataContainer *coutputDca  = mgr->CreateContainer(Form("%s_dca", name), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s.root", name));

  AliAnalysisDataContainer *coutputQA   = mgr->CreateContainer(Form("%sQA", name), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%sQA.root", name));
    
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput);
  mgr->ConnectOutput (task,  2, coutputEff);
  mgr->ConnectOutput (task,  3, coutputCont);
  mgr->ConnectOutput (task,  4, coutputDca);
  mgr->ConnectOutput (task,  5, coutputQA);

  return task;
}
