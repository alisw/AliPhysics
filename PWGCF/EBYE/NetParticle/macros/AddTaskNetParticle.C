/* *********************************************************************************
 * File    : AddTaskNetParticle.C  
 * Authors : Jochen Thaeder <jochen@thaeder.de>
 *           Michael Weber <m.weber@cern.ch>
 * *********************************************************************************
 * Configuring NetParticle Task:
 * - ARGUMENTS : 
 *     name           -> Name of the task, containing particle type :
 *                       Currently : Proton, Pion, Kaon, Charge      
 *     isModeDist     -> Fill Distributions
 *     isModeEff      -> Fill Efficiency/Contamination ThnSparse
 *     isModeDCA      -> Fill DCA ThnSparse
 *     useQAThnSparse -> Fill QA ThnSparse
 *     isCreateCSC    -> Prepare for CrossSectionCorrection 
 *                       - requires isModeEff to be set
 *                       - Proton only
 *     isModeAOD      -> Use AOD input 
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
  // -- Create task 
  // ----------------------------------------------
  AliAnalysisTaskNetParticle *task = new AliAnalysisTaskNetParticle("AliAnalysisTaskNetParticle");
  if (!task) {
    Error("AddTaskNetParticle", "Task could not be created.");
    return NULL;
  }

  // ----------------------------------------------
  // -- Configure flags
  // ----------------------------------------------

  // -- Enable QA plots
  Int_t isModeQA = 0;

  if (isMC) 
    task->SetIsMC();
  if (isModeEff) 
    task->SetModeEffCreation(1);     // => 1 = on    | 0 = off (default)
  if (isModeDCA)
    task->SetModeDCACreation(1);     // => 1 = on    | 0 = off (default)
  if (isModeDist)
    task->SetModeDistCreation(1);    // => 1 = on    | 0 = off (default)
  if (isModeAOD) {
    task->SetIsAOD(1);               // => 1 = AOD   | 0 = ESD (default)
    task->SetTrackFilterBit(1024);   // 1024 = RAA cuts
  }
  if (isModeQA)
    task->SetModeQACreation(1);      // => 1 = on    | 0 = off (default)

  // ----------------------------------------------
  // -- Create helper class
  // ----------------------------------------------
  AliAnalysisNetParticleHelper *helper = new AliAnalysisNetParticleHelper;
  if (!helper) {
    Error("AddTaskNetParticle", "Helper could not be created.");
    delete task;
    return NULL;
  }

  task->SetNetParticleHelper(helper);

  // ----------------------------------------------
  // -- Set particle type
  // ----------------------------------------------
  Float_t minPt, maxPt, minPtEff, maxPtEff, minPtForTOF, etaMax, etaMaxEff, nSigmaTPC, nSigmaTOF; 

  if (sName.Contains("Proton")) {
    helper->SetParticleSpecies(AliPID::kProton);
    minPt    = 0.4;    maxPt    = 2.2;
    minPtEff = 0.2;    maxPtEff = 2.6;
    minPtForTOF = 0.8;
    etaMax     = 99.99;    // 0.8 ->> eta cut off for now
    etaMaxEff  = 99.99;    // 0.9 ->> eta cut off for now
    nSigmaTPC = 2.5;   nSigmaTOF = 2.5;
    if (isCreateCSC) {
      minPtForTOF = maxPtEff;
    }
  }
  else if (sName.Contains("Pion")) {
    helper->SetParticleSpecies(AliPID::kPion);
    minPt    = 0.25;   maxPt    = 0.7;
    minPtEff = 0.2;    maxPtEff = 1.2;
    minPtForTOF = 0.8;
    etaMax     = 99.99;    // 0.8 ->> eta cut off for now
    etaMaxEff  = 99.99;    // 0.9 ->> eta cut off for now
    nSigmaTPC = 2.5;   nSigmaTOF = 2.5;
  }
  else if (sName.Contains("Kaon")) {
    helper->SetParticleSpecies(AliPID::kKaon);
    minPt    = 0.5;    maxPt    = 1.4;
    minPtEff = 0.1;    maxPtEff = 2.5;
    minPtForTOF = 0.5;
    etaMax     = 99.99;    // 0.8 ->> eta cut off for now
    etaMaxEff  = 99.99;    // 0.9 ->> eta cut off for now
    nSigmaTPC = 2.5;   nSigmaTOF = 2.5;
  }
  else if (sName.Contains("Charge")) {
    helper->SetUsePID(kFALSE);
    minPt    = 0.3;    maxPt    = 2.5;
    minPtEff = 0.1;    maxPtEff = 3.0;
    minPtForTOF = -1.;
    etaMax     = 0.8; 
    etaMaxEff  = 0.9; 
    nSigmaTPC = -1.;   nSigmaTOF = -1.;
  }
  else {
    Error("AddTaskNetParticle", "Unknown Particle type.");
    delete task;
    return NULL;
  }

  // ----------------------------------------------
  // -- Configure cuts 
  // ----------------------------------------------

  // -- Set cut flags 
  task->SetESDTrackCutMode(0);              // => 0 = clean | 1 = dirty

  // -- Set analysis ranges
  task->SetEtaMax(etaMax);                  // eta cut
  task->SetEtaMaxEff(etaMaxEff);            // eta cut for efficiency
  task->SetPtRange(minPt, maxPt);           // pt cut range for the analysis
  task->SetPtRangeEff(minPtEff, maxPtEff);  // pt cut range for the correction / efficiency / contamination creation

  // ----------------------------------------------
  // -- Configure cuts - helper class
  // ----------------------------------------------

  // -- Set standard event cuts
  helper->SetVertexZMax(10.);   
  helper->SetCentralityBinMax(7);

  // -- Set track event cuts
  helper->SetRapidityMax(0.5); 
  helper->SetMinTrackLengthMC(70.);  
  helper->SetNSigmaMaxCdd(0.);    //  3. ||   ->> Turn off sigmaDCA cuts for now
  helper->SetNSigmaMaxCzz(0.);    //  3. ||   ->> Turn off sigmaDCA cuts for now
  helper->SetPhiRange(0., 4.);

  // -- Set pid cuts
  helper->SetNSigmaMaxTPC(nSigmaTPC);
  helper->SetNSigmaMaxTOF(nSigmaTOF);
  helper->SetMinPtForTOFRequired(minPtForTOF);
  
  // ----------------------------------------------
  // -- Add task to the ANALYSIS manager
  // ----------------------------------------------
  mgr->AddTask(task);

  // ----------------------------------------------
  // -- data containers - input
  // ----------------------------------------------
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  // ----------------------------------------------
  // -- data containers - output
  // ----------------------------------------------
  TString outputFileName = "";
  TString outputQAFileName = "";
  if(isModeAOD){
    outputFileName = Form("%s:%s",AliAnalysisManager::GetCommonFileName(),name);
    outputQAFileName = Form("%s:%s",AliAnalysisManager::GetCommonFileName(),name);
  }
  else{
    outputFileName = Form("%s.root",name);
    outputQAFileName = Form("%sQA.root",name);
  }

  AliAnalysisDataContainer *coutput     = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputEff  = mgr->CreateContainer(Form("%s_eff", name), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputCont = mgr->CreateContainer(Form("%s_cont", name), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputDca  = mgr->CreateContainer(Form("%s_dca", name), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

  AliAnalysisDataContainer *coutputQA   = mgr->CreateContainer(Form("%sQA", name), TList::Class(), AliAnalysisManager::kOutputContainer, outputQAFileName);
    
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput);
  mgr->ConnectOutput (task,  2, coutputEff);
  mgr->ConnectOutput (task,  3, coutputCont);
  mgr->ConnectOutput (task,  4, coutputDca);
  mgr->ConnectOutput (task,  5, coutputQA);

  return task;
}
