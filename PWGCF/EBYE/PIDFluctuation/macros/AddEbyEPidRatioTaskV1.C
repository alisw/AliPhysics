/* *************************************************************************
 *             AliEbyE Analysis for Particle Ratio Fluctuation             *
 *                   Deepika Rathee  | Satyajit Jena                       *
 *                   drathee@cern.ch | sjena@cern.ch                       *
 *                  Date: Wed Jul  9 18:38:30 CEST 2014                    * 
 *          New approch to find particle ratio to reduce memory            *
 *                             (Test Only)                                 *
 ***************************************************************************/

AliAnalysisTask *AddEbyEPidRatioTaskV1(const Char_t *name = "NuDyn",   //  0                     
				       Int_t  isModeDist  = 1,         //  1
				       Bool_t isModeEff   = 0,         //  2
				       Bool_t isModeDCA   = 0,         //  3
				       Bool_t isModeQA    = 0,         //  4
				       Int_t  isRatio     = 3,         //  5
				       Bool_t isModeAOD   = 0,         //  6
				       Bool_t isSetExt    = 0,         //  7 
				       
				       Int_t   aodFilterBit = 1024,    //  8
				       Float_t gEta         = 0.8,     //  9
				       Int_t   modeCuts     = 0,       // 10
				       Int_t   modePID      =-1,       // 11 
				       Float_t gMinPt       = 0.3,     // 12
				       Float_t gMaxPt       = 2.5,     // 13
				       Float_t gMinPtForTof = 0.21,    // 14
				       Float_t gMaxPtForTPClow = 0.69, // 15 
				       Float_t gY           = 0.5,     // 16
				       Float_t gMaxPtEff    = 2.5,     // 17
				       Float_t gSigmaITS    = 4.0,     // 18
				       Float_t gSigmaTOF    = 4.0,     // 19
				       Float_t gSigmaTPC    = 4.0,     // 20
				       Float_t gSigmaTPClow = 3.0,
                                       Bool_t  isPer        = 0,
                                       Int_t   nSub         = 25) {   // 21

  
  TString sName(name);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNetParticle", "No analysis manager found.");
    return NULL;
  }
  
  Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);
  if (isMC)
    Info("AddTaskNetParticle", "This task has MC.");
  
  AliEbyEPidRatioTask *task = new AliEbyEPidRatioTask("EbyEPidRatio");
  if (!task) {
    Error("EbyEPidRatio", "Task could not be created.");
    return NULL;
  }

  if (isMC) 
    task->SetIsMC();
  if (isModeEff) 
    task->SetModeEffCreation(1);             // => 1 = on    | 0 = off (default)
  if (isModeDCA)
    task->SetModeDCACreation(1);             // => 1 = on    | 0 = off (default)
 
 
  if (isModeAOD) {
    task->SetIsAOD(1);                       // => 1 = AOD   | 0 = ESD (default)
    task->SetTrackFilterBit(aodFilterBit);   
  }
  if (isModeQA)
    task->SetModeQACreation(1);              // => 1 = on    | 0 = off (default)

  task->SetModeDistCreation(isModeDist);            // => 1 = on    | 0 = off (default)
  task->SetIsRatio(isRatio);

  Float_t minPt,     maxPt,     minPtEff,     maxPtEff,  minPtForTOF;
  Float_t nSigmaITS, nSigmaTPC, nSigmaTPClow, nSigmaTOF, maxPtForTPClow; 
  Float_t etaMax,    etaMaxEff, maxRap; 
  Int_t   pidStrategy;

  minPtForTOF    = 0.69;   
  maxPtForTPClow = 0.69;
  minPt          = 0.5;     
  maxPt          = 2.0; 
  minPtEff       = 0.3;    
  maxPtEff       = 2.5; 
  
  maxRap         = 0.5;
  etaMax         = 0.8;  
 
  etaMaxEff      = 0.8;  
  nSigmaITS      = 4.0;   
  nSigmaTPC      = 4.0;   
  nSigmaTPClow   = 3.0;   
  nSigmaTOF      = 4.0; 
    
  if (isSetExt) {
    minPt          = gMinPt;    
    maxPt          = gMaxPt;
    minPtForTOF    = gMinPtForTof;     
    maxPtForTPClow = gMaxPtForTPClow;
    minPtEff       = gMinPt;
    maxPtEff       = gMaxPtEff;
    
    maxRap         = gY;
    etaMax         = gEta;  

    nSigmaITS      = gSigmaITS;   
    nSigmaTPC      = gSigmaTPC;   
    nSigmaTPClow   = gSigmaTPClow;  
    nSigmaTOF      = gSigmaTOF;
  }
 
  if (modePID == -1) { // default
    pidStrategy   = 7;         // 7: ITS + TPC + TOF (using minPtForTOF)
    if (modeCuts == 1)
      pidStrategy = 5;       // 5: TPC + TOF (using minPtForTOF) 
  }
  else
    pidStrategy = modePID;
  
  AliEbyEPidRatioHelper *helper = new AliEbyEPidRatioHelper;
  if (!helper) {
    Error("AddTaskNetParticle", "Helper could not be created.");
    delete task;
    return NULL;
  }
  

  task->SetESDTrackCutMode(modeCuts);       // => 0 = normal | 1 = LF
  task->SetEtaMax(etaMax);                  // eta cut
  task->SetEtaMaxEff(etaMaxEff);            // eta cut for efficiency
  task->SetPtRange(minPt, maxPt);           // pt cut range for the analysis
  task->SetPtRangeEff(minPtEff, maxPtEff);  // pt cut range for the correction / efficiency / contamination creation
  if (isPer) task->SetIsPer();  // 
  helper->SetVertexZMax(10.);   
  helper->SetCentralityBinMax(11);
  helper->SetRapidityMax(maxRap); 
  helper->SetMinTrackLengthMC(70.);  
  helper->SetNSigmaMaxCdd(0.);    //  3. ||   ->> Turn off sigmaDCA cuts for now
  helper->SetNSigmaMaxCzz(0.);    //  3. ||   ->> Turn off sigmaDCA cuts for now
  helper->SetPhiRange(0., 3.88);  //  Only used if requested in task - default is TwoPi
  helper->SetPIDStrategy(pidStrategy);
  helper->SetNSigmaMaxITS(nSigmaITS);
  helper->SetNSigmaMaxTPC(nSigmaTPC);
  helper->SetNSigmaMaxTPClow(nSigmaTPClow);
  helper->SetNSigmaMaxTOF(nSigmaTOF);
  helper->SetMinPtForTOFRequired(minPtForTOF);
  helper->SetMaxPtForTPClow(maxPtForTPClow);
  helper->SetNSubSamples(nSub);
  task->SetNetParticleHelper(helper);
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString outputFileName   = Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name);
  TString outputQAFileName = Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name);
    
  AliAnalysisDataContainer *coutput     = mgr->CreateContainer(name, TList::Class(),  AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputEff  = mgr->CreateContainer(Form("%s_eff",  name), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputCont = mgr->CreateContainer(Form("%s_cont", name), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputDca  = mgr->CreateContainer(Form("%s_dca",  name), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

  AliAnalysisDataContainer *coutputQA   = mgr->CreateContainer(Form("%sQA",    name), TList::Class(), AliAnalysisManager::kOutputContainer, outputQAFileName);
    
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput);
  mgr->ConnectOutput (task,  2, coutputEff);
  mgr->ConnectOutput (task,  3, coutputCont);
  mgr->ConnectOutput (task,  4, coutputDca);
  mgr->ConnectOutput (task,  5, coutputQA);

  return task;
}
