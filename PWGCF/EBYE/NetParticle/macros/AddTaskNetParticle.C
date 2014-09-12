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
 *     isModeQA       -> Fill QA ThnSparse
 *     isModeAOD      -> Use AOD input 
 *     isCreateCSC    -> Prepare for CrossSectionCorrection 
 *                       - requires isModeEff to be set
 *                       - Proton only
 *     isSetExt       -> 1 if want to set pt, nSigma from arguments
 *                       0 takes default values 
 *
 *     modeCSC        -> Use differnt Pt cut for
 *                       0 : TPC+TOF
 *                       1 : TPC
 *     modeCuts       -> Different Cut scenarios
 *                       0 : Standard cuts
 *                       1 : LF cuts
 *     modePID        -> PID Strategy 
 *                      -1 :   Default -> 7 (modeCuts=0) and 5 (modeCuts=1)  
 *                       0 :   TPC(TPClow+TPCHigh)
 *                       1 :   ITS
 *                       2 :   TOF
 *                       3 :   ITS+TPC(TPClow+TPCHigh)
 *                       4 :   TPC(TPClow+TPCHigh)+TOF
 *                       5 :   TPC(TPClow+TPCHigh)+TOF for pT >= fMinPtForTOFRequired TOF is required, below, only used if there
 *                       6 :   TPC(TPClow+TPCHigh)+ITS+TOF with TOF only for those tracks which have TOF information
 *                       7 :   TPC(TPClow+TPCHigh)+ITS+TOF for pT >= fMinPtForTOFRequired TOF is required, below, only used if there
 *                       8 :   TPC(TPClow+TPCHigh)+ITS+TOF 
 * 
 * *********************************************************************************
 * - PID Strategy
 *
* *********************************************************************************
 * - OUTPUT CONTAINER : #N = 5
 *   (1) - Standard Output, Distributions
 *   (2) - Efficiency ThnSparse
 *   (3) - Contamination ThnSparse
 *   (4) - DCA ThnSparse
 *   (5) - QA ThnSparse
 *
 ********************************************************************************* */

AliAnalysisTask *AddTaskNetParticle(const Char_t *name      = "ITS_NetProton", 
				    Bool_t  isModeDist      = kTRUE, 
				    Bool_t  isModeEff       = kFALSE, 
				    Bool_t  isModeDCA       = kFALSE,
				    Bool_t  isModeQA        = kFALSE,
				    Bool_t  isCreateCSC     = kFALSE, 
				    Bool_t  isModeAOD       = kFALSE,
				    Bool_t  isSetExt        = kFALSE,
				    Int_t   aodFilterBit    = 1024, /* 1024 = RAA cuts */
				    Int_t   modeCSC         = 0, 
				    Int_t   modeCuts        = 0, 
				    Int_t   modePID         =-1,
				    Float_t gMinPt          = 0.3,
				    Float_t gMaxPt          = 2.5,
				    Float_t gMinPtForTof    = 0.21,
				    Float_t gMaxPtForTPClow = 0.69,
				    Float_t gMinPtEff       = 0.3,    
				    Float_t gMaxPtEff       = 2.5,
				    Float_t gSigmaITS       = 4.0,
				    Float_t gSigmaTOF       = 4.0,
				    Float_t gSigmaTPC       = 4.0,   
				    Float_t gSigmaTPClow    = 3.0) {

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
  if (isMC) 
    task->SetIsMC();
  if (isModeEff) 
    task->SetModeEffCreation(1);             // => 1 = on    | 0 = off (default)
  if (isModeDCA)
    task->SetModeDCACreation(1);             // => 1 = on    | 0 = off (default)
  if (isModeDist)
    task->SetModeDistCreation(1);            // => 1 = on    | 0 = off (default)
  if (isModeAOD) {
    task->SetIsAOD(1);                       // => 1 = AOD   | 0 = ESD (default)
    task->SetTrackFilterBit(aodFilterBit);   
  }
  if (isModeQA)
    task->SetModeQACreation(1);              // => 1 = on    | 0 = off (default)
  
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
  Float_t minPt, maxPt, minPtEff, maxPtEff, minPtForTOF, etaMax, etaMaxEff, nSigmaITS, nSigmaTPC, nSigmaTPClow, nSigmaTOF, maxPtForTPClow; 
  Int_t pidStrategy;


if (sName.Contains("Proton")) {
    helper->SetParticleSpecies(AliPID::kProton);
    minPtForTOF    = 0.69;                //    minPtForTOF = 0.21;
    maxPtForTPClow = 0.69;
    minPt          = 0.5;     maxPt   = 2.0;    //    minPt    = 0.22;     maxPt   = 4.5;
    minPtEff       = 0.3;    maxPtEff = 2.5;    //    minPtEff = 0.22;     maxPtEff = 4.5;
    
    etaMax         = 0.8;  
    etaMaxEff      = 0.8;  
    nSigmaITS      = 4.0;   nSigmaTPC = 4.0;   nSigmaTPClow = 3.0;   nSigmaTOF = 4.0; 

    // For TPC only case
    if (isCreateCSC && modeCSC == 1)
      minPtForTOF = maxPtEff;
  }
  else if (sName.Contains("Pion")) {
    helper->SetParticleSpecies(AliPID::kPion);
    minPtForTOF    = 0.3;
    maxPtForTPClow = 0.69;
    minPt          = 0.3;    maxPt    = 1.5;
    minPtEff       = 0.3;    maxPtEff = 2.5;
  
    etaMax         = 8.8;    
    etaMaxEff      = 8.8;  
    nSigmaITS      = 4.0;   nSigmaTPC = 4.0;   nSigmaTPClow = 3.0;   nSigmaTOF = 4.0;   
    pidStrategy    = 1;
  }
  else if (sName.Contains("Kaon")) {
    helper->SetParticleSpecies(AliPID::kKaon);
    minPtForTOF    = 0.3;
    maxPtForTPClow = 0.69;
    minPt          = 0.3;    maxPt    = 1.5;
    minPtEff       = 0.3;    maxPtEff = 2.5;
 
    etaMax         = 0.8;   
    etaMaxEff      = 0.8; 
    nSigmaITS      = 4.0;   nSigmaTPC = 4.0;   nSigmaTPClow = 3.0;   nSigmaTOF = 4.0;      
    pidStrategy    = 1;
  }
  else if (sName.Contains("Charge")) {
    helper->SetUsePID(kFALSE);
    minPtForTOF    = 0.1;
    maxPtForTPClow = 0.1;
    minPt          = 0.1;    maxPt    = 2.5;
    minPtEff       = 0.1;    maxPtEff = 3.0;

    etaMax         = 0.8; 
    etaMaxEff      = 0.8; 
    nSigmaITS      = -1.;   nSigmaTPC = -1.;   nSigmaTOF = -1.;
    pidStrategy    = 1;
  }
  else {
    Error("AddTaskNetParticle", "Unknown Particle type.");
    delete task;
    return NULL;
  }
 
// ----------------------------------------------
// -- use value arguments --
// ----------------------------------------------
 if (isSetExt) {
   minPt          = gMinPt;    
   maxPt          = gMaxPt;
   minPtForTOF    = gMinPtForTof;     
   maxPtForTPClow = gMaxPtForTPClow;
   minPtEff       = gMinPtEff;
   maxPtEff       = gMaxPtEff;
   
   nSigmaITS      = gSigmaITS;   
   nSigmaTPC      = gSigmaTPC;   
   nSigmaTPClow   = gSigmaTPClow;  
   nSigmaTOF      = gSigmaTOF;
    
 }
 // ----------------------------------------------
 // -- PID Strategy
 // ----------------------------------------------
 
  if (modePID     == -1) { // default
    pidStrategy   =   7;         // 7: ITS + TPC + TOF (using minPtForTOF)
    if (modeCuts  ==  1)
      pidStrategy =   5;       // 5: TPC + TOF (using minPtForTOF) 
  }
  else
    pidStrategy   = modePID;
  
  // ----------------------------------------------
  // -- Read Environment Variables 
  // ----------------------------------------------
  ifstream in;
  in.open("setRunENV.txt");
  
  TString current;    
  while(in.good()) {
    in >> current;
    
    TObjArray *arr = current.Tokenize('=');
    if (!arr) 
      continue;
    
    TObjString* oKey = dynamic_cast<TObjString*>(arr->At(0));
    TObjString* oValue = dynamic_cast<TObjString*>(arr->At(1));
    if (!oKey)
      continue;

    TString key(oKey->GetString());
    TString value(oValue->GetString());
   
    if (!key.CompareTo("NETPARTICLE_PID_STRATEGY")) {
      pidStrategy = value.Atoi();
      printf(">>>> USE NETPARTICLE_PID_STRATEGY %d\n", pidStrategy);
    }
    if (!key.CompareTo("NETPARTICLE_NSIGMAMAX_ITS")) {
      nSigmaITS = value.Atof();
      printf(">>>> USE NETPARTICLE_NSIGMAMAX_ITS %.2f\n", nSigmaITS);
    }
    if (!key.CompareTo("NETPARTICLE_NSIGMAMAX_TPC")) {
      nSigmaTPC = value.Atof();
      printf(">>>> USE NETPARTICLE_NSIGMAMAX_TPC %.2f\n", nSigmaTPC);
    }
    if (!key.CompareTo("NETPARTICLE_NSIGMAMAX_TOF")) {
      nSigmaTOF = value.Atof();
      printf(">>>> USE NETPARTICLE_NSIGMAMAX_TOF %.2f\n", nSigmaTOF);
    }

    arr->Clear();
    delete arr;
  }

  in.close();

  // ----------------------------------------------
  // -- Configure cuts 
  // ----------------------------------------------

  // -- Set cut flags 
  task->SetESDTrackCutMode(modeCuts);       // => 0 = normal | 1 = LF

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
  //  helper->SetRapidityMax(0.2); 
  helper->SetMinTrackLengthMC(70.);  
  helper->SetNSigmaMaxCdd(0.);    //  3. ||   ->> Turn off sigmaDCA cuts for now
  helper->SetNSigmaMaxCzz(0.);    //  3. ||   ->> Turn off sigmaDCA cuts for now
  helper->SetPhiRange(0., 3.88);  //  Only used if requested in task - default is TwoPi

  // -- Set pid cuts
  helper->SetPIDStrategy(pidStrategy);
  helper->SetNSigmaMaxITS(nSigmaITS);
  helper->SetNSigmaMaxTPC(nSigmaTPC);
  helper->SetNSigmaMaxTPClow(nSigmaTPClow);
  helper->SetNSigmaMaxTOF(nSigmaTOF);
  helper->SetMinPtForTOFRequired(minPtForTOF);
  helper->SetMaxPtForTPClow(maxPtForTPClow);
  
  // -- Set N sub samples
  helper->SetNSubSamples(20);

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
