Float_t kPtTrackMin = 0.15; // pt cut on tracks
Int_t kDetectorType = 1;    // 0 = MC, 1 = TPC

AliAnalysisDataContainer*  cout_jet_cont = 0x0;

void SetJetReaderExchangeContainer(AliAnalysisDataContainer* acont) {
  cout_jet_cont = acont;

}

AliAnalysisDataContainer* GetJetReaderExchangeContainer() { return cout_jet_cont;}
AliAnalysisDataContainer* AddJetExchangeContainer(const char* contname) { return AliAnalysisManager::GetAnalysisManager()->CreateContainer(contname, TTree::Class(),AliAnalysisManager::kExchangeContainer);}

AliJetReader *CreateJetReader(AliAnalysisDataContainer* contname,Char_t *jr,UInt_t filterMask); // Common config

AliAnalysisTaskJetsReader *AddTaskJetsReader(AliAnalysisDataContainer* contname,Char_t *jr, UInt_t filterMask = 0,Float_t ptTrackMin = 0.15 , Int_t dettype = 1); // for the new AF
Int_t AddTaskJetsReaderDelta(AliAnalysisDataContainer* contname,char *nonStdFile = "",UInt_t filterMask = 0,Bool_t kUseAODMC = kTRUE,UInt_t runFlag = 1|4|32|128|256);     
AliAnalysisTaskJetsReader *AddTaskJetsReader(AliAnalysisDataContainer* contname,UInt_t filterMask = 0);

AliAnalysisTaskJetsReader *AddTaskJetsReader(const char* contname, UInt_t filterMask = 0) { return AddTaskJetsReader(AddJetExchangeContainer(contname),filterMask); } // LEGO trains  
AliAnalysisTaskJetsReader *AddTaskJetsReader(const char* contname,Char_t *jr, UInt_t filterMask = 0,Float_t ptTrackMin = 0.15 , Int_t dettype = 1) {  
 return AddTaskJetsReader(AddJetExchangeContainer(contname),jr,filterMask,ptTrackMin,dettype); 

}

AliAnalysisTaskJetsReader *AddTaskJetsReader(AliAnalysisDataContainer* contname,UInt_t filterMask){
  // fills the standard input "jets" branch in the AOD
  // need the ESDFilter to run before, to access the AODtracks
  // Tracks selected by the first Filter (1<<0)
  // needs to be adapted for different cuts
  

  return AddTaskJetsReader(contname,"AOD",filterMask);

}



Int_t AddTaskJetsReaderDelta(AliAnalysisDataContainer* contname,char *nonStdFile,UInt_t filterMask,Bool_t kUseAODMC,UInt_t runFlag){

  // Adds a whole set of jet finders  all to be written
  // to a delta AOD
  
  // this can in principle be done also with on the fly 
  // if we have an ESD input jet fidner task does automatically fetch the ouput aod

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetsReaderDelta", "No analysis manager to connect to.");
    return -1;
  }  
  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetsReaderDelta", "This task requires an input event handler");
    return -1;
  }


  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskJetsReaderDelta", "This task needs an output event handler");
    return -1;
  }

  TString outFile(nonStdFile);


  AliAnalysisTaskJetsReader *taskjetsReader = 0;
  Int_t iCount = 0;

  const char *cReader[4] = {"AOD","AODMC","AODMC2","AODMC2b"};  

  

  // Ligne a revoir
  //  if(!kUseAODMC)ib=1; // switch OFF MC if we do not have it
  for(int ib = 0;ib<4;ib++){      
    if(1<<ib){
      cout << "cReader[" << ib << "]: " << cReader[ib] << ", filterMask: " << filterMask << endl;
      taskjetsReader = AddTaskJetsReader(contname,cReader[ib],filterMask);
      if(taskjetsReader){
	if(outFile.Length()>0)taskjetsReader->SetNonStdOutputFile(outFile.Data());
      }
    }
  }
  return 0;
}

AliAnalysisTaskJetsReader *AddTaskJetsReader(AliAnalysisDataContainer* contname,Char_t *jr, UInt_t filterMask,Float_t ptTrackMin,Int_t dettype)
{
  // Creates a jet finder task, configures it and adds it to the analysis manager.

  kPtTrackMin = ptTrackMin;
  kDetectorType = dettype;

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetsReader", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetsReader", "This task requires an input event handler");
    return NULL;
  }

  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskJetsReader", "This task needs an output event handler");
    return NULL;
  }   
  

  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskJetsReader *taskjetsReader;
  AliJetReader *er = CreateJetReader(contname,jr,filterMask);

  taskjetsReader = new AliAnalysisTaskJetsReader(Form("JetAnalysis%s",jr));

  // Connect jet finder to task.
  taskjetsReader->SetJetReader(er);
  taskjetsReader->SetConfigFile("");
  taskjetsReader->SetDebugLevel(-1);
  mgr->SetGlobalStr(cout_jet_cont->GetName(),Form("%s#_Filter%05d_Cut%05d",jr,filterMask,(int)((1000.*ptTrackMin))));
  mgr->AddTask(taskjetsReader);

  cout_jet_cont = GetJetReaderExchangeContainer();
  if(! cout_jet_cont) { ::Error("AddTaskJetsReader", "This task needs an Exchange container");  return NULL;}
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (taskjetsReader, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (taskjetsReader, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (taskjetsReader, 1, cout_jet_cont); 
   

  return taskjetsReader;
}


AliJetReader *CreateJetReader(AliAnalysisDataContainer* contname,Char_t *jr,UInt_t filterMask){
  AliJetReader *er = 0;
  SetJetReaderExchangeContainer(contname);
  Float_t ptCut  = kPtTrackMin ; // cut on track p_T
 
  switch (jr) {
  case "ESD":
    AliJetReaderHeader *jrh = new AliJetReaderHeader(kDetectorType);
    jrh->SetDataType("ESD");
    jrh->SetComment("ESD Reader");
    jrh->SetDebug(-1);
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-0.9,0.9);
    jrh->SetFiducialPhi(0,2*TMath::Pi());
    // Define reader and set its header                                     
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;
  case "AOD":
    AliJetReaderHeader *jrh = new AliJetReaderHeader(kDetectorType);
    jrh->SetDataType("AOD");
    jrh->SetComment("AOD Reader");
    jrh->SetDebug(-1);
    jrh->SetPtCut(ptCut);
    jrh->SetTestFilterMask(32); // Change this one for a different set of cuts
    if(filterMask>0)jrh->SetTestFilterMask(filterMask); 
    jrh->SetFilterType(0); // Filter type
    jrh->SetFiducialEta(-0.9,0.9);
    jrh->SetFiducialPhi(0,2*TMath::Pi());
    // Define reader and set its header
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMC":
    AliJetReaderHeader *jrh = new AliJetReaderHeader(/*kDetectorType =*/0);
    jrh->SetComment("AOD MC Reader");
    jrh->SetDataType("AODMC");
    jrh->SetDebug(-1);
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9
    // Define reader and set its header
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMC2":
    AliJetReaderHeader *jrh = new AliJetReaderHeader(/*kDetectorType =*/0);
    jrh->SetDataType("AODMC2");
    jrh->SetComment("AOD MC Reader");
    jrh->SetDebug(-1);
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9
    // Define reader and set its header
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMC2b":
    AliJetReaderHeader *jrh = new AliJetReaderHeader(/*kDetectorType =*/0);
    jrh->SetDataType("AODMC2b");
    jrh->SetComment("AOD MC Reader");
    jrh->SetDebug(-1);
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-0.9,0.9); 
    // Define reader and set its header
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODextra":
    AliJetReaderHeader *jrh = new AliJetReaderHeader(kDetectorType);
    jrh->SetDataType("AODextra");
    jrh->SetComment("AOD Reader with extra branch");
    jrh->SetDebug(-1);
    jrh->SetPtCut(ptCut);
    jrh->SetTestFilterMask(32); // Change this one for a different set of cuts
    if(filterMask>0)jrh->SetTestFilterMask(filterMask);
    // Define reader and set its header
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODextraonly":
    AliJetReaderHeader *jrh = new AliJetReaderHeader(kDetectorType);
    jrh->SetDataType("AODextraonly");
    jrh->SetComment("AOD Reader with extra branch");
    jrh->SetDebug(-1);
    jrh->SetPtCut(ptCut);
    jrh->SetTestFilterMask(32); // Change this one for a different set of cuts
    if(filterMask>0)jrh->SetTestFilterMask(filterMask);
    // Define reader and set its header
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;

  case "MC":
    AliJetKineReaderHeader *jrh = new AliJetKineReaderHeader();
    jrh->SetDataType("MC");
    jrh->SetComment("MC full Kinematics");
    jrh->SetDebug(-1);
    jrh->SetFastSimTPC(kFALSE);
    jrh->SetFastSimEMCAL(kFALSE);
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9                                                                             
    // Define reader and set its header                                     
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;
  case "MC2":
    AliJetKineReaderHeader *jrh = new AliJetKineReaderHeader();
    jrh->SetDataType("MC2");
    jrh->SetComment("MC full Kinematics spearate config charged only");
    jrh->SetDebug(-1);
    jrh->SetFastSimTPC(kFALSE);
    jrh->SetFastSimEMCAL(kFALSE);
    jrh->SetChargedOnly(kTRUE);
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0 .9                                                                             
    // Define reader and set its header                                     
    er = new AliJetReader();
    er->SetReaderHeader(jrh);
    break;

  default:
    ::Error("AddTaskJetsReader", Form("Wrong jet reader selected: %s\n",jr));
    return 0x0;
    break;
  }

  return er;

}
