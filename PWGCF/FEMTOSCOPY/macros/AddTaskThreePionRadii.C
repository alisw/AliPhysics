AliThreePionRadii *AddTaskThreePionRadii(
				 Bool_t LEGO=kTRUE, 
				 Bool_t MCcase=kFALSE, 
				 Bool_t PbPbcase=kTRUE, 
				 Bool_t GenerateSignal=kFALSE, 
				 Int_t CentBinLowLimit=0, 
				 Int_t CentBinHighLimit=1,
				 Int_t RMax=11,
				 UInt_t FilterBit=7,
				 Float_t MaxChi2NDF=10,
				 Int_t MinTPCncls=0,
				 Float_t MinSepPairEta=0.02,
				 Float_t MinSepPairPhi=0.045,
				 Float_t SigmaCutTPC=2.0,
				 Float_t SigmaCutTOF=2.0,
				 Int_t NumKt3bins=1,
				 Bool_t V0Mbinning=kFALSE,
				 Int_t TriggerType=0,
				 TString StKName="alien:///alice/cern.ch/user/d/dgangadh/KFile_TPR.root"
				 ) {
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskThreePionRadii", "No analysis manager to connect to.");
    return NULL;
  }
 
  
  //____________________________________________//
  // Create task
  AliThreePionRadii *ThreePionRadiiTask = new AliThreePionRadii("ThreePionRadiiTask");
  if(!ThreePionRadiiTask) return NULL;
  ThreePionRadiiTask->SetLEGOCase(LEGO);
  ThreePionRadiiTask->SetMCdecision(MCcase);
  ThreePionRadiiTask->SetPbPbCase(PbPbcase);
  ThreePionRadiiTask->SetGenerateSignal(GenerateSignal);
  ThreePionRadiiTask->SetCentBinRange(CentBinLowLimit, CentBinHighLimit);
  ThreePionRadiiTask->SetRMax(RMax);
  ThreePionRadiiTask->SetFilterBit(FilterBit);
  ThreePionRadiiTask->SetMaxChi2NDF(MaxChi2NDF);
  ThreePionRadiiTask->SetMinTPCncls(MinTPCncls);
  ThreePionRadiiTask->SetPairSeparationCutEta(MinSepPairEta);
  ThreePionRadiiTask->SetPairSeparationCutPhi(MinSepPairPhi);
  ThreePionRadiiTask->SetNsigmaTPC(SigmaCutTPC);
  ThreePionRadiiTask->SetNsigmaTOF(SigmaCutTOF);
  ThreePionRadiiTask->SetNumKt3Bins(NumKt3bins);
  ThreePionRadiiTask->SetV0Mbinning(V0Mbinning);
  ThreePionRadiiTask->SetTriggerType(TriggerType);
  mgr->AddTask(ThreePionRadiiTask);


  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCF.outputThreePionRadiiAnalysis.root";
  AliAnalysisDataContainer *coutThreePionRadii = mgr->CreateContainer("ThreePionRadiiOutput", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(ThreePionRadiiTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(ThreePionRadiiTask, 1, coutThreePionRadii);
  

  TFile *inputFileFSI = 0;

  ////////////////////////////////////////////////////
  // FSI File
  inputFileFSI = TFile::Open(StKName,"OLD");
  if (!inputFileFSI){
    cout << "Requested file:" << inputFileFSI << " was not opened. ABORT." << endl;
    return NULL;
  }  
  TH1D *FSI2SS[10];
  TH1D *FSI2OS[10];

  for(Int_t FSIindex=0; FSIindex<10; FSIindex++) {
    TString *nameSS=new TString("K2ss_");
    *nameSS += FSIindex;
    FSI2SS[FSIindex] = (TH1D*)inputFileFSI->Get(nameSS->Data());
    TString *nameOS=new TString("K2os_");
    *nameOS += FSIindex;
    FSI2OS[FSIindex] = (TH1D*)inputFileFSI->Get(nameOS->Data());
    //
    FSI2SS[FSIindex]->SetDirectory(0);
    FSI2OS[FSIindex]->SetDirectory(0);
  }
  
  ThreePionRadiiTask->SetFSICorrelations( kTRUE, FSI2SS , FSI2OS);
  ////////////////////////////////////////////////////
  
  // Return the task pointer
  return ThreePionRadiiTask;
}
