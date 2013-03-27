AliChaoticity *AddTaskChaoticity(
				 Bool_t LEGO=kTRUE, 
				 Bool_t MCcase=kFALSE, 
				 Bool_t PbPbcase=kTRUE, 
				 Bool_t GenerateSignal=kFALSE, 
				 Bool_t TabulatePairs=kFALSE, 
				 Int_t CentBinLowLimit=0, 
				 Int_t CentBinHighLimit=1,
				 Int_t RBinMax=5,
				 Int_t FixedLambdaBin=11,
				 UInt_t FilterBit=7,
				 Float_t MinSepPair=0.035,
				 Float_t SigmaCutTPC=2.0,
				 Float_t SigmaCutTOF=2.0,
				 TString StWeightName="alien:///alice/cern.ch/user/d/dgangadh/WeightFile.root",
				 TString StMomResName="alien:///alice/cern.ch/user/d/dgangadh/MomResFile.root",
				 TString StKName="alien:///alice/cern.ch/user/d/dgangadh/KFile.root"
				 ) {
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }
 

  //____________________________________________//
  // Create task
  AliChaoticity *ChaoticityTask = new AliChaoticity("ChaoticityTask");
  if(!ChaoticityTask) return NULL;
  ChaoticityTask->SetLEGOCase(LEGO);
  ChaoticityTask->SetMCdecision(MCcase);
  ChaoticityTask->SetPbPbCase(PbPbcase);
  ChaoticityTask->SetGenerateSignal(GenerateSignal);
  ChaoticityTask->SetTabulatePairs(TabulatePairs);
  ChaoticityTask->SetCentBinRange(CentBinLowLimit, CentBinHighLimit);
  ChaoticityTask->SetRBinMax(RBinMax);
  ChaoticityTask->SetFixedLambdaBin(FixedLambdaBin);
  ChaoticityTask->SetFilterBit(FilterBit);
  ChaoticityTask->SetPairSeparationCut(MinSepPair);
  ChaoticityTask->SetNsigmaTPC(SigmaCutTPC);
  ChaoticityTask->SetNsigmaTOF(SigmaCutTOF);
  mgr->AddTask(ChaoticityTask);


  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCF.outputChaoticityAnalysis.root";
  AliAnalysisDataContainer *coutChaoticity = mgr->CreateContainer("ChaoticityOutput", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(ChaoticityTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(ChaoticityTask, 1, coutChaoticity);
  

  TFile *inputFileWeight = 0;
  TFile *inputFileMomRes = 0;
  TFile *inputFileFSI = 0;

 
  if(!TabulatePairs){
    inputFileWeight = TFile::Open(StWeightName,"OLD");
    if (!inputFileWeight){
      cout << "Requested file:" << inputFileWeight << " was not opened. ABORT." << endl;
      return NULL;
    }
    ////////////////////////////////////////////////////
    // C2 Weight File
    const Int_t ktbins_temp = ChaoticityTask->GetNumKtBins();
    const Int_t cbins_temp = ChaoticityTask->GetNumCentBins(); 
    const Int_t ktbins = ktbins_temp;
    const Int_t cbins = cbins_temp;

    TH3F *weightHisto[ktbins][cbins];
    for(Int_t i=0; i<ktbins; i++){
      for(Int_t j=0; j<cbins; j++){
	TString name = "Weight_Kt_";
	name += i;
	name += "_Ky_0_M_";
	name += j;
	name += "_ED_0";
	
	weightHisto[i][j] = (TH3F*)inputFileWeight->Get(name);
      }
    }
    ChaoticityTask->SetWeightArrays( kTRUE, weightHisto );
    ////////////////////////////////////////////////////
  }// TabulatePairs check
  
  if(!MCcase && !TabulatePairs){
    
    inputFileMomRes = TFile::Open(StMomResName,"OLD");
    if (!inputFileMomRes){
      cout << "Requested file:" << inputFileMomRes << " was not opened. ABORT." << endl;
      return NULL;
    }
    ////////////////////////////////////////////////////
    // Momentum Resolution File
    TH2D *momResHisto2D = 0;
    momResHisto2D = (TH2D*)inputFileMomRes->Get("MomResHisto_pp");
    ChaoticityTask->SetMomResCorrections( kTRUE, momResHisto2D);
    ////////////////////////////////////////////////////
  }// MCcase and TabulatePairs check
  

  ////////////////////////////////////////////////////
  // FSI File
  inputFileFSI = TFile::Open(StKName,"OLD");
  if (!inputFileFSI){
    cout << "Requested file:" << inputFileFSI << " was not opened. ABORT." << endl;
    return NULL;
  }  
  TH2D *FSI2gaus[2];
  TH2D *FSI2therm[2];
  TH3D *FSI3ss[6];
  TH3D *FSI3os[6];
  FSI2gaus[0] = (TH2D*)inputFileFSI->Get("K2ssG");
  FSI2gaus[1] = (TH2D*)inputFileFSI->Get("K2osG");
  FSI2therm[0] = (TH2D*)inputFileFSI->Get("K2ssT");
  FSI2therm[1] = (TH2D*)inputFileFSI->Get("K2osT");
  for(Int_t CB=0; CB<6; CB++) {
    TString *nameSS=new TString("K3ss_");
    *nameSS += CB;
    FSI3ss[CB] = (TH3D*)inputFileFSI->Get(nameSS->Data());
    TString *nameOS=new TString("K3os_");
    *nameOS += CB;
    FSI3os[CB] = (TH3D*)inputFileFSI->Get(nameOS->Data());
  }
  //
  FSI2gaus[0]->SetDirectory(0);
  FSI2gaus[1]->SetDirectory(0);
  FSI2therm[0]->SetDirectory(0);
  FSI2therm[1]->SetDirectory(0);
  for(Int_t CB=0; CB<6; CB++) {
    FSI3ss[CB]->SetDirectory(0);
    FSI3os[CB]->SetDirectory(0);
  }
  ChaoticityTask->SetFSICorrelations( kTRUE, FSI2gaus, FSI2therm , FSI3os, FSI3ss);
  ////////////////////////////////////////////////////
  
  
  // Return the task pointer
  return ChaoticityTask;
}
