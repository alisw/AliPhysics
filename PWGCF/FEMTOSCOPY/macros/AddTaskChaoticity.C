AliChaoticity *AddTaskChaoticity(TString ParListName = "alien:///alice/cern.ch/user/d/dgangadh/ParListLego_def.txt") {
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }
  
  Char_t dummy[100];
  Bool_t lego;
  Bool_t MCcase;
  Bool_t PbPbcase;
  Bool_t GenSignal;
  Bool_t TabulatePairs;
  Int_t CentBinLow;
  Int_t CentBinHigh;
  Int_t RBinMax;
  Int_t FixedLambdaBin;
  Char_t WeightName[500];
  Char_t MomResName[500];
  Char_t KName[500];
  UInt_t FilterBit;
  Float_t MinPairSep;
  Float_t NsigmaTPC;
  Float_t NsigmaTOF;
  
  ifstream ParList;
  ParList.open(ParListName.Data());
  
  ParList >> dummy >> lego;
  ParList >> dummy >> MCcase;
  ParList >> dummy >> PbPbcase;
  ParList >> dummy >> GenSignal;
  ParList >> dummy >> TabulatePairs;
  ParList >> dummy >> CentBinLow;
  ParList >> dummy >> CentBinHigh;
  ParList >> dummy >> RBinMax;
  ParList >> dummy >> FixedLambdaBin;
  ParList >> dummy >> WeightName;
  ParList >> dummy >> MomResName;
  ParList >> dummy >> KName;
  ParList >> dummy >> FilterBit;
  ParList >> dummy >> MinPairSep;
  ParList >> dummy >> NsigmaTPC;
  ParList >> dummy >> NsigmaTOF;
  //
  ParList.close();

  TString StWeightName(WeightName);
  TString StMomResName(MomResName);
  TString StKName(KName);

  //____________________________________________//
  // Create task
  AliChaoticity *ChaoticityTask = new AliChaoticity("ChaoticityTask");
  if(!ChaoticityTask) return NULL;
  ChaoticityTask->SetLEGOCase(lego);
  ChaoticityTask->SetMCdecision(MCcase);
  ChaoticityTask->SetPbPbCase(PbPbcase);
  ChaoticityTask->SetGenerateSignal(GenSignal);
  ChaoticityTask->SetTabulatePairs(TabulatePairs);
  ChaoticityTask->SetCentBinRange(CentBinLow, CentBinHigh);
  ChaoticityTask->SetRBinMax(RBinMax);
  ChaoticityTask->SetFixedLambdaBin(FixedLambdaBin);
  ChaoticityTask->SetFilterBit(FilterBit);
  ChaoticityTask->SetPairSeparationCut(MinPairSep);
  ChaoticityTask->SetNsigmaTPC(NsigmaTPC);
  ChaoticityTask->SetNsigmaTOF(NsigmaTOF);
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
