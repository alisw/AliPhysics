AliFourPion *AddTaskFourPion(
			         Bool_t MCcase=kFALSE,
			         Bool_t TabulatePairs=kFALSE,
				 Int_t CentBinLowLimit=0, 
				 Int_t CentBinHighLimit=1,
				 TString StWeightName="alien:///alice/cern.ch/user/d/dgangadh/WeightFile_FourPion.root",
				 TString StMomResName="alien:///alice/cern.ch/user/d/dgangadh/MomResFile_FourPion.root",
				 TString StKName="alien:///alice/cern.ch/user/d/dgangadh/KFile_FourPion.root",
				 TString StMuonName="alien:///alice/cern.ch/user/d/dgangadh/MuonCorrection_FourPion.root"
			     ) {
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFourPion", "No analysis manager to connect to.");
    return NULL;
  }
 

  //____________________________________________//
  // Create task
  AliFourPion *FourPionTask = new AliFourPion("FourPionTask");
  if(!FourPionTask) return NULL;
  FourPionTask->SetLEGOCase(kTRUE);
  FourPionTask->SetMCdecision(MCcase);
  FourPionTask->SetCollisionType(0);
  FourPionTask->SetGenerateSignal(kFALSE);
  FourPionTask->SetTabulatePairs(TabulatePairs);
  FourPionTask->SetInterpolationType(kFALSE);
  FourPionTask->SetCentBinRange(CentBinLowLimit, CentBinHighLimit);
  FourPionTask->SetRMax(11);
  FourPionTask->SetfcSq(0.7);
  FourPionTask->SetFilterBit(7);
  FourPionTask->SetMaxChi2NDF(10);
  FourPionTask->SetMinTPCncls(0);
  FourPionTask->SetPairSeparationCutEta(0.02);
  FourPionTask->SetPairSeparationCutPhi(0.045);
  FourPionTask->SetNsigmaTPC(2.0);
  FourPionTask->SetNsigmaTOF(2.0);
  //
  FourPionTask->SetMixedChargeCut(kFALSE);
  FourPionTask->SetMinPt(0.16);
  FourPionTask->SetMaxPt(1.0);
  FourPionTask->SetKT3transition(0.3);
  FourPionTask->SetKT4transition(0.3);
  mgr->AddTask(FourPionTask);


  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCF.outputFourPionAnalysis.root";
  AliAnalysisDataContainer *coutFourPion = mgr->CreateContainer("FourPionOutput", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(FourPionTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(FourPionTask, 1, coutFourPion);
  

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
    const Int_t ktbins_temp = FourPionTask->GetNumKtBins();
    const Int_t cbins_temp = FourPionTask->GetNumCentBins(); 
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
    FourPionTask->SetWeightArrays( kTRUE, weightHisto );
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
    TH2D *momResHisto2DSC = 0;
    TH2D *momResHisto2DMC = 0;
    momResHisto2DSC = (TH2D*)inputFileMomRes->Get("MRC_C2_SC");
    momResHisto2DMC = (TH2D*)inputFileMomRes->Get("MRC_C2_MC");
    FourPionTask->SetMomResCorrections( kTRUE, momResHisto2DSC, momResHisto2DMC );
    ////////////////////////////////////////////////////

    // Muon corrections
    inputFileMuon = TFile::Open(StMuonName,"OLD");
    if (!inputFileMuon){
      cout << "Requested file:" << inputFileMuon << " was not opened. ABORT." << endl;
      return NULL;
    }
    TH2D *muonHisto2D = 0;
    muonHisto2D = (TH2D*)inputFileMuon->Get("WeightmuonCorrection");
    FourPionTask->SetMuonCorrections( kTRUE, muonHisto2D);

  }// MCcase and TabulatePairs check
  

  ////////////////////////////////////////////////////
  // FSI File
  inputFileFSI = TFile::Open(StKName,"OLD");
  if (!inputFileFSI){
    cout << "Requested file:" << inputFileFSI << " was not opened. ABORT." << endl;
    return NULL;
  }  
  TH1D *FSIss[12];
  TH1D *FSIos[12];
  for(Int_t index=0; index<12; index++) {
    TString *nameSS=new TString("K2ss_");
    *nameSS += index;
    FSIss[index] = (TH1D*)inputFileFSI->Get(nameSS->Data());
    TString *nameOS=new TString("K2os_");
    *nameOS += index;
    FSIos[index] = (TH1D*)inputFileFSI->Get(nameOS->Data());
    //
    FSIss[index]->SetDirectory(0);
    FSIos[index]->SetDirectory(0);
  }
  //
  FourPionTask->SetFSICorrelations( kTRUE, FSIss, FSIos );
  ////////////////////////////////////////////////////
  
  
  // Return the task pointer
  return FourPionTask;
}
