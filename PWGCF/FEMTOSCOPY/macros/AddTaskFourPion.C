AliFourPion *AddTaskFourPion(
				 Bool_t LEGO=kTRUE, 
				 Bool_t MCcase=kFALSE, 
				 Bool_t PbPbcase=kTRUE, 
				 Bool_t GenerateSignal=kFALSE, 
				 Bool_t TabulatePairs=kFALSE,
				 Bool_t LinearInterpolate=kTRUE,
				 Int_t CentBinLowLimit=0, 
				 Int_t CentBinHighLimit=1,
				 Int_t RMax=11,
				 Float_t fcSq=0.7,
				 UInt_t FilterBit=7,
				 Float_t MaxChi2NDF=10,
				 Int_t MinTPCncls=0,
				 Float_t MinSepPairEta=0.02,
				 Float_t MinSepPairPhi=0.045,
				 Float_t SigmaCutTPC=2.0,
				 Float_t SigmaCutTOF=2.0,
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
  FourPionTask->SetLEGOCase(LEGO);
  FourPionTask->SetMCdecision(MCcase);
  FourPionTask->SetPbPbCase(PbPbcase);
  FourPionTask->SetGenerateSignal(GenerateSignal);
  FourPionTask->SetTabulatePairs(TabulatePairs);
  FourPionTask->SetInterpolationType(LinearInterpolate);
  FourPionTask->SetCentBinRange(CentBinLowLimit, CentBinHighLimit);
  FourPionTask->SetRMax(RMax);
  FourPionTask->SetfcSq(fcSq);
  FourPionTask->SetFilterBit(FilterBit);
  FourPionTask->SetMaxChi2NDF(MaxChi2NDF);
  FourPionTask->SetMinTPCncls(MinTPCncls);
  FourPionTask->SetPairSeparationCutEta(MinSepPairEta);
  FourPionTask->SetPairSeparationCutPhi(MinSepPairPhi);
  FourPionTask->SetNsigmaTPC(SigmaCutTPC);
  FourPionTask->SetNsigmaTOF(SigmaCutTOF);
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
