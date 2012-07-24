AliChaoticity *AddTaskChaoticity(bool MCcase=kFALSE, bool Tabulatecase=kFALSE, bool PbPbcase=kTRUE, int CentLow=0, int CentHigh=1, TString inputFileNameWeight = "alien:///alice/cern.ch/user/d/dgangadh/WeightFile.root", TString inputFileNameMomRes = "alien:///alice/cern.ch/user/d/dgangadh/MomResFile.root", TString inputFileNameCoulomb = "alien:///alice/cern.ch/user/d/dgangadh/cc2_l100_all.txt") {
 
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }
  
  //____________________________________________//
  // Create tasks
  AliChaoticity *ChaoticityTask = new AliChaoticity("ChaoticityTask", MCcase, Tabulatecase, PbPbcase, CentLow, CentHigh);
  if(!XiStarTask) exit(-1);
  mgr->AddTask(ChaoticityTask);


  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCF.outputChaoticityAnalysis.root";
  AliAnalysisDataContainer *coutXiStar = mgr->CreateContainer("ChaoticityOutput", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(ChaoticityTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(ChaoticityTask, 1, coutChaoticity);
  

  TFile *inputFileWeight = 0;
  TFile *inputFileMomRes = 0;
  TFile *inputFileCoulomb = 0;

  if(!MCcase && !Tabulatecase){
   
    TGrid::Connect("alien:");
    inputFileWeight = TFile::Open(inputFileNameWeight,"OLD");
    inputFileMomRes = TFile::Open(inputFileNameMomRes,"OLD");
    if (!inputFileWeight){
      cout << "Requested file:" << inputFileWeight << " was not opened. ABORT." << endl;
      return;
    }
    if (!inputFileMomRes){
      cout << "Requested file:" << inputFileMomRes << " was not opened. ABORT." << endl;
      return;
    }
    
    //

    Int_t ktbins = ChaoticityTask->GetNumKtbins();
    Int_t cbins = ChaoticityTask->GetNumCentbins();
    TH3F *weightHisto[ktbins][cbins] = 0;
    for(int i=0; i<ktbins; i++){
      for(int j=0; j<cbins; j++){
	TString name = "Weight_Kt_";
	name += i;
	name += "_0_Ky_0_M_";
	name += j;
	name += "_ED_0";
	
	weightHisto[i][j] = (TH3F*)inputFileWeight->Get(name);
      }
    }
    ChaoticityTask->SetWeightArraysLEGO( weightHisto );
    
    //
    
    TH2D *momResHisto = 0;
    TString name = "MomResHisto_pp";
    momResHisto = (TH2D*)inputFileMomRes->Get(name);
    ChaoticityTask->SetMomResCorrectionsLEGO( momResHisto );
    
    //

    Int_t lines = ChaoticityTask->GetNumCoulLines();
    Int_t rbins = ChaoticityTask->GetNumRValues();
    Float_t qCoul[lines];
    Float_t coulSS[rbins][lines];
    Float_t coulOS[rbins][lines];
    // Set default values
    for(Int_t i=0; i<lines; i++) {
      qCoul[i]=0;
      for(Int_t j=0; j<rbins; j++) {// radii columns
	coulSS[j][i]=-1;
	coulOS[j][i]=-1;
      }
    }
  
    ifstream mystream(inputFileNameCoulomb);
    for(Int_t j=0; j<rbins; j++) {// radii columns (3-10fm in cc2_l100_all.txt)
      for(Int_t i=0; i<lines; i++) {
	mystream >> qCoul[i];
	//
	mystream >> coulSS[j][i];
	mystream >> coulOS[j][i];
      }
    }
    ChaoticityTask->SetCoulCorrectionsLEGO( qCoul, coulSS, coulOS );
    
    
    
  }// MCcase and Tabulatecase



  

  // Return the task pointer
  return ChaoticityTask;
}
