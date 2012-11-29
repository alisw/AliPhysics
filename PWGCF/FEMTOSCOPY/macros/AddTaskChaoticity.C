AliChaoticity *AddTaskChaoticity(bool MCcase=kFALSE, bool Tabulatecase=kFALSE, bool PbPbcase=kTRUE, int CentLow=0, int CentHigh=1, TString inputFileNameWeight = "alien:///alice/cern.ch/user/d/dgangadh/WeightFile.root", TString inputFileNameMomRes = "alien:///alice/cern.ch/user/d/dgangadh/MomResFile.root", TString inputFileNameFSI = "alien:///alice/cern.ch/user/d/dgangadh/KFile.root") {
 
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }
  
  //____________________________________________//
  // Create tasks
  AliChaoticity *ChaoticityTask = new AliChaoticity("ChaoticityTask", MCcase, Tabulatecase, PbPbcase, CentLow, CentHigh, kTRUE);
  if(!ChaoticityTask) exit(-1);
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
    ////////////////////////////////////////////////////
    // C2 Weight File
    Int_t ktbins = ChaoticityTask->GetNumKtbins();
    Int_t cbins = ChaoticityTask->GetNumCentbins();
    TH3F *weightHisto[ktbins][cbins] = 0;
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

    //
    ////////////////////////////////////////////////////
    // Momentum Resolution File
    TH2D *momResHisto2D = 0;
    TH3D *momResHisto3D[5] = 0;
    momResHisto2D = (TH2D*)inputFileMomRes->Get("MomResHisto_pp");
    momResHisto3D[0] = (TH3D*)inputFileMomRes->Get("MomResHisto_3d_ppp_term1");
    momResHisto3D[1] = (TH3D*)inputFileMomRes->Get("MomResHisto_3d_ppp_term2");
    momResHisto3D[2] = (TH3D*)inputFileMomRes->Get("MomResHisto_3d_ppp_term3");
    momResHisto3D[3] = (TH3D*)inputFileMomRes->Get("MomResHisto_3d_ppp_term4");
    momResHisto3D[4] = (TH3D*)inputFileMomRes->Get("MomResHisto_3d_ppp_term5");
    ChaoticityTask->SetMomResCorrections( kTRUE, momResHisto2D, momResHisto3D );
    ////////////////////////////////////////////////////
  }// !MCcase and !Tabulatecase
  

  ////////////////////////////////////////////////////
  // FSI File
  inputFileFSI = TFile::Open(inputFileNameFSI,"OLD");
  if (!inputFileFSI){
    cout << "Requested file:" << inputFileFSI << " was not opened. ABORT." << endl;
    return;
  }  
  TH2D *FSI2D[2];
  TH3D *FSI3D[2];
  FSI2D[0] = (TH2D*)inputFileFSI->Get("K2ss");
  FSI2D[1] = (TH2D*)inputFileFSI->Get("K2os");
  FSI3D[0] = (TH3D*)inputFileFSI->Get("K3ss");
  FSI3D[1] = (TH3D*)inputFileFSI->Get("K3os");
  FSI2D[0]->SetDirectory(0);
  FSI2D[1]->SetDirectory(0);
  FSI3D[0]->SetDirectory(0);
  FSI3D[1]->SetDirectory(0);
  ChaoticityTask->SetFSICorrelations( kTRUE, FSI2D , FSI3D);
  ////////////////////////////////////////////////////
   

  // Return the task pointer
  return ChaoticityTask;
}
