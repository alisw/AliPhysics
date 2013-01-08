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

 
  if(!Tabulatecase){
    inputFileWeight = TFile::Open(inputFileNameWeight,"OLD");
    if (!inputFileWeight){
      cout << "Requested file:" << inputFileWeight << " was not opened. ABORT." << endl;
      return;
    }
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
  }// Tabulatecase check
  
  if(!MCcase && !Tabulatecase){
    
    inputFileMomRes = TFile::Open(inputFileNameMomRes,"OLD");
    if (!inputFileMomRes){
      cout << "Requested file:" << inputFileMomRes << " was not opened. ABORT." << endl;
      return;
    }
    ////////////////////////////////////////////////////
    // Momentum Resolution File
    TH2D *momResHisto2D = 0;
    momResHisto2D = (TH2D*)inputFileMomRes->Get("MomResHisto_pp");
    ChaoticityTask->SetMomResCorrections( kTRUE, momResHisto2D);
    ////////////////////////////////////////////////////
  }// MCcase and Tabulatecase check
  

  ////////////////////////////////////////////////////
  // FSI File
  inputFileFSI = TFile::Open(inputFileNameFSI,"OLD");
  if (!inputFileFSI){
    cout << "Requested file:" << inputFileFSI << " was not opened. ABORT." << endl;
    return;
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
