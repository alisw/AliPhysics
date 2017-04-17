  ///////////////////////////////////////////////////////////////////////////
  ///\file AddTaskEMCALTrackIsolation.C
  ///\brief Configuration of AliAnalysisTaskEMCALTrackIsolation
  ///
  /// Version to be used in lego train for testing on pp@7TeV
  ///
  /// \author Davide Francesco Lodato <davide.francesco.lodato@cern.ch>, Utrecht University
  ///////////////////////////////////////////////////////////////////////////

AliAnalysisTaskEMCALTrackIsolation* AddTaskEMCALTrackIsolation(
                                                               const char*     periodstr      = "LHC11c",
                                                               const char*     ntracks        = "EmcalTracks",
                                                               const char*     nclusters      = "EmcCaloClusters",
                                                               const UInt_t    pSel           = AliVEvent::kEMC7,
                                                               const TString   dType          = "ESD",
                                                               const Bool_t	   bHisto         = kTRUE,
                                                               const Bool_t	   bIsMC  	      = kFALSE,
                                                               const Double_t  PtIso          = 2.,
                                                               const Int_t     iIsoMethod     = 0,
                                                               const Int_t     iPtIsoMethod   = 0,
                                                               const Int_t     iUEMethod      = 1,
                                                               const Bool_t    bUseofTPC      = kFALSE,
                                                               const Bool_t    isLTAnalysis   = kTRUE,
                                                               const Float_t   iIsoConeRadius = 0.4,
                                                               const Bool_t    i_pPb          = kFALSE,
                                                               const Bool_t    isQA           = kTRUE,
                                                               TString         configBasePath = "",
                                                               const Bool_t    bmcTruth       = kTRUE
                                                               )
{
  
  Printf("Preparing neutral cluster analysis\n");
  
    // #### Define manager and data container names
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    ::Error("AddTaskEMCALTrackIsolation", "No analysis manager to connect to.");
    return NULL;
  }
  
  
  printf("Creating container names for Tracks analysis\n");
  TString myContName("");
  if(bIsMC)
    myContName = Form("Analysis_Tracks_MC");
  else
    myContName = Form("Analysis_Tracks");
  
  myContName.Append(Form("_IsoMethod%d_PtIsoMet%d_UEMet%d_TPCbound_%s_IsoConeR%.1f_LeadTrack_%s",iIsoMethod,iPtIsoMethod,iUEMethod,bUseofTPC ? "Yes" : "No",iIsoConeRadius,isLTAnalysis?"yes":"no"));
  
  cout<<myContName<<endl;
    // #### Define analysis task
  AliAnalysisTaskEMCALTrackIsolation* task = new AliAnalysisTaskEMCALTrackIsolation("AnalysisTrack",bHisto);
  
  TString configFile("config_TrackIsolation.C"); //Name of config file
  
  cout<<configBasePath<<endl;
  Printf("original BasePath for config file: %s",configBasePath.Data());
  
  if(configBasePath.IsNull()){ //Check if a specific config should be used and copy appropriate file
    
    configBasePath="$ALICE_PHYSICS/PWGGA/EMCALTasks/macros";
    gSystem->Exec(Form("cp %s/%s .",configBasePath.Data(),configFile.Data()));
  }
  else if(configBasePath.Contains("alien:///")){
    gSystem->Exec(Form("alien_cp %s/%s .",configBasePath.Data(),configFile.Data()));
  }
  else{
    gSystem->Exec(Form("cp %s/%s .",configBasePath.Data(),configFile.Data()));
  }
  
  configBasePath=Form("%s/",gSystem->pwd());
  ifstream configIn; //Open config file for hash calculation
  configIn.open(configFile);
  TString configStr;
  configStr.ReadFile(configIn);
  TString configMD5 = configStr.MD5();
  configMD5.Resize(5); //Short hash value for usable extension
  TString configFileMD5 = configFile;
  TDatime time; //Get timestamp
  Int_t timeStamp = time.GetTime();
  configFileMD5.ReplaceAll(".C",Form("\_%s_%i.C",configMD5.Data(),timeStamp));
  
  if(gSystem->AccessPathName(configFileMD5.Data())){ //Add additional identifier if file exists
    gSystem->Exec(Form("mv %s %s",configFile.Data(),configFileMD5.Data()));
  }
  else{
    while(!gSystem->AccessPathName(configFileMD5.Data())){
      configFileMD5.ReplaceAll(".C","_1.C");
    }
    gSystem->Exec(Form("mv %s %s",configFile.Data(),configFileMD5.Data()));
  }
  
  TString configFilePath(configBasePath+"/"+configFileMD5);
  gROOT->LoadMacro(configFilePath.Data());
  printf("Path of config file: %s\n",configFilePath.Data());
  
  cout<<task<<endl;
    // #### Task preferences
    //  task->SetLCAnalysis(kFALSE);
  cout<<iIsoConeRadius<<"\t"<<PtIso<<"\t"<<iIsoMethod<<"\t"<<iPtIsoMethod<<"\t"<<iUEMethod<<"\t"<<bUseofTPC<<"\t"<<isLTAnalysis<<"\t"<<bIsMC<<"\t"<<i_pPb<<"\t"<<bmcTruth<<endl;
  task->SetIsoConeRadius(iIsoConeRadius);
  task->SetPtIsoThreshold(PtIso); // after should be replace by EtIso
  task->SetIsoMethod(iIsoMethod);
  task->SetPtIsoMethod(iPtIsoMethod);
  task->SetUEMethod(iUEMethod);
  task->SetUSEofTPC(bUseofTPC);
  task->SetLTAnalysis(isLTAnalysis);
  task->SetMC(bIsMC);
  task->SetAnalysispPb(i_pPb);
  task->SetQA(isQA);
  task->SetPtBinning(ptBin);
  task->SetPtisoBinning(PtisoBin);
  task->SetPtueBinning(PtueBin);
  task->SetEtaBinning(EtaBin);
  task->SetPhiBinning(PhiBin);
  task->SetLabelBinning(LabelBin);
  task->SetPDGBinning(PDGBin);
  task->SetMomPDGBinning(MomPDGBin);
  task->SetTrackPDGBinning(TrackPDGBin);
  task->SetDxBinning(DxBin);
  task->SetDzBinning(DzBin);
  task->SetDecayBinning(DecayBin);
  task->SetMCtruth(bmcTruth);
  
  TString name(Form("TrackIsolation_%s", ntracks));
  cout<<"name of the containers  "<<name.Data()<<endl;
  if(iIsoMethod==1)
    AliClusterContainer *clustForAnalysis = task->AddClusterContainer(nclusters);
  
    // tracks to be used in the analysis (Hybrid tracks)
  AliTrackContainer * tracksForAnalysis = task->AddTrackContainer("tracks");
  if(!tracksForAnalysis)
    Printf("Error with Hybrids!!");
  tracksForAnalysis->SetName("filterTracksAna");
  tracksForAnalysis->SetFilterHybridTracks(kTRUE);
  if(!bIsMC){
    tracksForAnalysis->SetTrackCutsPeriod(periodstr);
    tracksForAnalysis->SetDefTrackCutsPeriod(periodstr);
  }
  Printf("Name for Tracks for Isolation: %s",tracksForAnalysis->GetName());
  
  printf("Task for isolated track analysis created and configured, pass it to AnalysisManager\n");
    // #### Add analysis task
  manager->AddTask(task);
  
  
  AliAnalysisDataContainer *cinput  = manager->GetCommonInputContainer();
  
  AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:IsolatedTracks",AliAnalysisManager::GetCommonFileName()));
  manager->ConnectInput(task, 0, cinput);
  manager->ConnectOutput(task, 1, contHistos);
  return task;
}
