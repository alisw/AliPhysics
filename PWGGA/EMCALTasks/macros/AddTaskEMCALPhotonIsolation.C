  ///////////////////////////////////////////////////////////////////////////
  ///\file AddTaskEMCALPhotonIsolation.C
  ///\brief Configuration of AliAnalysisTaskEMCALPhotonIsolation
  ///
  /// Version to be used in lego train for testing on pp@7TeV
  ///
  /// \author Lucile Ronflette <lucile.ronflette@cern.ch>, SUBATECH, Nantes
  /// \author Davide Francesco Lodato <davide.francesco.lodato@cern.ch>, Utrecht University
  /// \author Marco Marquard <marco.marquard@cern.ch>, University Frankfurt am Main
  ///////////////////////////////////////////////////////////////////////////

AliAnalysisTaskEMCALPhotonIsolation* AddTaskEMCALPhotonIsolation(
                                                                 const char*            periodstr                 = "LHC11c",
                                                                 const char*            ntracks                   = "EmcalTracks",
                                                                 const char*            nclusters                 = "EmcCaloClusters",
                                                                 const UInt_t           pSel                      = AliVEvent::kEMC7,
                                                                 const TString          dType                     = "ESD",
                                                                 const Bool_t	        bHisto                    = kTRUE,
                                                                 const Int_t	        iOutput	  	          = 0,
                                                                 const Bool_t	        bIsMC  	                  = kFALSE,
                                                                 const Bool_t           bMCNormalization          = kFALSE,
                                                                 const Bool_t           bNLMCut                   = kFALSE,
                                                                 const Int_t            NLMCut                    = 0,
                                                                 const Double_t         minPtCutCluster           = 0.3,
                                                                 const Double_t         EtIso                     = 2.,
                                                                 const Int_t            iIsoMethod                = 1,
                                                                 const Int_t            iEtIsoMethod              = 0,
                                                                 const Int_t            iUEMethod                 = 1,
                                                                 const Bool_t           bUseofTPC                 = kFALSE,
                                                                 const Double_t         TMdeta                    = 0.02,
                                                                 const Double_t         TMdphi                    = 0.03,
                                                                 const Bool_t           bTMClusterRejection       = kTRUE,
                                                                 const Bool_t           bTMClusterRejectionInCone = kTRUE,
                                                                 const Float_t          iIsoConeRadius            = 0.4,
                                                                 const Bool_t           iSmearingSS               = kFALSE,
                                                                 const Float_t          iWidthSSsmear             = 0.,
                                                                 const Float_t          iMean_SSsmear             = 0.,
                                                                 const Bool_t           iExtraIsoCuts             = kFALSE,
                                                                 const Bool_t           i_pPb                     = kFALSE,
                                                                 const Bool_t           isQA                      = kFALSE,
                                                                 TString                configBasePath            = "",
                                                                 const Int_t            bWhichToSmear             = 0,
                                                                 const Int_t            minNLM                    = 1,
                                                                 const Double_t         TMdetaIso                 = 0.02,
                                                                 const Double_t         TMdphiIso                 = 0.03,
                                                                 const Bool_t           bmcTruth                  = kTRUE,
                                                                 const Bool_t           isLCAnalysis              = kFALSE,
								 TString                L1triggerName             = ""
                                                                 )
{
  Printf("Preparing neutral cluster analysis\n");
  
  // #### Define manager and data container names
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if(!manager){
    ::Error("AddTaskEMCALPhotonIsolation", "No analysis manager to connect to.");
    return NULL;
  }

  // Use the input containers naming convention with "usedefault" (already used in several EMCal tasks and new correction framework)
  TString trackName(ntracks);
  TString clusName(nclusters);

  if(trackName == "usedefault"){
    if(dType == "ESD"){
      trackName = "Tracks";
    }
    else if(dType == "AOD"){
      trackName = "tracks";
    }
    else{
      trackName = "";
    }
  }

  if(clusName == "usedefault"){
    if(dType == "ESD"){
      clusName = "CaloClusters";
    }
    else if(dType == "AOD"){
      clusName = "caloClusters";
    }
    else{
      clusName = "";
    }
  }

  // Set the task output container name
  Printf("Creating container name for cluster analysis\n");

  TString myContName("");
  if(bIsMC){
    myContName = Form("Analysis_Neutrals_MC");
  }
  else{
    myContName = Form("Analysis_Neutrals");
  }
  
  if(L1triggerName.Contains("EG1") || L1triggerName.Contains("EGA1")){
    L1triggerName = "EG1";
  }
  else if(L1triggerName.Contains("EG2") || L1triggerName.Contains("EGA2")){
    L1triggerName = "EG2";
  }
  else{
    L1triggerName = "";
  }

  if(L1triggerName != ""){
    myContName.Append(Form("%s_TM_%s_CPVe%.2lf_CPVp%.2lf_IsoMet%d_EtIsoMet%d_UEMet%d_TPCbound_%s_IsoConeR%.1f_NLMCut_%s_minNLM%d_maxNLM%d_SSsmear_%s_Width%.3f_Mean_%.3f_PureIso_%s_WhichSmear_%d_Trigger_%s",isLCAnalysis?"_LC_Yes":"",bTMClusterRejection? "On" :"Off", TMdeta , TMdphi ,iIsoMethod,iEtIsoMethod,iUEMethod,bUseofTPC ? "Yes" : "No",iIsoConeRadius,bNLMCut ? "On": "Off",minNLM, NLMCut, iSmearingSS ? "On":"Off",iWidthSSsmear,iMean_SSsmear,iExtraIsoCuts?"On":"Off",bWhichToSmear,L1triggerName.Data()));
  }
  else{
    myContName.Append(Form("%s_TM_%s_CPVe%.2lf_CPVp%.2lf_IsoMet%d_EtIsoMet%d_UEMet%d_TPCbound_%s_IsoConeR%.1f_NLMCut_%s_minNLM%d_maxNLM%d_SSsmear_%s_Width%.3f_Mean_%.3f_PureIso_%s_WhichSmear_%d",isLCAnalysis?"_LC_Yes":"",bTMClusterRejection? "On" :"Off", TMdeta , TMdphi ,iIsoMethod,iEtIsoMethod,iUEMethod,bUseofTPC ? "Yes" : "No",iIsoConeRadius,bNLMCut ? "On": "Off",minNLM, NLMCut, iSmearingSS ? "On":"Off",iWidthSSsmear,iMean_SSsmear,iExtraIsoCuts?"On":"Off",bWhichToSmear));
  }

  // #### Define analysis task
  AliAnalysisTaskEMCALPhotonIsolation* task = new AliAnalysisTaskEMCALPhotonIsolation("Analysis",bHisto);
  
  TString configFile("config_PhotonIsolation.C"); // Name of config file
  //  if(gSystem->AccessPathName(configFile.Data())){ // Check for exsisting file and delete it
  //    gSystem->Exec(Form("rm %s",configFile.Data()));
  //  }

  if(configBasePath.IsNull()){ // Check if a specific config should be used and copy appropriate file
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
  ifstream configIn; // Open config file for hash calculation
  configIn.open(configFile);
  TString configStr;
  configStr.ReadFile(configIn);
  TString configMD5 = configStr.MD5();
  configMD5.Resize(5); // Short hash value for usable extension
  TString configFileMD5 = configFile;
  TDatime time; // Get timestamp
  Int_t timeStamp = time.GetTime();
  configFileMD5.ReplaceAll(".C",Form("\_%s_%i.C",configMD5.Data(),timeStamp));

  if(gSystem->AccessPathName(configFileMD5.Data())){ // Add additional identifier if file exists
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
  Printf("Path of config file: %s\n",configFilePath.Data());
  
  // #### Task preferences
  task->SetOutputFormat(iOutput);
  task->SetLCAnalysis(isLCAnalysis);
  task->SetIsoConeRadius(iIsoConeRadius);
  task->SetEtIsoThreshold(EtIso);
  task->SetCTMdeltaEta(TMdeta);
  task->SetCTMdeltaPhi(TMdphi);
  task->SetCTMdeltaEtaIso(TMdetaIso);
  task->SetCTMdeltaPhiIso(TMdphiIso);
  task->SetQA(isQA);
  task->SetIsoMethod(iIsoMethod);
  task->SetEtIsoMethod(iEtIsoMethod);
  task->SetUEMethod(iUEMethod);
  task->SetUSEofTPC(bUseofTPC);
  task->SetMC(bIsMC);
  task->SetM02Smearing(iSmearingSS);
  task->SetWidth4Smear(iWidthSSsmear);
  task->SetMean4Smear(iMean_SSsmear);
  task->SetExtraIsoCuts(iExtraIsoCuts);
  task->SetAnalysispPb(i_pPb);
  task->SetNLMCut(bNLMCut,NLMCut,minNLM);
  task->SetPtBinning(ptBin);
  task->SetM02Binning(M02Bin);
  task->SetEtisoBinning(EtisoBin);
  task->SetEtueBinning(EtueBin);
  task->SetEtaBinning(EtaBin);
  task->SetPhiBinning(PhiBin);
  task->SetLabelBinning(LabelBin);
  task->SetPDGBinning(PDGBin);
  task->SetMomPDGBinning(MomPDGBin);
  task->SetClustPDGBinning(ClustPDGBin);
  task->SetDxBinning(DxBin);
  task->SetDzBinning(DzBin);
  task->SetDecayBinning(DecayBin);
  task->SetSmearForClusters(bWhichToSmear);
  task->SetNeedEmcalGeom(kTRUE);
  task->SetMCtruth(bmcTruth);
  
  if(bIsMC && bMCNormalization) task->SetIsPythia(kTRUE);
  
  TString name(Form("PhotonIsolation_%s_%s", trackName.Data(), clusName.Data()));
  cout<<"Name of the container "<<name.Data()<<endl;
  
  // Tracks to be used for the track matching (already used in TM task, TPC only tracks)
  AliTrackContainer *trackCont = task->AddTrackContainer(trackName);
  if(!trackCont) Printf("Error with TPC-Only Tracks!!");
  trackCont->SetName("tpconlyMatch");
  trackCont->SetTrackFilterType(AliEmcalTrackSelection::kTPCOnlyTracks);
  
  // Tracks to be used in the analysis (Hybrid Tracks)
  AliTrackContainer * tracksForAnalysis = task->AddTrackContainer(trackName);
  if(!tracksForAnalysis) Printf("Error with Hybrids Tracks!!");
  tracksForAnalysis->SetName("filterTracksAna");
  tracksForAnalysis->SetFilterHybridTracks(kTRUE);
  if(!bIsMC){
    tracksForAnalysis->SetTrackCutsPeriod(periodstr);
    tracksForAnalysis->SetDefTrackCutsPeriod(periodstr);
  }

  // Clusters to be used in the analysis (already filtered)
  AliClusterContainer *clusterCont = task->AddClusterContainer(clusName);

  Printf("Name of track container for matching: %s \nName of track container for isolation: %s \nName of cluster container: %s",trackCont->GetName(),tracksForAnalysis->GetName(),clusterCont->GetName());
  
  // #### Add analysis task
  Printf("Task for neutral cluster analysis created and configured, pass it to AnalysisManager\n");
  manager->AddTask(task);
  
  AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:NeutralClusters",AliAnalysisManager::GetCommonFileName()));

  AliAnalysisDataContainer *cinput = manager->GetCommonInputContainer();
  manager->ConnectInput(task, 0, cinput);
  manager->ConnectOutput(task, 1, contHistos);
  
  return task;
}
