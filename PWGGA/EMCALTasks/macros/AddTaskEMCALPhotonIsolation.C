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
                                                                 const char*            nclusters                 = "EmcalClusters",
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
                                                                 )
{
  
  Printf("Preparing neutral cluster analysis\n");
  
    // #### Define manager and data container names
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    ::Error("AddTaskEMCALPhotonIsolation", "No analysis manager to connect to.");
    return NULL;
  }
  
  
  printf("Creating container names for cluster analysis\n");
  TString myContName("");
  if(bIsMC)
    myContName = Form("Analysis_Neutrals_MC");
  else
    myContName = Form("Analysis_Neutrals");
  
  myContName.Append(Form("_TM_%s_CPVe%.2lf_CPVp%.2lf_IsoMet%d_EtIsoMet%d_UEMet%d_TPCbound_%s_IsoConeR%.1f_NLMCut_%s_nNLM%d_SSsmear_%s_Width%.3f_Mean_%.3f",bTMClusterRejection? "On" :"Off", TMdeta , TMdphi ,iIsoMethod,iEtIsoMethod,iUEMethod,bUseofTPC ? "Yes" : "No",iIsoConeRadius,bNLMCut ? "On": "Off",NLMCut, iSmearingSS ? "On":"Off",iWidthSSsmear,iMean_SSsmear));
  
    // #### Define analysis task
  AliAnalysisTaskEMCALPhotonIsolation* task = new AliAnalysisTaskEMCALPhotonIsolation("Analysis",bHisto);
  
    // #### Task preferences
  task->SetOutputFormat(iOutput);
  task->SetLCAnalysis(kFALSE);
  task->SetIsoConeRadius(iIsoConeRadius);
  task->SetEtIsoThreshold(EtIso); // after should be replace by EtIso
  task->SetCTMdeltaEta(TMdeta); // after should be replaced by TMdeta
  task->SetCTMdeltaPhi(TMdphi); // after should be replaced by TMdphi
  task->SetQA(kTRUE);
  task->SetIsoMethod(iIsoMethod);
  task->SetEtIsoMethod(iEtIsoMethod);
  task->SetUEMethod(iUEMethod);
  task->SetUSEofTPC(bUseofTPC);
  task->SetMC(bIsMC);
  task->SetM02Smearing(iSmearingSS);
  task->SetWidth4Smear(iWidthSSsmear);
  task->SetMean4Smear(iMean_SSsmear);
  
  if(bIsMC && bMCNormalization) task->SetIsPythia(kTRUE);
  
  task->SetNLMCut(bNLMCut,NLMCut);
  
  
  
  TString name(Form("PhotonIsolation_%s_%s", ntracks, nclusters));
  cout<<"name of the containers  "<<name.Data()<<endl;
  
    // tracks to be used for the track matching (already used in TM task, TPC only tracks)
  AliTrackContainer *trackCont  = task->AddTrackContainer("tracks");
  if(!trackCont)
    Printf("Error with TPCOnly!!");
  trackCont->SetName("tpconlyMatch");
  trackCont->SetTrackFilterType(AliEmcalTrackSelection::kTPCOnlyTracks);
    // clusters to be used in the analysis already filtered
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  
    // tracks to be used in the analysis (Hybrid tracks)
  AliTrackContainer * tracksForAnalysis = task->AddTrackContainer("tracks");
  if(!tracksForAnalysis)
    Printf("Error with Hybrids!!");
  tracksForAnalysis->SetName("filterTracksAna");
  tracksForAnalysis->SetFilterHybridTracks(kTRUE);
  tracksForAnalysis->SetTrackCutsPeriod(periodstr);
  tracksForAnalysis->SetDefTrackCutsPeriod(periodstr);
  
  
  Printf("Name of Tracks for matching: %s \n Name for Tracks for Isolation: %s",trackCont->GetName(),tracksForAnalysis->GetName());
  
  printf("Task for neutral cluster analysis created and configured, pass it to AnalysisManager\n");
    // #### Add analysis task
  manager->AddTask(task);
  
  
  AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:NeutralClusters",AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *cinput  = manager->GetCommonInputContainer();
  manager->ConnectInput(task, 0, cinput);
  manager->ConnectOutput(task, 1, contHistos);
  
  
  
  return task;
}
