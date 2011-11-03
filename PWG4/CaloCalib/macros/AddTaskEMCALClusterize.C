AliAnalysisTaskEMCALClusterize* AddTaskEMCALClusterize(
                                                       const Int_t   bMC        = kFALSE,
                                                       const TString name       = "V1Unfold", 
                                                       const TString trigger    = "", 
                                                       const Int_t   run        = 0, 
                                                       const TString pass       = "pass3",
                                                       const Bool_t  tm         = kTRUE, 
                                                       const Int_t   minEcell   = 50,
                                                       const Int_t   minEseed   = 100,
                                                       const Int_t   maxDeltaT  = 250,
                                                       const Int_t   timeWindow = 1000,
                                                       const Int_t   minEUnf    = 15, 
                                                       const Int_t   minFrac    = 1
                                                       )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALClusterize", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALClusterize", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskEMCALClusterize* clusterize = new AliAnalysisTaskEMCALClusterize(Form("EMCALClusterize%s_Ecell%d_Eseed%d_DT%d_WT%d",
                                                                                       name.Data(),minEcell,minEseed,maxDeltaT,timeWindow));

  clusterize->SetAODBranchName(Form("%s_Ecell%d_Eseed%d_DT%d_WT%d",
                                    name.Data(),minEcell,minEseed,maxDeltaT,timeWindow));

  printf("Created Branch Name: %s_Ecell%d_Eseed%d_DT%d_WT%d\n",name.Data(),minEcell,minEseed, maxDeltaT,timeWindow);
  //clusterize->SetOCDBPath("local://$ALICE_ROOT/OCDB");

  clusterize->SwitchOffFillAODCaloCells();
  clusterize->SwitchOffFillAODHeader();
  clusterize->FillAODFile(kFALSE); // fill aod.root with clusters?, not really needed for analysis.

  if(tm) clusterize->SwitchOnTrackMatching();
  else   clusterize->SwitchOffTrackMatching();
  
  // Settings for LHC11a
  if(run > 140000 && run <= 146860) {
    clusterize->SwitchOnLEDEventsRemoval() ;
    printf("Clusterizer: Reject LED events\n");
  }
  else clusterize->SwitchOffLEDEventsRemoval() ;
  
  printf(" ---- Clusterize RUN >%d< ---- \n",run);
  
  if(run > 140000)  clusterize->SetGeometryName("EMCAL_COMPLETEV1");
  else              clusterize->SetGeometryName("EMCAL_FIRSTYEARV1");

  AliEMCALRecParam * params = clusterize->GetRecParam();

  params->SetW0(4.5);

  //printf("**** InputHandler %s ***\n",(mgr->GetInputEventHandler())->ClassName());
  TString sHandler((mgr->GetInputEventHandler())->ClassName());
  if(sHandler.Contains("AOD")){
    printf("AliAnalysisTaskEMCALClusterize - Open time cuts for AODs\n");
    params->SetTimeCut(1e6);//Open this cut for AODs
    params->SetTimeMin(-1); //Open this cut for AODs
    params->SetTimeMax(1e6);//Open this cut for AODs    
  }
  else{
    printf("AliAnalysisTaskEMCALClusterize - Set time cuts for ESDs\n");
    params->SetTimeCut(maxDeltaT*1.e-9);
    params->SetTimeMin(-1*timeWindow*1.e-9);
    params->SetTimeMax(timeWindow*1.e-9);
  }

  params->SetClusteringThreshold(minEseed/1.e3);                                          
  params->SetMinECut            (minEcell/1.e3); 

  // Clusterizer type
  if(name.Contains("V2")) params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2);
  if(name.Contains("V1")) params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv1);
  if(name.Contains("NxN"))params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerNxN);


  //Unfolding
  if(name.Contains("JustUnfold"))
    clusterize->JustUnfold(kTRUE); // if TRUE, do just unfolding, do not recluster cells
  else  
    clusterize->JustUnfold(kFALSE); 

  if (name.Contains("Unfold")){

    clusterize->SwitchOnCellEnergySelection();
    clusterize->SetCellCuts(minEUnf/1000.,minFrac/10000.);
    printf("AliAnalysisTaskEMCALClusterize - Cuts: min E %f, frac %f\n",minEUnf/1000.,minFrac/10000.);
    //clusterize->SwitchOffCellEnergySelection(); 

    if(!name.Contains("Just"))
    params->SetUnfold(kTRUE);
  else 
    params->SetUnfold(kFALSE);

  } // unfold
  
  TGeoHMatrix* matrix[10];
  AliEMCALRecoUtils * reco = clusterize->GetRecoUtils();
  
  gROOT->LoadMacro("ConfigureEMCALRecoUtils.C"); // $ALICE_ROOT/PWG4/CaloCalib/macros
  ConfigureEMCALRecoUtils(
                          reco,
                          bMC, 
                          matrix,
                          "",//AODB path, default
                          run, 
                          pass
                          );
  
  //Alignment matrices
  
  for (Int_t mod=0;mod<10;mod++)
  {
    //((TGeoHMatrix*) mobj->At(mod))->Print();
    clusterize->SetGeometryMatrixInSM(matrix[mod],mod);
  }
  
  clusterize->SwitchOnLoadOwnGeometryMatrices();
    
  if(trigger=="EMC7"){
    printf("AliAnalysisTaskEMCALClusterize ---Data with kEMC7 Trigger and clusterizer %s\n",name.Data());
    clusterize->SelectCollisionCandidates(AliVEvent::kEMC7);
  }
  else if (trigger=="INT7"){
    printf("AliAnalysisTaskEMCALClusterize ---Data with kINT7 Trigger and clusterizer %s\n",name.Data());
    clusterize->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  if(trigger=="EMC1"){
    printf("AliAnalysisTaskEMCALClusterize ---Data with kEMC1 Trigger and clusterizer %s\n",name.Data());
    clusterize->SelectCollisionCandidates(AliVEvent::kEMC1);
  }
  else if(trigger=="MB"){
    printf("AliAnalysisTaskEMCALClusterize ---Data with kMB Trigger and clusterizer %s\n",name.Data());
    clusterize->SelectCollisionCandidates(AliVEvent::kMB);
  }

  mgr->AddTask(clusterize);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer() ;
  
  mgr->ConnectInput  (clusterize, 0,  cinput1 );
  mgr->ConnectOutput (clusterize, 0, coutput1 );
  
  return clusterize;
  
}
