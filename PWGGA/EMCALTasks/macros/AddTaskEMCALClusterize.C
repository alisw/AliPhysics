AliAnalysisTaskEMCALClusterize* AddTaskEMCALClusterize(
						       TString & arrayName,
                                                       const Int_t   bMC        = kFALSE,
                                                       const Bool_t  exotic     = kTRUE,
                                                       const TString name       = "V1Unfold", 
                                                       const TString trigger    = "", 
                                                       const Bool_t  tm         = kTRUE, 
                                                       const Int_t   minEcell   = 50,
                                                       const Int_t   minEseed   = 100,
                                                       const Int_t   maxDeltaT  = 250,
                                                       const Int_t   timeWindow = 1000,
                                                       const Int_t   minEUnf    = 15, 
                                                       const Int_t   minFrac    = 1,
                                                       const Bool_t  bRecalE    = kTRUE,
                                                       const Bool_t  bBad       = kTRUE,
                                                       const Bool_t  bRecalT    = kTRUE,
                                                       const Bool_t  bNonLine   = kFALSE,
                                                       const Int_t   nRowDiff   = 1,
                                                       const Int_t   nColDiff   = 1
                                                       )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEMCALClusterize", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEMCALClusterize", "This clusterize requires an input event handler");
    return NULL;
  }
  
  printf("Passed Settings : mc %d, exo %d, name %s, trigger %s, tm %d\n",bMC,exotic,name.Data(),trigger.Data(),tm);
  printf("                  Ecell %d, Eseed %d, dT %d, wT %d, minUnf %d, minFrac %d \n",minEcell, minEseed,maxDeltaT,timeWindow,minEUnf,minFrac);
  printf("                  recalE %d, bad %d, recalT %d, nonlin %d, rowDiff %d, colDiff %d \n",bRecalE,bBad,bRecalT,bNonLine,nRowDiff,nColDiff);

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliAnalysisTaskEMCALClusterize* clusterize = new AliAnalysisTaskEMCALClusterize(Form("EMCALClusterize%s_Ecell%d_Eseed%d_DT%d_WT%d",
                                                                                       name.Data(),minEcell,minEseed,maxDeltaT,timeWindow));


  //clusterize->SetOCDBPath("local://$ALICE_ROOT/OCDB");

  // Some general settings to create AOD file in case we want to keep it
  clusterize->SwitchOffFillAODCaloCells();
  clusterize->SwitchOffFillAODHeader();
  clusterize->FillAODFile(kFALSE); // fill aod.root with clusters?, not really needed for analysis.

  // Do track matching after clusterization
  if(tm) clusterize->SwitchOnTrackMatching();
  else   clusterize->SwitchOffTrackMatching();
  
  //-------------------------------------------------------
  // Set clusterization parameters via rec param
  //-------------------------------------------------------

  AliEMCALRecParam * params = clusterize->GetRecParam();

  // Position and SS weight parameter
  params->SetW0(4.5);

  // Time cuts, depend on data type (no cells time in AODs)
  TString sHandler((mgr->GetInputEventHandler())->ClassName());
  if(sHandler.Contains("AOD"))
  {
    printf("AliAnalysisTaskEMCALClusterize - Open time cuts for AODs\n");
    params->SetTimeCut(1e6);//Open this cut for AODs
    params->SetTimeMin(-1); //Open this cut for AODs
    params->SetTimeMax(1e6);//Open this cut for AODs    
  }
  else
  {
    printf("AliAnalysisTaskEMCALClusterize - Set time cuts for ESDs\n");
    if(maxDeltaT > 1) params->SetTimeCut(maxDeltaT*1.e-9);
    else            { params->SetTimeCut(250*1.e-9); printf("default maxDeltaT = 250 ns\n"); }// Same as in reco
    
    if(timeWindow > 1)
    {
      params->SetTimeMin(-1*timeWindow*1.e-9);
      params->SetTimeMax(timeWindow*1.e-9);
    }
    else
    { 
      if(bRecalT && !bMC)
      {
        params->SetTimeMin(-250*1.e-9);
        params->SetTimeMax( 250*1.e-9);
        printf("default time window for calibrated time -250 ns < T < 250 ns\n");
      }
      else 
      {
        // same as in reco, USE IF NO TIME RECALIBRATION
        params->SetTimeMin(425*1.e-9);
        params->SetTimeMax(825*1.e-9);
        printf("default time window 425 ns < T < 825 ns\n");
      }
    }
  }

  // Energy cuts
  params->SetClusteringThreshold(minEseed/1.e3);                                          
  params->SetMinECut            (minEcell/1.e3); 

  // Clusterizer type
  if(name.Contains("V2")) params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2);
  if(name.Contains("V1")) params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv1);
  if(name.Contains("NxN"))
  {
    params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerNxN);
    printf("Set NxN cluster size to %dx%d (row diff %d, col diff %d)\n",2*nRowDiff+1,2*nColDiff+1,nRowDiff,nColDiff);
    params->SetNxM(nRowDiff, nColDiff);
  }

  //-------------------------------------------------------
  // Unfolding, 2 options :
  //-------------------------------------------------------

  //    1) Just unfold existing clusters
  if(name.Contains("JustUnfold"))
    clusterize->JustUnfold(kTRUE); // if TRUE, do just unfolding, do not recluster cells
  else  
    clusterize->JustUnfold(kFALSE); 

  //   2) Unfold clusters created in the clusterize (revise settings)
  if (name.Contains("Unfold"))
  {
    clusterize->SwitchOnCellEnergySelection();
    clusterize->SetCellCuts(minEUnf/1000.,minFrac/10000.);
    printf("AliAnalysisTaskEMCALClusterize - Cuts: min E %f, frac %f\n",minEUnf/1000.,minFrac/10000.);
    //clusterize->SwitchOffCellEnergySelection(); 
    
    if(!name.Contains("Just"))
      params->SetUnfold(kTRUE);
    else 
      params->SetUnfold(kFALSE);
    
  } // unfold
  
  //-------------------------------------------------------
  // Configure AliEMCALRecoUtils
  //-------------------------------------------------------

  AliEMCALRecoUtils * reco = clusterize->GetRecoUtils();
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigureEMCALRecoUtils.C"); // $ALICE_ROOT/PWGGA/EMCALTasks/macros
  ConfigureEMCALRecoUtils(reco,bMC,exotic,bNonLine,bRecalE,bBad,bRecalT);
  
  //-------------------------------------------------------
  //Alignment matrices
  //-------------------------------------------------------

  clusterize->SetImportGeometryFromFile(kTRUE); // import geometry.root file
  
  if(!bMC)
  {    
    clusterize->SwitchOnLoadOwnGeometryMatrices();
  }
  
  //-------------------------------------------------------
  // Trigger options
  //-------------------------------------------------------

  if(trigger=="EMC7")
  {
    printf("Clusterizing task trigger EMC7\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kEMC7);
  }
  else if (trigger=="INT7")
  {
    printf("Clusterizing task trigger INT7\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  else if(trigger=="EMC1")
  {
    printf("Clusterizing task trigger EMC1\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kEMC1);
  }
  else if(trigger=="MB")
  {
    printf("Clusterizing task trigger MB\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kMB);
  }  
  else if(trigger=="PHOS")
  {
    printf("Clusterizing task trigger PHOS\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kPHI7);
  }  
  else if(trigger=="PHOSPb")
  {
    printf("Clusterizing task trigger PHOSPb\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kPHOSPb);
  }
  else if(trigger=="AnyINT")
  {
    printf("Clusterizing task trigger AnyINT\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kAnyINT);
  }  
  else if(trigger=="INT")
  {
    printf("Clusterizing task trigger AnyINT\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kAny);
  }
  else if(trigger=="EMCEGA")
  {
    printf("Clusterizing task trigger EMC Gamma\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  } 
  else if(trigger=="EMCEJE")
  {
    printf("Clusterizing task trigger EMC Jet\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kEMCEJE);
  }
  else if(trigger=="Central")
  {
    printf("Clusterizing task trigger Central\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kCentral);
  } 
  else if(trigger=="SemiCentral")
  {
    printf("Clusterizing task trigger SemiCentral\n");
    clusterize->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  }
  
  // Set clusters branch name, make sure the analysis after this one reads this name
  
  if(name.Contains("NxN")) arrayName = Form("%dx%d_Ecell%d_Eseed%d_DT%d_WT%d",2*nRowDiff+1,2*nColDiff+1,minEcell,minEseed,maxDeltaT,timeWindow);
  else                     arrayName = Form(   "%s_Ecell%d_Eseed%d_DT%d_WT%d",              name.Data(),minEcell,minEseed,maxDeltaT,timeWindow);
  
  clusterize->SetAODBranchName(arrayName);
  
  printf("Created Branch Name: \n",arrayName.Data());
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(clusterize);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer() ;
  
  mgr->ConnectInput  (clusterize, 0,  cinput1 );
  mgr->ConnectOutput (clusterize, 0, coutput1 );
  
  return clusterize;
  
}
