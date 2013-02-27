AliAnalysisTask *AddTaskHFEpPb(Bool_t isAOD = kFALSE){
  // Switches for the TOF-TPC analysis
  Bool_t kTPC_Only                   = kTRUE;
  Bool_t kTPCTOF_Ref                 = kTRUE;
  Bool_t kTPCTOF_Sys                 = kFALSE;
  Bool_t kTRD_Test                   = kFALSE;
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_hfe_HFE", "No analysis manager found.");
    return 0;
  }
  
  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }
  
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  
  //@@ 0 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  // Default settings 
  const int kDefITScl = 3;
  const int kDefTPCcl = 110;
  const int kDefTPCclPID = 80;
  const double kDefDCAr = 1.;
  const double kDefDCAz = 2.;
  // For pPb, new fits by Martin, 07.02.2013
  const double kDefTPCs = 0.045;
  // it seems the sigma is 1 and not 0.88 (3 sigma = 2.685)
  const double kDefTPCu = 3.045;   
  const double kDefTOFs = 3.;

  Double_t dEdxlm[12] = {0.045,0.045,0.045,0.045,0.045,0.045,0.045,0.045,0.045,0.045,0.045,0.045};  //  50%
  Double_t dEdxhm[12] = {2.685,2.685,2.685,2.685,2.685,2.685,2.685,2.685,2.685,2.685,2.685,2.685};
  // For systematics:
  Double_t tpcl0[12]  = {-0.835,-0.835,-0.835,-0.835,-0.835,-0.835,-0.835,-0.835,-0.835,-0.835,-0.835,-0.835};  //  84%
  Double_t tpcl1[12]  = {-0.184,-0.184,-0.184,-0.184,-0.184,-0.184,-0.184,-0.184,-0.184,-0.184,-0.184,-0.184};  //  60%
  Double_t tpcl2[12]  = {0.265,0.265,0.265,0.265,0.265,0.265,0.265,0.265,0.265,0.265,0.265,0.265};  //   40%

  
  // For TPC only
  const double kDefTPCs1 = 0.034;
  const double kDefTPCu1 = 2.674;
  
  if(kTPC_Only){
    // for the moment (09.02.2013) use the same as TOF-TPC. Refine later on
    RegisterTaskPID2(MCthere,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth);
  }
  
  //------------------------------//
  //   TPC-TOF analysis
  //------------------------------//
  
  if(kTPCTOF_Ref){
    // Reference task
     // Fits by Martin 07.02.2013
    // mean = 0.045
    // sigma = 0.88
    // 50%
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,0.295,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);

    // add centrality
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth); // 50%
    //if (!MCthere){
    //  RegisterTask(MCthere,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl0[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);  // 84%
    //}
  }
  if(kTPCTOF_Sys){

    // TPC PID
    if (!MCthere){
      // 84%
      RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,-0.835,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
      // 60%
      RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,-0.184,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
      // 40%
      RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,0.265,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
    }

    // TOF PID
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,2.0,0,AliHFEextraCuts::kBoth);
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,4.0,0,AliHFEextraCuts::kBoth);
    // TOF latest mismatch - helps nothing
    //RegisterTaskPID2(MCthere,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,1,AliHFEextraCuts::kBoth);

    // ITS hits and SPD request
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,4,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,5,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,3,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kFirst);
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,4,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kFirst);

    // TPC clusters
    RegisterTaskPID2(MCthere,isAOD,100,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTaskPID2(MCthere,isAOD,120,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTaskPID2(MCthere,isAOD,130,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
    // TPC clusters PID
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,100,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);

  }
  
  return NULL;
}


//===============================================================================

AliAnalysisTask *RegisterTask(Bool_t useMC, Bool_t isAOD, Int_t tpcCls=120, Int_t tpcClsPID = 80, 
			      Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
			      Double_t *tpcdEdxcutlow=NULL, 
			      Double_t *tpcdEdxcuthigh=NULL, 
			      //Double_t tpcs=-0.0113, Double_t tpcu=3.09, 
			      Double_t tofs=3., Int_t tofm=0,
			      Int_t itshitpixel = AliHFEextraCuts::kBoth, 
			      Bool_t withetacorrection = kFALSE){

  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t tpclow = 0;
  if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;
  Int_t ietacorr = 0;
  if(withetacorrection) ietacorr = 1;
  
  TString appendix(TString::Format("centTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%d",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,tofm,ipixelany));
  printf("Add macro appendix %s\n", appendix.Data());


  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEpPb.C");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEpPb(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					  tpcdEdxcutlow,tpcdEdxcuthigh,
					  tofs,tofm,itshitpixel);

  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();
  if (useMC)
    task->SetHasMCData(kTRUE);
  else{
    task->SetHasMCData(kFALSE);
  }
  task->SelectCollisionCandidates(AliVEvent::kINT7);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFEpPbcent";
  containerName += appendix.Data();
  printf("container name: %s\n", containerName.Data());
 
  //create data containers
  task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), 
					      TList::Class(), AliAnalysisManager::kOutputContainer, 
					      containerName.Data()));
  task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(), 
					      AliAnalysisManager::kOutputContainer, 
					      containerName.Data()));
  mgr->ConnectInput(task,  0, cinput );

  mgr->AddTask(task);
  return NULL;
}

//===============================================================================

AliAnalysisTask *RegisterTaskPID2(Bool_t useMC, Bool_t isAOD, Int_t tpcCls=120, Int_t tpcClsPID = 80, 
				  Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
				  Double_t tpcs=-0.0113, Double_t tpcu=3.09, Double_t tofs=3., 
				  Int_t tofm=0,
				  Int_t itshitpixel = AliHFEextraCuts::kBoth, 
				  Bool_t withetacorrection = kTRUE){
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEmbpPb.C");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEmbpPb(useMC, isAOD, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					     tpcs,tpcu,tofs,tofm,3.,kFALSE,kTRUE,kFALSE,itshitpixel,withetacorrection);
  if()
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();
  if (useMC)
    task->SetHasMCData(kTRUE);
  else{
    task->SetHasMCData(kFALSE);
  }
  task->SelectCollisionCandidates(AliVEvent::kINT7);
 
  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t itpcs = (Int_t)(tpcs*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;
  Int_t ietacorr = 0;
  if(withetacorrection) ietacorr = 1;
  
  TString appendix(TString::Format("TPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%d",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,itpcs,itofs,tofm,ipixelany));
  printf("Add macro appendix %s\n", appendix.Data());

  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFEpPbmb";
  containerName += appendix.Data();
  printf("container name: %s\n", containerName.Data());

  //create data containers
  task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), 
					      TList::Class(), AliAnalysisManager::kOutputContainer, 
					      containerName.Data()));
  task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(), 
					      AliAnalysisManager::kOutputContainer, 
					      containerName.Data()));
  mgr->ConnectInput(task,  0, cinput );

  mgr->AddTask(task);
  return NULL;
}


//=========================================================================
