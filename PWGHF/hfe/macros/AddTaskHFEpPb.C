AliAnalysisTask *AddTaskHFEpPb(Bool_t isMC = kFALSE,
                               Bool_t isAOD = kFALSE,
                               Bool_t kTPCTOF_Ref = kTRUE,
                               Bool_t kTPC_Only = kFALSE,
                               Bool_t kTPCTOF_Cent = kFALSE,
                               Bool_t kTPCTOF_Sys = kFALSE,
			       Bool_t kTPCTOFTRD_Ref = kFALSE,
			       Bool_t kTPCTOFTRD_mbPID = kFALSE,
                               Bool_t kTPCTOFTRD_PID = kFALSE,
                               int TRDtrigger = 0
  ){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEpPb", "No analysis manager found.");
    return 0;
  }
  
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  
  //@@ 0 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  // Default settings 
  const int kDefITScl = 4;
  const int kDefTPCcl = 110;
  const int kDefTPCclPID = 80;
  const double kDefDCAr = 1.;
  const double kDefDCAz = 2.;
  
  // 24.02.2013 After the latest fits, mean is 0 and sigma is 1
  const double kDefTPCs = 0;
  const double kDefTPCu = 3;   
  const double kDefTOFs = 3.;
  //const int TRDtrigger = 1; // trd trigger type
  if(!kTPCTOFTRD_Ref) TRDtrigger = 0;

  // NEW SPLINES + CORRECTIONS: mean 0, width 1
  // mean is actually 0.06 (abs(eta)<0.6)
  Double_t dEdxhm[12] = {3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06};  // Above 3 sigma we neglect 0.13%

  Double_t tpcl0[12]  = {-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94};  // 84%
  Double_t tpcl1[12]  = {-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44};  // 69%
  Double_t tpcl2[12]  = {0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06};  // 50%

  Double_t tpcl3[12]  = {-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69};  // 77.2%
  Double_t tpcl4[12]  = {-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19};  // 59.7%
  Double_t tpcl5[12]  = {0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186};  // 45%
  Double_t tpcl6[12]  = {0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31};  // 40%
  Double_t tpcl7[12]  = {0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56};  // 30.7%

  // For TPC only
  const double kDefTPCs1 = 0.0;
  const double kDefTPCu1 = 3.0;

  // eta cut
  Int_t etacut=0; // eta cut off
  
  if(kTPC_Only){
    // Reference, 50% TPC PID, with centrality V0A
    printf("Adding TPC-only task\n");
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1);
  }

  //------------------------------//
  //   TPC-TOF analysis
  //------------------------------//

  if(kTPCTOF_Ref){
    // Reference task
    // with centrality V0A
    printf("Adding TPC-TOF task\n");
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1);
  }

  if (kTPCTOF_Cent){
    // For the moment we set V0A as the reference centrality estimator, and test the other ones:
    // 1: V0A, 2: V0M, 3: CL1, 4: ZNA
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,2); // V0M
    RegisterTask(isMC,isAOD, kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,3); // CL1
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,4); // ZNA

  }

  if(kTPCTOF_Sys){

    // TPC PID
    if (!isMC){
      // 84%
      //RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl0[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
      //RegisterTaskPID2(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,-0.835,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
      // 70%
      RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
      // 60%
      RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl2[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
      //RegisterTaskPID2(isMC,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,-0.184,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
      // 40%
      RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl3[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
      //RegisterTaskPID2(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,0.265,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth);
    }

    // TOF PID
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],2.0,0,AliHFEextraCuts::kBoth);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],4.0,0,AliHFEextraCuts::kBoth);
    // TOF latest mismatch - helps nothing
    //RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,1,AliHFEextraCuts::kBoth);

    // ITS hits and SPD request
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,4,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,5,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,3,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,4,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,5,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst);

    // TPC clusters
    RegisterTask(isMC,isAOD,100,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTask(isMC,isAOD,120,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
    // TPC clusters PID
    RegisterTask(isMC,isAOD,kDefTPCcl,100,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
  }

   if(kTPCTOFTRD_Ref){

      RegisterTaskPID2(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger);
      RegisterTaskPID2(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger);
   }

  if(kTPCTOFTRD_mbPID){
      // without TOF
      RegisterTaskPID2mbTRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,6,4,"TRD,TPC",etacut);
      RegisterTaskPID2mbTRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,6,4,"TPC,TRD",etacut);
      RegisterTaskPID2mbTRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,6,2,"TRD,TPC",etacut);
      RegisterTaskPID2mbTRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,5,4,"TRD,TPC",etacut);
      RegisterTaskPID2mbTRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,5,4,"TPC,TRD",etacut);
      RegisterTaskPID2mbTRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,5,2,"TRD,TPC",etacut);
  }

  if(kTPCTOFTRD_PID){
      // without TOF
      RegisterTaskPID2TRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,0,5,4,"TRD,TPC");
      RegisterTaskPID2TRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,0,5,4,"TPC,TRD");
      RegisterTaskPID2TRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,0,6,4,"TRD,TPC");
      RegisterTaskPID2TRD(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,0,6,4,"TPC,TRD");

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
			      Int_t itshitpixel = AliHFEextraCuts::kBoth, Int_t icent=1, 
            Double_t EtaMin=-0.8, Double_t EtaMax=0.8,
			      Bool_t withetacorrection = kFALSE){

  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t tpclow = 0;
  if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;
  Int_t ietacorr = 0;
  if(withetacorrection) ietacorr = 1;
  Int_t iEtaMin = (Int_t)(EtaMin*10.);
  Int_t iEtaMax = (Int_t)(EtaMax*10.);

  printf("Argument passed to function to determine the centrality estimator %i\n", icent);

  if (icent == 2) TString cesti("V0M");
  else if (icent == 3) TString cesti("CL1");
  else if (icent == 4) TString cesti("ZNA");
  else TString cesti("V0A");

  TString appendix(TString::Format("centTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%dce%s",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,tofm,ipixelany,cesti.Data()));
  printf("Add macro appendix %s\n", appendix.Data());


  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEpPb.C");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEpPb(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					  tpcdEdxcutlow,tpcdEdxcuthigh,
					  tofs,tofm,itshitpixel,icent);

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
  containerName += ":HFEtask";
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
				  Bool_t withetacorrection = kTRUE,
          Int_t TRDtrigger=0){
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEmbpPb.C");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEmbpPb(useMC,isAOD, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					     tpcs,tpcu,tofs,tofm,3.,kFALSE,kTRUE,kFALSE,itshitpixel,withetacorrection,0,TRDtrigger);
  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();
  if (useMC)
    task->SetHasMCData(kTRUE);
  else{
    task->SetHasMCData(kFALSE);
  }

  if(TRDtrigger==0) task->SelectCollisionCandidates(AliVEvent::kINT7);
  else task->SelectCollisionCandidates(AliVEvent::kTRD);
 
  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t itpcs = (Int_t)(tpcs*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;
  Int_t ietacorr = 0;
  if(withetacorrection) ietacorr = 1;
  
  TString appendix(TString::Format("mbTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%dtrdtrg%d",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,itpcs,itofs,tofm,ipixelany,TRDtrigger));
  printf("Add macro appendix %s\n", appendix.Data());

  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFEtask";
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

AliAnalysisTask *RegisterTaskPID2mbTRD(Bool_t useMC, Bool_t isAOD, Int_t tpcCls=120, Int_t tpcClsPID = 80,
				  Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
				  Double_t tpcs=-0.0113, Double_t tpcu=3.09, Double_t tofs=3., 
				  Int_t tofm=0,
				  Int_t itshitpixel = AliHFEextraCuts::kBoth, 
				  Bool_t withetacorrection = kTRUE,
				     Int_t TRDtrigger=0,Int_t trdl=6, Int_t trde=4,TString detector="TRD,TPC",Int_t etacut=0){

    TString detused=detector.Copy();
    TPRegexp(",").Substitute(detused,"","g");
    printf("detectors in use %s %s \n",detector.Data(),detused.Data());

    Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t itpcs = (Int_t)(tpcs*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;
  Int_t ietacorr = 0;
  if(withetacorrection) ietacorr = 1;
  
  TString appendix(TString::Format("mbTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%dtrdtrg%dtrdl%itrde%iPID%setacut%i",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,itpcs,itofs,tofm,ipixelany,TRDtrigger,trdl,trde,detused.Data(),etacut));
  printf("Add macro appendix %s\n", appendix.Data());



  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEmbpPbTRD.C");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEmbpPbTRD(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz,
					     tpcs,tpcu,tofs,tofm,3.,kFALSE,kTRUE,kFALSE,itshitpixel,withetacorrection,0,TRDtrigger,trdl,trde,detector,etacut);

  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();

  if(useMC)
    task->SetHasMCData(kTRUE);
  else{
    task->SetHasMCData(kFALSE);
  }

  if(TRDtrigger==0) task->SelectCollisionCandidates(AliVEvent::kINT7);
  else task->SelectCollisionCandidates(AliVEvent::kTRD);


  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFEtask";
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
AliAnalysisTask *RegisterTaskPID2TRD(Bool_t useMC, Bool_t isAOD, Int_t tpcCls=120, Int_t tpcClsPID = 80,
				     Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0,
				     Double_t *tpcdEdxcutlow=NULL,
				     Double_t *tpcdEdxcuthigh=NULL,
				     Double_t tofs=3., Int_t tofm=0,
				     Int_t itshitpixel = AliHFEextraCuts::kBoth, Int_t icent=1,
				     Int_t TRDtrigger=0,Int_t trdl=6, Int_t trde=4,TString detector){

    gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEpPbTRD.C");
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    TString detused=detector.Copy();
    TPRegexp(",").Substitute(detused,"","g");
    printf("detectors in use %s %s \n",detector.Data(),detused.Data());

    if (icent == 2) TString cesti("V0M");
    else if (icent == 3) TString cesti("CL1");
    else if (icent == 4) TString cesti("ZNA");
    else TString cesti("V0A");

    Int_t idcaxy = (Int_t)(dcaxy*10.);
    Int_t idcaz = (Int_t)(dcaz*10.);
    Int_t tpclow = 0;
    if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
    Int_t itofs = (Int_t)(tofs*10.);
    Int_t ipixelany = itshitpixel;

    TString appendix(TString::Format("mbTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%dtrdtrg%dtrdl%itrde%iPID%scent%s",tpcCls,
				     tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,tofm,ipixelany,TRDtrigger,trdl,trde,detused.Data(),cesti.Data()));
    printf("Add macro appendix %s\n", appendix.Data());

   AliAnalysisTaskHFE *task = ConfigHFEpPbTRD(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz,
					      tpcdEdxcutlow,tpcdEdxcuthigh, tofs, tofm, itshitpixel, icent, -0.8, 0.8,
					      TRDtrigger,trdl,trde,detector);
  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();

  if (useMC)
    task->SetHasMCData(kTRUE);
  else{
    task->SetHasMCData(kFALSE);
  }

  if(TRDtrigger<2) task->SelectCollisionCandidates(AliVEvent::kINT7);
  else task->SelectCollisionCandidates(AliVEvent::kTRD);

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
