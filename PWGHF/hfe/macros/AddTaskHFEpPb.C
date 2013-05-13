AliAnalysisTask *AddTaskHFEpPb(Bool_t isMC = kFALSE,
                               Bool_t isAOD = kFALSE,
                               Bool_t kTPCTOF_Ref = kTRUE,
                               Bool_t kTPC_Only = kFALSE,
                               Bool_t kTPCTOF_Cent = kFALSE,
                               Bool_t kTPCTOF_Sys = kFALSE,
			       Bool_t kTPCTOFTRD_Ref = kFALSE,
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
  
  Double_t dEdxlm[12] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};  //  50%
  Double_t dEdxhm[12] = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
  // For systematics:
  Double_t tpcl0[12]  = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};  //  84%
  Double_t tpcl1[12]  = {-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53};  //  70%
  Double_t tpcl2[12]  = {-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26};  //  60%
  Double_t tpcl3[12]  = {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};  //   40%/   40%

  
  // For TPC only
  const double kDefTPCs1 = 0.0;
  const double kDefTPCu1 = 3.0;

  // eta cut
  Int_t etacut=0; // eta cut off
  
  if(kTPC_Only){
    // for the moment (09.02.2013) use the same as TOF-TPC. Refine later on
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,-0.8,0.8); // 50%
  }
  
  //------------------------------//
  //   TPC-TOF analysis
  //------------------------------//
  
  if(kTPCTOF_Ref){
    // Reference task
    // with centrality V0A
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1,-0.8,0.8); // 50%
  }

  if (kTPCTOF_Cent){
    // For the moment we set V0A as the reference centrality estimator, and test the other ones:
    // 1: V0A, 2: V0M, 3: CL1, 4: ZNA
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,2); // V0M
    RegisterTask(isMC,isAOD, kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,3); // CL1
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,4); // ZNA

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
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,4,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,5,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,3,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,4,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst);
    RegisterTask(isMC,isAOD,kDefTPCcl,kDefTPCclPID,5,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst);

    // TPC clusters
    RegisterTask(isMC,isAOD,100,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
    RegisterTask(isMC,isAOD,120,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth);
    // TPC clusters PID
    RegisterTask(isMC,isAOD,kDefTPCcl,100,kDefITScl,kDefDCAr,kDefDCAz,&dEdxlm[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth); 
  }

   if(kTPCTOFTRD_Ref){

      RegisterTaskPID2(isMC,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger);
   }

  if(kTPCTOFTRD_PID){
      // without TOF
      RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,0,4,"TPC",etacut);
      RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,6,4,"TRD,TPC",etacut);
      RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,6,4,"TPC,TRD",etacut);
      RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,6,2,"TRD,TPC",etacut);
      RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,5,4,"TRD,TPC",etacut);
      RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,5,4,"TPC,TRD",etacut);
      RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,0,5,2,"TRD,TPC",etacut);
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

  
  TString appendix(TString::Format("centTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%dce%seta%d%d",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,tofm,ipixelany,cesti.Data(),iEtaMin,iEtaMax));
  printf("Add macro appendix %s\n", appendix.Data());


  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEpPb.C");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEpPb(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					  tpcdEdxcutlow,tpcdEdxcuthigh,
					  tofs,tofm,itshitpixel,icent,EtaMin,EtaMax);

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

AliAnalysisTask *RegisterTaskPID2TRD(Bool_t useMC, Bool_t isAOD, Int_t tpcCls=120, Int_t tpcClsPID = 80,
				  Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
				  Double_t tpcs=-0.0113, Double_t tpcu=3.09, Double_t tofs=3., 
				  Int_t tofm=0,
				  Int_t itshitpixel = AliHFEextraCuts::kBoth, 
				  Bool_t withetacorrection = kTRUE,
				     Int_t TRDtrigger=0,Int_t trdl=6, Int_t trde=4,TString detector,Int_t etacut){

    TString detused=detector.Copy();
    TPRegexp(",").Substitute(detused,"","g");
    printf("detectors in use %s %s \n",detector.Data(),detused.Data());


  gROOT->LoadMacro("$TRAIN_ROOT/util/hfe/configs/ConfigHFEmbpPbTRD.C");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEmbpPbTRD(useMC, isAOD, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz,
					     tpcs,tpcu,tofs,tofm,3.,kFALSE,kTRUE,kFALSE,itshitpixel,withetacorrection,0,TRDtrigger,trdl,trde,detector,etacut);

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
  
  TString appendix(TString::Format("mbTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%dtrdtrg%dtrdl%itrde%iPID%setacut%i",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,itpcs,itofs,tofm,ipixelany,TRDtrigger,trdl,trde,detused.Data(),etacut));
  printf("Add macro appendix %s\n", appendix.Data());

  //create data containers
  task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), 
					      TList::Class(), AliAnalysisManager::kOutputContainer, 
					      Form("HFEtask%s.root", appendix.Data())));
  task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(), 
					      AliAnalysisManager::kOutputContainer, 
					      Form("HFEtask%s.root", appendix.Data())));
  mgr->ConnectInput(task,  0, cinput );

  mgr->AddTask(task);
  return NULL;
}