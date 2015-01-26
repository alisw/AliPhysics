AliAnalysisTask *AddTaskHFEpPb(Bool_t MCthere = kFALSE,
                               Bool_t isAOD = kFALSE,
                               Bool_t kTPCTOF_Ref = kTRUE,
                               Bool_t kTPC_Ref = kFALSE,
                               Bool_t kTPCTOF_Cent = kFALSE,
                               Bool_t kTPCTOF_Sys = kFALSE,
                               Bool_t kTPCTOF_TPCPID = kFALSE,
                               Bool_t kITS_Sys = kFALSE,
                               Bool_t kTPCTOFTRD_Ref = kFALSE,
                               Bool_t kTPCTOFTRD_mbPID = kFALSE,
			       Bool_t kTPCTOFTRD_PID1 = kFALSE,
                               Bool_t kTPCTOFTRD_PID2 = kFALSE,
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
  // 15.04.2013 The ITS cut changes now officially from 3 to 4
  const int kDefITScl = 4;
  const int kDefTPCcl = 110;
  const int kDefTPCclPID = 80;
  const double kDefDCAr = 1.;
  const double kDefDCAz = 2.;

  // 24.02.2013 After the latest fits, mean is 0 and sigma is 1
  const double kDefTPCs = 0.0;
  const double kDefTPCu = 3.0;   
  const double kDefTOFs = 3.0;

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

  // 24.02.2013 From now on, use task with centrality by default and keep one without for cross check
  // 15.04.2013 We do not need the minimum bias task any more, take out by default

  if(kTPC_Ref){
    // Reference, 69% TPC PID, with centrality V0A
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1);
  }

  //------------------------------//
  //   TPC-TOF analysis
  //------------------------------//

  if(kTPCTOF_Ref){
    // Reference task, 69% TPC PID, with centrality V0A
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1);
  }

  if (kTPCTOF_TPCPID){
    // TPC PID
    if (!MCthere){
      //TPC only
      //RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl2[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1);
      //RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl0[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1);

      //TPC-TOF
      //RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
      //RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl2[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
      // new pT bins
      RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl0[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
      RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl2[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
      RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl3[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
      RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl4[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
      RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl5[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
      RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl6[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
      RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl7[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    }
  }

  if (kITS_Sys){
    // ITS hits and SPD request
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,3,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,5,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,3,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,4,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,5,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kFirst,1); 
  }

  if (kTPCTOF_Cent){
    // For the moment we set V0A as the reference centrality estimator, and test the other ones:
    // 1: V0A, 2: V0M, 3: CL1, 4: ZNA
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,2); // V0M
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,3); // CL1
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,4); // ZNA
  }

  if(kTPCTOF_Sys){
    // TOF PID
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],2.0,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],2.5,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],3.5,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],4.0,0,AliHFEextraCuts::kBoth,1); 

    // TPC clusters
    RegisterTask(MCthere,isAOD,80,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,90,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,100,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,120,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,130,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,140,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 

    // TPC clusters PID
    RegisterTask(MCthere,isAOD,kDefTPCcl,60,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,90,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,100,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
    RegisterTask(MCthere,isAOD,kDefTPCcl,110,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],kDefTOFs,0,AliHFEextraCuts::kBoth,1); 
  }
  
  //------------------------------//
  //   TPC-TOF-TRD analysis
  //------------------------------//

  if(kTPCTOFTRD_Ref){
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,kDefTOFs,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger);
    RegisterTaskPID2(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger);
  }

  if(kTPCTOFTRD_mbPID){
    // without TOF
    RegisterTaskPID2mbTRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger,6,4,"TRD,TPC",etacut);
    RegisterTaskPID2mbTRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger,6,4,"TPC,TRD",etacut);
    RegisterTaskPID2mbTRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger,6,2,"TRD,TPC",etacut);
    RegisterTaskPID2mbTRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger,5,4,"TRD,TPC",etacut);
    RegisterTaskPID2mbTRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger,5,4,"TPC,TRD",etacut);
    RegisterTaskPID2mbTRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,kDefTPCs,kDefTPCu,0,0,AliHFEextraCuts::kBoth,kFALSE,TRDtrigger,5,2,"TRD,TPC",etacut);
  }

  if(kTPCTOFTRD_PID1){
    // without TOF
    RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,TRDtrigger,1,5,4,"TRD,TPC");
    RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,TRDtrigger,1,5,4,"TPC,TRD");
    RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,TRDtrigger,1,6,4,"TRD,TPC");
    RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,TRDtrigger,1,6,4,"TPC,TRD");
  }

  if(kTPCTOFTRD_PID2){
    // without TOF
    RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,TRDtrigger,2,5,4,"TRD,TPC");
    RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,TRDtrigger,2,5,4,"TPC,TRD");
    RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,TRDtrigger,2,6,4,"TRD,TPC");
    RegisterTaskPID2TRD(MCthere,isAOD,kDefTPCcl,kDefTPCclPID,kDefITScl,kDefDCAr,kDefDCAz,&tpcl1[0],&dEdxhm[0],0,0,AliHFEextraCuts::kBoth,1,TRDtrigger,2,6,4,"TPC,TRD");
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
                  Int_t itshitpixel = AliHFEextraCuts::kBoth, Int_t icent=1){

  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t tpclow = 0;
  if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;

  printf("Argument passed to function to determine the centrality estimator %i\n", icent);

  if (icent == 2) TString cesti("V0M");
  else if (icent == 3) TString cesti("CL1");
  else if (icent == 4) TString cesti("ZNA");
  else TString cesti("V0A");

  TString appendix(TString::Format("centTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%dce%s",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,tofm,ipixelany,cesti.Data()));
  printf("Add macro appendix %s\n", appendix.Data());
  printf("Centrality estimator %s\n", cesti.Data());

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigHFEpPb.C");
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
                  Double_t tpcs=-0.0113, Double_t tpcu=3.09, 
                  Double_t tofs=3., Int_t tofm=0,
                  Int_t itshitpixel = AliHFEextraCuts::kBoth, 
                  Bool_t withetacorrection = kTRUE,
                  Int_t TRDtrigger=0){

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigHFEmbpPb.C");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEmbpPb(useMC, isAOD, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz,
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

  if(TRDtrigger<2) task->SelectCollisionCandidates(AliVEvent::kINT7);
  else task->SelectCollisionCandidates(AliVEvent::kTRD);

  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t itpcs = (Int_t)(tpcs*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;
  Int_t ietacorr = 0;
  if(withetacorrection) ietacorr = 1;
  
  TString appendix(TString::Format("mbTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%ipa%dtrdt%d",tpcCls,
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

//===============================================================================
// TRD PID
//===============================================================================

AliAnalysisTask *RegisterTaskPID2mbTRD(Bool_t useMC, Bool_t isAOD, Int_t tpcCls=120, Int_t tpcClsPID = 80,
                  Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0,
                  Double_t tpcs=-0.0113, Double_t tpcu=3.09,
                  Double_t tofs=3., Int_t tofm=0,
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

  TString appendix(TString::Format("mbTPCc%dp%dITS%dr%dz%dTPCs%dTOFs%dm%ipa%dTRDt%dl%ie%i%se%i",tpcCls,
				   tpcClsPID,itsCls,idcaxy,idcaz,itpcs,itofs,tofm,ipixelany,TRDtrigger,trdl,trde,detused.Data(),etacut));
  printf("Add macro appendix %s\n", appendix.Data());

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigHFEmbpPbTRD.C");
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

  if(TRDtrigger<2) task->SelectCollisionCandidates(AliVEvent::kINT7);
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
				     Int_t TRDtrigger=0,Int_t trdpidmethod=1, Int_t trdl=6, Int_t trde=4,TString detector){

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigHFEpPbTRD.C");
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

  TString appendix(TString::Format("TPCc%dp%dITS%dr%dz%dTPCs%dTOFs%dm%ipa%dTRDt%dm%il%ie%i%sc%s",tpcCls,
                   tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,tofm,ipixelany,TRDtrigger,trdpidmethod,trdl,trde,detused.Data(),cesti.Data()));
  printf("Add macro appendix %s\n", appendix.Data());

  AliAnalysisTaskHFE *task = ConfigHFEpPbTRD(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz,
                             tpcdEdxcutlow,tpcdEdxcuthigh, tofs, tofm, itshitpixel, icent, -0.8, 0.8, TRDtrigger,trdpidmethod,trdl,trde,detector);

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
