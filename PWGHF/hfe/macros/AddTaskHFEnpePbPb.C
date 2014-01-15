AliAnalysisTask *AddTaskHFEnpePbPb(Bool_t MCthere, 
                                   Bool_t isAOD,
                                   Bool_t kNPERef = kTRUE,
                                   Bool_t kNPEkAny = kFALSE,
				   Bool_t kNPERefMCf =  kTRUE,
				   Bool_t kNPERefTPConly =  kTRUE)
{
  // Default settings (TOF-TPC PbPb)
  const int	kDefTPCcl	= 130;
  const int	kDefTPCclPID	=  80;
  const int kDefTPCclshared = 1.1;
  const int	kDefITScl	=   4;
  const int kDefITSchi2percluster = -1; // cleanup removes badly matching tracks - effects high pt  (cut value = 36)
  const double	kDefDCAr	=   1.;
  const double	kDefDCAz	=   2.;
  const double	kDefTOFs	=   3.;
  const double  kDefEtaIncMin = -0.8;
  const double  kDefEtaIncMax = 0.8;
  const Bool_t   etacorrection   = kFALSE;
  const Bool_t   multicorrection = kTRUE;

  Double_t dEdxhm[12] = {3.11,3.11,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};  
  Double_t tpcl1[12]  = {-0.14,-0.14,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3};  

  // Default setting for the associated electron for the NonPhotonic Analysis
  const double	kassETAm = -0.8;
  const double	kassETAp = 0.8;
  const int	kassITS		=   2;
  const int	kassTPCcl	= 60;
  const int	kassTPCPIDcl	=  60;
  const double	kassDCAr	=   1.0;
  const double	kassDCAz	=   2.0;
  const double	kassTPCSminus	=  -3.0;
  const double	kassTPCSplus	=   3.0;

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_hfe_HFE", "No analysis manager found.");
    return 0;
  }

  //mgr->AddClassDebug("AliAnalysisTaskHFE",12);
 

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //@@ 0 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Double_t dEdxaclm[12], dEdxachm[12];
  for(int icent = 0; icent < 12; icent++){
    dEdxaclm[icent] = kassTPCSminus;
    dEdxachm[icent] = kassTPCSplus;
  }

  if(kNPERef){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
    RegisterTaskNPEPbPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kFALSE);
  }

  if(kNPEkAny){
    // **************************************************************
    // 
    // task for kAny instead of kBoth
    //
    // **************************************************************
    RegisterTaskNPEPbPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kAny, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kFALSE);
  }
  if(kNPERefMCf){
    // **************************************************************
    // 
    // Reference task + MC fake rejected
    //
    // **************************************************************
    RegisterTaskNPEPbPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE, kTRUE);
  }

  if(kNPERefTPConly){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
    RegisterTaskNPEPbPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, 0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE, kFALSE);
  }

  return NULL;
}

//===============================================================================
AliAnalysisTask *RegisterTaskNPEPbPb(Bool_t useMC, Bool_t isAOD, 
               Int_t tpcCls=120, Int_t tpcClsPID=80, 
               Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
               Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL, 
               Double_t tofs=3., Int_t itshitpixel =AliHFEextraCuts::kBoth,
               Double_t itschi2percluster = -1, Double_t tpcsharedcluster = 1.1,
               Bool_t etacorr=kFALSE, Bool_t multicorr = kFALSE,
               Double_t etaIncMin = -0.8, Double_t etaIncMax = 0.8,
               Double_t assETAm=-0.8, Double_t assETAp=0.8, Int_t assITS=2, Int_t assTPCcl=100,
               Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
               Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
				     Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE, Bool_t rejectMCFake = kFALSE)
{

  //
  // Cuts on the inclusive leg
  //
  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t tpclow = 0;
  if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;
  Int_t imult = multicorr ? 1 : 0;

  //
  // Cuts on the associated leg
  //
  Int_t iassDCAr = (Int_t)(assDCAr*10);
  Int_t iassDCAz = (Int_t)(assDCAz*10);
  Int_t iassTPCSplus  = assTPCSplus ? (Int_t)(assTPCSplus[0]*1000) : 0;
  Int_t icat1 = useCat1Tracks ? 1 : 0;
  Int_t icat2 = useCat2Tracks ? 1 : 0;
  Int_t irejectMCFake = rejectMCFake ? 1 : 0;

  TString appendix(TString::Format("SPD%d_incTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dm%d_photTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%dMCf%d",ipixelany,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,imult,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,irejectMCFake));

  printf("Add macro appendix %s\n", appendix.Data());

  if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpePbPb"))gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFEnpePbPb.C");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEnpePbPb(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, itschi2percluster, tpcsharedcluster, etacorr, multicorr, etaIncMin, etaIncMax,
					      assETAm, assETAp, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, assTPCSplus, useCat1Tracks, useCat2Tracks,rejectMCFake);
  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();

  if (useMC)	task->SetHasMCData(kTRUE);
  else		task->SetHasMCData(kFALSE);

  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFE";
  containerName += appendix.Data();
  printf("container name: %s\n", containerName.Data());

  //create data containers
  task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
  task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput(task,  0, cinput );

  mgr->AddTask(task);

  return NULL;
}
