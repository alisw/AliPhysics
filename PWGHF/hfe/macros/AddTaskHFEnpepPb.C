AliAnalysisTask *AddTaskHFEnpepPb(Bool_t isAOD, Bool_t useMC, Bool_t kNPERef = kFALSE, Bool_t kNPECent = kFALSE){

  // Default settings (TOF-TPC pPb)
  const int	kDefTPCcl	= 110;
  const int	kDefTPCclPID	=  80;
  const int	kDefITScl	=   4;
  const double	kDefDCAr	=   1.;
  const double	kDefDCAz	=   2.;
  const double	kDefTPCs	=   0.0;
  const double	kDefTPCu	=   3.0;
  const double	kDefTOFs	=   3.;

  // Default setting for the associated electron for the NonPhotonic Analysis
  const double	kassETA		=   0.8;
  const int	kassITS		=   2;
  const int	kassTPCcl	= 100;
  const int	kassTPCPIDcl	=  80;
  const double	kassDCAr	=   1.0;
  const double	kassDCAz	=   2.0;
  const double	kassTPCSminus	=  -3.0;
  const double	kassTPCSplus	=   3.0;

  enum{
          kHFEV0A = 1,
          kHFEV0M = 2,
          kHFECL1 = 3,
          kHFEZNA = 4
  };

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_hfe_HFE", "No analysis manager found.");
    return 0;
  }

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //@@ 0 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Double_t dEdxlm[12], dEdxhm[12], dEdxaclm[12], dEdxachm[12];
  for(int icent = 0; icent < 12; icent++){
    dEdxlm[icent] = kDefTPCs;
    dEdxhm[icent] = kDefTOFs;
    dEdxaclm[icent] = kassTPCSminus;
    dEdxachm[icent] = kassTPCSplus;
  }

  if(kNPERef){
    RegisterTaskNPEpPb(isAOD, useMC, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, dEdxlm, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassETA, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm);
  }
  if(kNPECent){
    RegisterTaskNPEpPb(isAOD, useMC, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, dEdxlm, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0M, kassETA, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm);
    RegisterTaskNPEpPb(isAOD, useMC, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, dEdxlm, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFECL1, kassETA, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm);
    RegisterTaskNPEpPb(isAOD, useMC, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, dEdxlm, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEZNA, kassETA, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm);
  }
  return NULL;
}

//===============================================================================

AliAnalysisTask *RegisterTaskNPEpPb(Bool_t isAOD, Bool_t useMC, 
               Int_t tpcCls=120, Int_t tpcClsPID=80, 
               Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
               Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL, 
               Double_t tofs=3., Int_t itshitpixel =AliHFEextraCuts::kBoth, Int_t icent=1,
               Double_t assETA=0.8, Int_t assITS=2, Int_t assTPCcl=100, 
               Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0, 
               Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL )
{

  Int_t iassETA  = (Int_t)(assETA*10);
  Int_t iassDCAr = (Int_t)(assDCAr*10);
  Int_t iassDCAz = (Int_t)(assDCAz*10);
  Int_t iassTPCSminus = assTPCSminus ? (Int_t)(assTPCSminus[0]*10) : 0;
  Int_t iassTPCSplus  = assTPCSplus ? (Int_t)(assTPCSplus[0]*10) : 0;

  if (icent == 2) TString cesti("V0M");
  else if (icent == 3) TString cesti("CL1");
  else if (icent == 4) TString cesti("ZNA");
  else TString cesti("V0A");

  TString appendix(TString::Format("PhotonicEta%dITS%dTPCcl%dTPCPIDcl%dDCAr%dDCAz%dTPCSminus%dTPCSplus%d%s", iassETA, assITS, assTPCcl, assTPCPIDcl, iassDCAr, iassDCAz, iassTPCSminus, iassTPCSplus, cesti.Data()));
  printf("Add macro appendix %s\n", appendix.Data());

  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEnpepPb.C");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEnpepPb(isAOD, useMC, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, icent, assETA, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, assTPCSplus);
  task->SetESDAnalysis();

  if (useMC)	task->SetHasMCData(kTRUE);
  else		task->SetHasMCData(kFALSE);

  task->SelectCollisionCandidates(AliVEvent::kINT7);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFEtask";
  containerName += appendix.Data();
  printf("container name: %s\n", containerName.Data());

  //create data containers
  task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data() ));
  task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput(task,  0, cinput );

  mgr->AddTask(task);

  return NULL;
}
