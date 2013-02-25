AliAnalysisTask *AddTaskHFEnpepPb(){
  // Switches
  Bool_t kPhotonic                   = kTRUE;

  // Default settings (TOF-TPC 7 TeV pp paper)
  const int	kDefTPCcl	= 120;
  const int	kDefTPCclPID	=  80;
  const int	kDefITScl	=   4;
  const double	kDefDCAr	=   1.;
  const double	kDefDCAz	=   2.;
  const double	kDefTPCs	=   0.026;
  const double	kDefTPCu	=   2.8;
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

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_hfe_HFE", "No analysis manager found.");
    return 0;
  }

  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH) MCthere=kFALSE;

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //@@ 0 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // Remember: (e.g. for 120 clusters) This is only an example to show how to calculate!
  // effective mean m0=-0.113
  // width w=1.03
  // cut --> nsigma*w + m0
  // e.g. 45% 0.122 sigma --> 0.122*1.03-0.113 = 0.013

  if(kPhotonic)
  {
    // Task for pp @ 7 TeV
    //RegisterTaskNPPID2( MCthere, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, kDefTPCs, kDefTPCu, kDefTOFs, AliHFEextraCuts::kFirst, kTRUE, kassETA, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, kassTPCSminus, kassTPCSplus);

    // task for pPb @ 5.023 TeV
    //RegisterTaskNPPID2( MCthere, 110, 80, 3, 1, 2, 0.045, 2.685, 3, AliHFEextraCuts::kBoth, kTRUE, kassETA, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, kassTPCSminus, kassTPCSplus);

    //task for pp @ 2.76 TeV
    RegisterTaskNPPID2( MCthere, 110, 80, 3, 1, 2, -1.35, 2.45, 3, AliHFEextraCuts::kBoth, kTRUE, kassETA, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, kassTPCSminus, kassTPCSplus);
  }

  return NULL;
}

//===============================================================================

AliAnalysisTask *RegisterTaskNPPID2( Bool_t useMC, Int_t tpcCls=120, Int_t tpcClsPID=80, Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, Double_t tpcs=-0.0113, Double_t tpcu=3.09, Double_t tofs=3., Int_t itshitpixel=0, Bool_t kNoPhotonic=kFALSE, Double_t assETA=0.8, Int_t assITS=2, Int_t assTPCcl=100, Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0, Double_t assTPCSminus=-3.0, Double_t assTPCSplus=3.0 )
{
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEnpepPb.C");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEnpepPb(useMC, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, tpcs, tpcu, tofs, itshitpixel, kNoPhotonic, assETA, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, assTPCSplus);
  task->SetESDAnalysis();

  if (useMC)	task->SetHasMCData(kTRUE);
  else		task->SetHasMCData(kFALSE);

  task->SelectCollisionCandidates(AliVEvent::kINT7);

  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz  = (Int_t)(dcaz*10.);
  Int_t itpcs  = (Int_t)(tpcs*1000.);
  Int_t itofs  = (Int_t)(tofs*10.);
  Int_t ipixelany   = itshitpixel;

  Int_t iNoPhotonic = 0;
  if(kNoPhotonic) iNoPhotonic = 1;

  Int_t iassETA  = (Int_t)(assETA*10);
  Int_t iassDCAr = (Int_t)(assDCAr*10);
  Int_t iassDCAz = (Int_t)(assDCAz*10);
  Int_t iassTPCSminus = (Int_t)(assTPCSminus*10);
  Int_t iassTPCSplus  = (Int_t)(assTPCSplus*10);

  TString appendix(TString::Format("Photonic%dEta%dITS%dTPCcl%dTPCPIDcl%dDCAr%dDCAz%dTPCSminus%dTPCSplus%d", iNoPhotonic, iassETA, assITS, assTPCcl, assTPCPIDcl, iassDCAr, iassDCAz, iassTPCSminus, iassTPCSplus));
  //TString appendix(TString::Format("TPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dpa%dNoPhotonic%d",tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,itpcs,itofs,ipixelany,iNoPhotonic));
  printf("Add macro appendix %s\n", appendix.Data());

  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFEpPb";
  containerName += appendix.Data();
  printf("container name: %s\n", containerName.Data());

  //create data containers
  task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data() ));
  task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput(task,  0, cinput );

  mgr->AddTask(task);

  return NULL;
}
