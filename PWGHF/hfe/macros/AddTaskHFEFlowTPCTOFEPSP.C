AliAnalysisTask *AddTaskHFEFlowTPCTOFEPSP(UInt_t trigger=131073,Int_t aodfilter=16,Bool_t scalarProduct=kFALSE,Bool_t cutPileup=kFALSE,Int_t tpcCls=110, Double_t tpcClsr=60, Int_t tpcClspid=80, Int_t itsCls=4, Int_t pixellayer=2, Double_t dcaxy=100,Double_t dcaz=200, Double_t tofsig=30., Double_t tpceff=50., Int_t vzero=1,Int_t debuglevel=0,Double_t etarange=80,Double_t ITSclustersback=0,Double_t minTPCback=-2.0,Double_t maxTPCback=5.0){

  //
  // Define TPC cut for 2011 data
  //
  Double_t tpcdedx[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-100
  // -0.2 0-5%
  // -0.15 5-10%
  // -0.1 10-20%
  // -0.0 20-30%, sigma=1.17
  // 0.156 30-40%, sigma=1.2
  // 0.19 40-50%, sigma=1.2
  // 0.2 50-60%
  // 0.2 60-80% list_t65536f16TPC110r60p80s11km0ITS4C36Pi2DCAr100z200TOF30TPCe50V1D0er8i0t-20t50
  tpcdedx[0]=-0.2;
  tpcdedx[1]=-0.15;
  tpcdedx[2]=-0.1;
  tpcdedx[3]=0.0;
  tpcdedx[4]=0.156;
  tpcdedx[5]=0.19;
  tpcdedx[6]=0.2;
  tpcdedx[7]=0.2;
  if(TMath::Abs(tpceff-55)<0.01) {
    tpcdedx[0]=-0.365;
    tpcdedx[1]=-0.314;
    tpcdedx[2]=-0.267;
    tpcdedx[3]=-0.165;
    tpcdedx[4]=-0.022;
    tpcdedx[5]= 0.01;
    tpcdedx[6]= 0.018;
    tpcdedx[7]= 0.018;
  }
  if(TMath::Abs(tpceff-45)<0.01) {
    tpcdedx[0]=-0.062;
    tpcdedx[1]=-0.015;
    tpcdedx[2]=0.035;
    tpcdedx[3]=0.131;
    tpcdedx[4]=0.278;
    tpcdedx[5]=0.32;
    tpcdedx[6]=0.32;
    tpcdedx[7]=0.32;
  }
  if(TMath::Abs(tpceff-60)<0.01) {
    tpcdedx[0]=-0.518;
    tpcdedx[1]=-0.47;
    tpcdedx[2]=-0.42;
    tpcdedx[3]=-0.315;
    tpcdedx[4]=-0.178;
    tpcdedx[5]=-0.145;
    tpcdedx[6]=-0.135;
    tpcdedx[7]=-0.135;
  }
  if(TMath::Abs(tpceff-40)<0.01) {
    tpcdedx[0]=0.09;
    tpcdedx[1]=0.14;
    tpcdedx[2]=0.188;
    tpcdedx[3]=0.28;
    tpcdedx[4]=0.43;
    tpcdedx[5]=0.462;
    tpcdedx[6]=0.473;
    tpcdedx[7]=0.473;
  }

  // Name
  TString appendixx(TString::Format("t%df%ds%dp%dTPC%dr%dp%dITS%dPi%dDCAr%dz%dTOF%dTPCe%dV%dD%der%di%dt%dt%d",(Int_t)trigger,aodfilter,(Int_t)scalarProduct,(Int_t)cutPileup,tpcCls,(Int_t)tpcClsr,tpcClspid,itsCls,(Int_t) pixellayer,(Int_t) dcaxy,(Int_t)dcaz,(Int_t) tofsig,(Int_t)tpceff,vzero,debuglevel,(Int_t)(etarange*0.1),(Int_t)ITSclustersback,(Int_t)(minTPCback*10.0),(Int_t)(maxTPCback*10.0)));
  //TString appendixx("tpctofv2");
  

  //set config file name
  TString configFile("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TOFTPC.C");
  //TString configFile("/d/alice12/bailhache/AliRootInstallations/29_11_2012/AliRoot/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TOFTPC.C");
  TString checkconfig="ConfigHFE_FLOW_TOFTPC";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskFlowTPCTOFEPSP *task = ConfigHFE_FLOW_TOFTPC(kFALSE,appendixx,trigger,aodfilter,scalarProduct,cutPileup,tpcCls, tpcClsr, tpcClspid, itsCls, pixellayer, dcaxy, dcaz,tofsig,&tpcdedx[0],vzero,debuglevel,etarange,kFALSE,ITSclustersback,minTPCback,maxTPCback);  

  
  task->SetNbBinsCentralityQCumulant(4);
  //task->SetBinCentralityLess(0,0.0);
  task->SetBinCentralityLess(0,0.0);
  task->SetBinCentralityLess(1,10.0);
  task->SetBinCentralityLess(2,20.0);
  task->SetBinCentralityLess(3,40.0);
  task->SetBinCentralityLess(4,50.0);
  //task->SetBinCentralityLess(5,60.0);
  //task->SetBinCentralityLess(7,80.0);

  task->SetHFEVZEROEventPlane(0x0);
  //AliLog::SetClassDebugLevel("AliAnalysisTaskHFEFlow",3);

  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendixx.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer(Form("list_%s",appendixx.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return NULL;

  
}
