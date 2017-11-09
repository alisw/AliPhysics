AliAnalysisTask *AddTaskHFEnpePbPb(Bool_t MCthere, 
                                   Bool_t isAOD = kTRUE,
				   Bool_t kNPERef = kTRUE,
				   Bool_t kNPEkAny = kFALSE,
                                   Bool_t newCentralitySelection = kTRUE, // kTRUE: new framework used; kFALSE: old framework used 
				   Bool_t kNPERefTPConly = kFALSE,
				   Bool_t kNPETOFITS = kTRUE,
				   Bool_t kNPETOFlast = kFALSE,
				   Bool_t kNPEw      = kFALSE,
				   Bool_t kNPEkf  = kFALSE)		   
  
{
  // Default settings (TOF-TPC PbPb)
  const int	kDefTPCcl	= 120;
  const int	kDefTPCclPID	=  80;
  const int kDefTPCclshared = 1.1;
  const int	kDefITScl	=   4;
  const int kDefITSchi2percluster = -1; // cleanup removes badly matching tracks - effects high pt  (cut value = 36)
  const double	kDefDCAr	=   1.;
  const double	kDefDCAz	=   2.;
  const double	kDefTOFs	=   3.;
  const double	kDefITSs	=   2.;
  const double  kDefEtaIncMin = -0.8;
  const double  kDefEtaIncMax = 0.8;
  const Bool_t   etacorrection   = kFALSE;
  const Bool_t   multicorrection = kFALSE;

  // --- TPC nsigma max and min cut ---
  Double_t dEdxhm[12] = {3.11,3.11,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};  
  Double_t tpcl1[12]  = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};  

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
  // --- Centrality selection ---
  const int     centrMin        = 0;
  const int     centrMax        = 10;

  Int_t kWei = -1;
  /*
  if (MCthere) kWei = 9;    // default Pb-Pb
  enum {

    k11a10abisweiData = 6,  // LHC11a10abis weights 
    k11a10bplusweiData = 7,     // LHC11a10b_plus weights 
    k11a10bpluswei11a10abis   = 8,  // LHC11a10b_plus weights for LHC11a10abis
    k11a10bbisweiData = 9,  // LHC11a10bbis weights 
  };
  */

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


  const Bool_t isBeauty = kFALSE; // should be false to prevent inclusive analysis

  if(kNPERef){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
  }

  if(kNPETOFlast){
    // **************************************************************
    // 
    // Apply TOF after TPC for mismatch background studies
    //
    // **************************************************************
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kTRUE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
  }

  if(kNPEkAny){
    // **************************************************************
    // 
    // task for kAny instead of kBoth
    //
    // **************************************************************
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kAny, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
  }

  if(kNPEw && MCthere){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE, 0,kWei);
  }

  if(kNPERefTPConly){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, 0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
  }


 if(kNPETOFITS){
   // **************************************************************
   // 
   // Reference task
   //
   // **************************************************************
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm,  kDefTOFs,  kDefITSs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
 }

 
  if(kNPEkf){
    // **************************************************************
    // 
    // Use KF particle
    //
    // **************************************************************
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1,2,kFALSE,kFALSE,kFALSE,kTRUE);
  }

  
  return NULL;

}

//===============================================================================

//===============================================================================
AliAnalysisTask *RegisterTaskNPEPbPb(
                                     Int_t centrMin = 0, Int_t centrMax = 100,
                                     Bool_t newCentralitySelection = kTRUE, // kTRUE: new framework used; kFALSE: old framework used 
                                     Bool_t useMC, Bool_t isAOD, Bool_t beauty,
				     Int_t tpcCls=120, Int_t tpcClsPID=80, 
				     Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
				     Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL, 
				     Double_t tofs=3., Double_t itss=0., Int_t itshitpixel =AliHFEextraCuts::kBoth,
				     Double_t itschi2percluster = -1, Double_t tpcsharedcluster = 1.1,
				     Bool_t etacorr=kFALSE, Bool_t multicorr = kFALSE, Bool_t toflast = kFALSE,
				     Double_t etaIncMin = -0.8, Double_t etaIncMax = 0.8,
				     Double_t assETAm=-0.8, Double_t assETAp=0.8, Int_t assITS=2, Int_t assTPCcl=100,
				     Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
				     Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
				     Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE,
				     Int_t weightlevelback = -1,Int_t wei = 2,
				     Bool_t releasemcvx = kFALSE,
				     Bool_t ipCharge = kFALSE,
				     Bool_t ipOpp = kFALSE,
				     Bool_t usekfparticle = kFALSE)
{

  //
  // Cuts on the inclusive leg
  //
  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t tpclow = 0;
  if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t iitss = (Int_t)(itss*10.);
  Int_t ipixelany = itshitpixel;
  Int_t imult = multicorr ? 1 : 0;
  Int_t itofpos = toflast ? 1 : 0;

  //
  // Cuts on the associated leg
  //
  Int_t iassDCAr = (Int_t)(assDCAr*10);
  Int_t iassDCAz = (Int_t)(assDCAz*10);
  Int_t iassTPCSplus  = assTPCSplus ? (Int_t)(assTPCSplus[0]*1000) : 0;
  Int_t icat1 = useCat1Tracks ? 1 : 0;
  Int_t icat2 = useCat2Tracks ? 1 : 0;
  
  Bool_t nondefaultcentr = kFALSE;

  TString cweightsback("");
  if(weightlevelback>=0) {
    cweightsback += "Wa";
    cweightsback += weightlevelback;
  }

  TString cmvx("");
  if(releasemcvx) {
    cmvx += "MCVR";
  }

  TString kfp("");
  if(usekfparticle) {
    kfp += "kf";
  }
  
  if(beauty) {
    if(ipCharge && ipOpp) TString cbeauty("BeautyIPopp");
    else if(ipCharge) TString cbeauty("BeautyIPcrg");
    else if(!ipCharge) TString cbeauty("Beauty");
    else TString cbeauty("BeautyWrong");
  }
  else TString cbeauty("");
  
  TString appendix(TString::Format("SPD%d_incTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dITSs%dm%dt%d_photTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%d%s%s%s%s",ipixelany,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,iitss,imult,itofpos,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,cweightsback.Data(),cmvx.Data(),cbeauty.Data(),kfp.Data()));
  
  printf("Add macro appendix %s\n", appendix.Data());
  
  if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors")) gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors.C");
  if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpePbPb"))gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFEnpePbPb.C");
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEnpePbPb(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itss, itshitpixel, itschi2percluster, tpcsharedcluster, etacorr, multicorr, toflast, etaIncMin, etaIncMax,
					      assETAm, assETAp, assITS, assTPCcl, assTPCPIDcl,
					      assDCAr, assDCAz, assTPCSminus, assTPCSplus,
					      useCat1Tracks, useCat2Tracks, weightlevelback,usekfparticle);
  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();
  
  if (useMC)	task->SetHasMCData(kTRUE);
  else		task->SetHasMCData(kFALSE);

  if(useMC&&(beauty || (weightlevelback>=0))) ConfigWeightFactors(task,kFALSE,wei);//2;For default PbPb

  // ----- centrality selection -----
  task->SetCentralityCheck(newCentralitySelection,"V0M");
  task->SetCentralityInterval(centrMin,centrMax);               // all events outside the desired centrality interval are rejected
  // --------------------------------
  // ----- trigger selecton ---------
  if(!newCentralitySelection)   task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral); // old framework
  if(newCentralitySelection)    task->SelectCollisionCandidates(AliVEvent::kINT7);                                               // new framework
  // --------------------------------
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFEtask";
  containerName += appendix.Data();
  printf("container name: %s\n", containerName.Data());

  //create data containers
  task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
  task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput(task,  0, cinput );
  
  mgr->AddTask(task);

  return NULL;
}
