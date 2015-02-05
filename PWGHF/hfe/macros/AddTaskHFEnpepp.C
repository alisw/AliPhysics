AliAnalysisTask *AddTaskHFEnpepp(Bool_t MCthere,
				 Bool_t isAOD, 
                                 Bool_t kNPERef = kTRUE,
				 Bool_t kNPERefTPConly = kFALSE
){

 // Default settings (TOF-TPC pPb)
  const int	kDefTPCcl	= 120;
  const int	kDefTPCclPID	=  80;
  const int	kDefITScl	=   4;
  const double	kDefDCAr	=   1.;
  const double	kDefDCAz	=   2.;
  const double	kDefTOFs	=   3.;
  const double  kDefEtaIncMin   = -0.8;
  const double  kDefEtaIncMax   =  0.8;
  const double  kDefPhiMin      = -1.; // by default no phi cut, otherwise units of 2.*TMath::Pi()/18.
  const double  kDefPhiMax      = -1.; // by default no phi cut
 

  // TPC PID Cuts Inclusive leg:
  // General, if mean=0 and sigma=1: 
  // Above 3 sigma we neglect 0.13%
  // cut in sigma (effective efficiency from cut to 3 sigma)
  // -1 (84%), -0.75 (77.2%), -0.5 (69%), -0.25 (59.7%), -0.129 (55%)
  //  0 (49.9%), 0.122 (45%), 0.25 (40%), 0.5 (30.7%)

  // ESDs: mean 0.06, sigma 1    --> -0.94, -0.69, -0.44,  -0.19, -0.009 ,0.06, 0.182, 0.31, 0.56
  // AODs: mean 0.09, sigma 1    --> -0.91, -0.66, -0.41,  -0.16, -0.039 ,0.09, 0.212, 0.34, 0.59
  // AODs: mean 0.09, sigma 1.03 --> -0.94, -0.68, -0.425, -0.17, -0.043 ,0.09, 0.216, 0.35, 0.60

  // On ESD:
  // mean is actually 0.06 (abs(eta)<0.6)
  Double_t dEdxhm[12] = {3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2};  // Above 3 sigma we neglect 0.13%
  Double_t tpcl0[12]  = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};  // 50%
  Double_t tpcl1[12]  = {0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};  // 48.34%
  Double_t tpcl2[12]  = {0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08};  // 46.74%
  Double_t tpcl3[12]  = {0.122,0.122,0.122,0.122,0.122,0.122,0.122,0.122,0.122,0.122,0.122,0.122};  // 45%
  Double_t tpcl4[12]  = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2};  // 42%
  Double_t tpcl5[12]  = {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};  // 40%
  Double_t tpcl6[12]  = {-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04};  // 51.5%
  Double_t tpcl7[12]  = {-0.08,-0.08,-0.08,-0.08,-0.08,-0.08,-0.08,-0.08,-0.08,-0.08,-0.08,-0.08};  // 53.1%
  Double_t tpcl8[12]  = {-0.122,-0.122,-0.122,-0.122,-0.122,-0.122,-0.122,-0.122,-0.122,-0.122,-0.122,-0.122};  // 54.8%
  Double_t tpcl9[12]  = {-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2};  // 57.9%
  Double_t tpcl10[12]  = {-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.25};  // 59.8%


  /*
  Double_t dEdxhm[12] = {3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06};  // Above 3 sigma we neglect 0.13%
  Double_t tpcl0[12]  = {-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94};  // 84%
  Double_t tpcl1[12]  = {-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44};  // 69%
  Double_t tpcl2[12]  = {0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06};  // 50%
  Double_t tpcl3[12]  = {-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69};  // 77.2%
  Double_t tpcl4[12]  = {-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19};  // 59.7%
  Double_t tpcl5[12]  = {0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186};  // 45%
  Double_t tpcl6[12]  = {0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31};  // 40%
  Double_t tpcl7[12]  = {0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56};  // 30.7%
  */

  // Default setting for the associated electron for the NonPhotonic Analysis
  const double	kassETAm        = -0.8;
  const double	kassETAp        =  0.8;
  const double  kassMinPt       =  0.1;
  const int	kassITS		=    2;
  const int	kassTPCcl	=   60;
  const int	kassTPCPIDcl	=   60;
  const double	kassDCAr	=  1.0;
  const double	kassDCAz	=  2.0;
  const double	kassTPCSminus	= -3.0;
  const double	kassTPCSplus	=  3.0;
  const double  kassITSpid      =  3.0;
  const double  kassTOFpid      =  0.0;
  
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_hfe_HFE", "No analysis manager found.");
    return 0;
  }

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //@@ 0 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Double_t dEdxaclm[12], dEdxachm[12],dEdxaclm1[12], dEdxachm1[12],dEdxaclm2[12], dEdxachm2[12];
  for(int icent = 0; icent < 12; icent++){
    dEdxaclm[icent] = kassTPCSminus;
    dEdxachm[icent] = kassTPCSplus;
    dEdxaclm1[icent] = -2;
    dEdxachm1[icent] = 2;
    dEdxaclm2[icent] = -4;
    dEdxachm2[icent] = 4;
  }

  Int_t kWei = -1;
  if (MCthere) kWei = 0;
  enum {
    kWeiLHC10f6a = 10,  // weights for published pp @ 7 TeV: LHC10f6a, Pythia low-pt (Perugia0)
    kWeiLHC10f6 = 11,  // weights for published pp @ 7 TeV: LHC10f6, Pythia low-pt (Perugia0)
    kWeiLHC10f7a_d = 12,  // weights for published pp @ 7 TeV: LHC10f7a_d, Pythia low-pt (Perugia0)
  };
  int kWeiData; 
  kWeiData = kWeiLHC10f6a;
  /*
  TString list = gSystem->Getenv("LIST");
  if(list.Contains("LHC10f6a")){
	kWeiData = kWeiLHC10f6a;
  } else if (list.Contains("LHC10f6")){
	kWeiData = kWeiLHC10f6;
  } else if (list.Contains("LHC10f7a_d")){
	kWeiData = kWeiLHC10f7a_d;
  } else {
	kWei = -1;
  	printf("no weighting found");
  } 
  */

  if(kNPERef){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
    if(isAOD==1){ 
      // Reference
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
    }
    else {
      // Reference
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
     
      // different re-weighting (only for MC)
      if (MCthere){
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f7a_d);
      }
    }
  }
  
  if(kNPERefTPConly){
    // **************************************************************
    // 
    // Reference task for TPC-only on the inclusive leg
    //
    // **************************************************************
    if(isAOD==1){
      // Reference
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, 0., AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
    }
    else
      { 
	// Reference
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, 0., AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
	
	// different re-weighting (only for MC)
	if (MCthere){
	  RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			     dEdxhm, 0., AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			     kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6);
	  RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			     dEdxhm, 0., AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			     kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f7a_d);
	}
	
      }
  }
  
  
  return NULL;
}

//===============================================================================
AliAnalysisTask *RegisterTaskNPEpp(Bool_t useMC, Bool_t isAOD,
				   Int_t tpcCls=120, Int_t tpcClsPID=80, 
				   Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
				   Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL, 
				   Double_t tofs=3., Int_t itshitpixel =AliHFEextraCuts::kBoth, 
				   Int_t iKink = 0,
				   Int_t assITS=2, Int_t assTPCcl=100,
				   Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
				   Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
				   Double_t assITSpid = 3.0, Double_t assTOFpid = 0.0,
				   Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE,
				   Int_t weightlevelback = -1, Int_t WhichWei = 0)
{
  // Fixed values
  Double_t etaIncMin = -0.8; Double_t etaIncMax = 0.8;
  Double_t phimi = -1.; Double_t phima = -1.;
  Double_t assETAm=-0.8; Double_t assETAp=0.8;
  Double_t assMinPt = 0.1;

  //
  // Cuts on the inclusive leg
  //
  Int_t idcaxy = (Int_t)(dcaxy*10.);
  Int_t idcaz = (Int_t)(dcaz*10.);
  Int_t tpclow = 0;
  if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t ipixelany = itshitpixel;
  TString phirange("");
  if (phimi >= 0. && phima >= 0.){ 
    phirange += "Phi";
    phirange += phimi;
    phirange += "-";
    phirange += phima;
  } 

  //
  // Cuts on the associated leg
  //
  Int_t iassDCAr = (Int_t)(assDCAr*10);
  Int_t iassDCAz = (Int_t)(assDCAz*10);
  Int_t iassTPCSminus = assTPCSplus ? (Int_t)(assTPCSplus[0]*1000.) : 0;
  Int_t iassTOF = (Int_t)(assTOFpid*10);
  Int_t iassITS = (Int_t)(assITSpid * 10.);
  Int_t phoTrack = 0;
  if (useCat1Tracks) phoTrack = 1;
  if (useCat2Tracks) phoTrack = 2;

  TString cweightsback("");
  if(weightlevelback>=0) {
    cweightsback += "Wa";
    if (WhichWei>0){
      cweightsback += WhichWei;
      //cweightsback += weightlevelback;
    }
  }

  TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCs%dTOFs%dK%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%d%s",
				   tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,itofs,iKink,assTPCcl,assTPCPIDcl,assITS,iassDCAr,
				   iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,cweightsback.Data()));

  printf("Add macro appendix %s\n", appendix.Data());

  if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors")) 
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigWeightFactors.C");
  if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpepp")) 
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHFEnpepp.C");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  //mgr->AddClassDebug("AliHFENonPhotonicElectron", 1);
  AliAnalysisTaskHFE *task = ConfigHFEnpepp(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					    tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, iKink,etaIncMin, etaIncMax,
					    phimi*TMath::Pi()/9., phima*TMath::Pi()/9., 
					    assETAm, assETAp, assMinPt, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, 
					    assTPCSplus,assITSpid,assTOFpid,
					    useCat1Tracks, useCat2Tracks, weightlevelback);

  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();
  
  if (useMC)	task->SetHasMCData(kTRUE);
  else		task->SetHasMCData(kFALSE);

  task->SelectCollisionCandidates();

  if(useMC && weightlevelback>=0) {
    ConfigWeightFactors(task,kFALSE,WhichWei);
  }

  TString containerName = mgr->GetCommonFileName();
  containerName += ":HFEtask";
  containerName += appendix.Data();
  printf("container name: %s\n", containerName.Data());

  //create data containers
  task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(), 
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(), 
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput(task,  0, cinput );

  mgr->AddTask(task);

  return NULL;
}

