AliAnalysisTask *AddTaskHFEnpepp(Bool_t MCthere,
				 Bool_t isAOD, 
                                 Bool_t kNPERef = kTRUE,
				 Bool_t kNPERefTPConly = kFALSE
){

  // Default settings (TOF-TPC pp)
  // ESD analysis of LHC10d pass2
 
  const int	kDefTPCcl	= 120;
  const int	kDefTPCclPID	=  80;
  const int	kDefITScl	=   4;
  const double	kDefDCAr	=   1.;
  const double	kDefDCAz	=   2.;
  const double	kDefTOFs	=   3.;
  
  // TPC PID Cuts Inclusive leg:
  // General, if mean=0 and sigma=1: 
  // Above 3 sigma we neglect 0.13%.
  // Cut in sigma (effective efficiency from cut to 3 sigma)
  // -1 (84%), -0.75 (77.2%), -0.5 (69%), -0.25 (59.7%), -0.129 (55%)
  //  0 (49.9%), 0.122 (45%), 0.25 (40%), 0.5 (30.7%)

  // On ESD:
  Double_t dEdxhm[12] = {3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2};  // Above 3 sigma we neglect 0.13%
  Double_t dEdxhm20[12] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0};
  Double_t dEdxhm21[12] = {1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};
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
  Double_t tpcl11[12] = {-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2}; 

  Double_t tpcl20[12]  = {0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002};  // 
  Double_t tpcl21[12]  = {0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};  // 


  // Default setting for the associated electron for the NonPhotonic Analysis
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
    kWeiLHC10f6 = 10,      // weights for published pp @ 7 TeV: LHC10f6, Phojet
    kWeiLHC10f6a = 11,     // weights for published pp @ 7 TeV: LHC10f6a, Pythia
    kWeiLHC10f7a_d = 12,  // weights for published pp @ 7 TeV: LHC10f7a_d, Pythia hf enhanced
    kWeiLHC10f6af7a_d = 13,   // new for MC clouser test
    kWeiLHC10f7a_df6a = 14,   
    kWeiLHC10f6f6a = 15,
    kWeiLHC10f6af6 = 16,   
    kWeiLHC10f6f7a_d = 17,
    kWeiLHC10f7a_df6 = 18
  };
  int kWeiData; 
  //kWeiData = kWeiLHC10f6a; changed on September 15, again aware about re-weighting
  kWeiData = kWeiLHC10f7a_d;

 
  if(kNPERef){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
    if(isAOD==1){ 
      // AOD analysis currently not available
    }
    else {
      // Reference
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      /*
      // A few systematic checks
      // include kink mothers
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 1, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      // SPD selection
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kAny, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      // ITS hits
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, 3, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      */
      /*
      // TPC clusters
      RegisterTaskNPEpp( MCthere, isAOD, 90, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      RegisterTaskNPEpp( MCthere, isAOD, 100, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      RegisterTaskNPEpp( MCthere, isAOD, 110, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      RegisterTaskNPEpp( MCthere, isAOD, 130, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      // with 130 TPC clusters, vary also TPC PID clusters
      RegisterTaskNPEpp( MCthere, isAOD, 130, 100, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      */
      // different re-weighting (only for MC)
      if (MCthere){
	/*
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6);
	*/

	// with no weights
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);
	// different weights
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6a);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6af7a_d);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f7a_df6a);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6f6a);	
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6af6);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f6f7a_d);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC10f7a_df6);

    }
      // different TOF PID cuts
      /*
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, 2., AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, 3.5, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, 4., AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			 dEdxhm, 5., AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      */
      /*
      // different TPC PID cuts (only for data)
      if (!MCthere){
	// vary top dE/dx cut
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl20, 
			   dEdxhm20, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl21, 
			   dEdxhm21, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl8, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl10, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl11, 
			   dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
      }
      */
    }
  }
  
  if(kNPERefTPConly){
    // **************************************************************
    // 
    // Reference task for TPC-only on the inclusive leg
    //
    // **************************************************************
    if(isAOD==1){
      // AOD analysis currently not available
    }
    else
      { 
	RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, 
			   dEdxhm, 0., AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl, 
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);	
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
	       Int_t iKink = 0, Int_t assITS=2, Int_t assTPCcl=100,
				   Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
				   Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
				   Double_t assITSpid = 3.0, Double_t assTOFpid = 0.0,
				   Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE,
				   Int_t weightlevelback = -1, Int_t WhichWei = 0)
{
  // Fixed values
  Double_t etaIncMin = -0.8; Double_t etaIncMax = 0.8;
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
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors.C");
  if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpepp")) 
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHFEnpepp.C");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  //mgr->AddClassDebug("AliHFENonPhotonicElectron", 1);
  AliAnalysisTaskHFE *task = ConfigHFEnpepp(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					    tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, iKink,etaIncMin, etaIncMax,
					    assETAm, assETAp, assMinPt, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, 
					    assTPCSplus,assITSpid,assTOFpid, useCat1Tracks, useCat2Tracks, weightlevelback);

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

