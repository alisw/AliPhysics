AliAnalysisTask *AddTaskHFEnpepPb(Bool_t MCthere,
                                  Bool_t isAOD, 
                                  Bool_t kNPERef = kTRUE,
                                  Bool_t kNPERefTPConly = kFALSE,
                                  Bool_t kNPEkAny = kFALSE,
                                  Bool_t kNPEsystematics = kFALSE,
                                  Bool_t kNPEpidsys = kFALSE,
				  Bool_t kNPEw = kFALSE,
                                  Bool_t kBeauty = kFALSE
){

  // Default settings (TOF-TPC pPb)
  const int	kDefTPCcl	= 100;
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
  Double_t dEdxhm[12] = {3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06,3.06};  // Above 3 sigma we neglect 0.13%
  Double_t tpcl0[12]  = {-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94};  // 84%
  Double_t tpcl1[12]  = {-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44,-0.44};  // 69%
  Double_t tpcl2[12]  = {0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06};  // 50%
  Double_t tpcl3[12]  = {-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69};  // 77.2%
  Double_t tpcl4[12]  = {-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19};  // 59.7%
  Double_t tpcl5[12]  = {0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186,0.186};  // 45%
  Double_t tpcl6[12]  = {0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31,0.31};  // 40%
  Double_t tpcl7[12]  = {0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56,0.56};  // 30.7%

  // AOD 139 with TPC multi cor; w/o TOF PID
  Double_t dEdxhmAOD[12] = {3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18};  // upper cut
  Double_t tpcl0AOD[12]  = {-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94};  // 84% <-- hadron contamination
  Double_t tpcl1AOD[12]  = {-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68};  // 77.2%
  Double_t tpcl2AOD[12]  = {-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425};  // 69% <-- had cont
  Double_t tpcl3AOD[12]  = {-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17};  // 59.7%           <-- had cont
  Double_t tpcl4AOD[12]  = {-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043};  // 55%
  Double_t tpcl5AOD[12]  = {0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09};  // 49.9%                       <-- had cont
  Double_t tpcl6AOD[12]  = {0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216};  // 45%
  Double_t tpcl7AOD[12]  = {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35};  // 40%
  Double_t tpcl8AOD[12]  = {0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60};  // 30.7%

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

  enum{
          kHFEV0A = 1,
          kHFEV0M = 2,
          kHFECL1 = 3,
          kHFEZNA = 4
  };

  const Bool_t isBeauty = kFALSE;
  
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
    // updated after changes by Jan on 17.04, commit 21648
    kb2weiData = 3,  // DPMJET weights 
    kd3weiData = 4,  // HIJING weights 
    kd3weib2   = 5,  // HIJING weights for DPMJET
  };
  
  

  if(kNPERef){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
    if(isAOD==1){ 
      // Reference
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      /*
      if (MCthere){  // run without weights
	// this is an exceptional case, so I handle it by hand for the moment
	RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			    dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			    kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData);
	RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			    dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			    kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weib2);
	RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			    dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			    kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE);
      }
      // no ITS SA
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */
    }
    else {
      // Reference
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      }
  }
  
  if(kNPERefTPConly){
    // **************************************************************
    // 
    // Reference task for TPC-only on the inclusive leg
    //
    // **************************************************************
    if(isAOD==1){
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE);
    }
    else
      {
	RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, 
			    dEdxhm, 0., AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			    kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE);
      }
  }
  
  if(kNPEkAny){
    // **************************************************************
    // 
    // task for kAny instead of kBoth
    //
    // **************************************************************

    // coming soon
  }
  
  if (kNPEpidsys){
    if(isAOD==1){
      // Inclusive: TPC pid cuts
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      // Inclusive: TOF PID cuts
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 2.0, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 2.5, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 4.0, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
    }
    else
      {
	RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5, 
			    dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			    kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE);
      }

  }

  if(kNPEsystematics){
    // **************************************************************
    // 
    // Cut systematics 
    //
    // **************************************************************
    if(isAOD==1){

      // TPC clusters for tracking
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 80, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 90, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty,100, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty,120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty,130, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      // TPC clusters for PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, 60, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl,100, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl,110, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);

      /* SYST 2
      // Cat1 + cat2 tracks: vary ITS PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, 2.0, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, 4.0, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      // Cat1 + cat2 tracks: vary TPC PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm1, dEdxachm1, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm2, dEdxachm2, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      // no ITS SA: add TOF PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, 2.0, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, 3.0, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, 4.0, kTRUE, kFALSE, kWei, kd3weiData);
      // no ITS SA: add TOF PID, vary TPC PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm1, dEdxachm1, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm2, dEdxachm2, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */
      /* SYST 1
      // Inclusive DCA cuts
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, 0.5, 1.0, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, 2.0, 5.0, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      // Associated DCA cuts
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  0.5, 1.0, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  2.0, 5.0, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei);
      */
    }
  }

  if(kBeauty){
    // **************************************************************
    //
    // Beauty task
    //
    // **************************************************************
    if(isAOD==1){
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, dEdxhmAOD, kDefTOFs,
                          AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl,
                          kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE, -1, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, dEdxhmAOD, kDefTOFs,
                          AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl,
                          kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE,-1, kTRUE, kTRUE, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, dEdxhmAOD, kDefTOFs,
                          AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl,
                          kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE,-1, kTRUE, kTRUE, kTRUE);
    }
    else{
      // reference
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      // nonphotonic bg 
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiDataSys, kTRUE); 
      //ip sys
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE, kFALSE, kFALSE, 1);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE, kFALSE, kFALSE, 2);

      // tpc cluster
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 100, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);

      // tpc pid cluster
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, 60, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, 100, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);

      // dca xy
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, 2, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, 0.5, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, 0.1, 0.5, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);

      // dca z 
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, 4, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, 1, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);

      // ITS cluster
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, 2, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, 3, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);

      // TOF PID
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, 2, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1,
                          dEdxhm, 4, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);

      // TPC PID
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0, //84%
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, //50%
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData, kTRUE);
    }
  }
  
  return NULL;
}

//===============================================================================
AliAnalysisTask *RegisterTaskNPEpPb(Bool_t useMC, Bool_t isAOD, Bool_t beauty,
               Int_t tpcCls=120, Int_t tpcClsPID=80, 
               Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
               Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL, 
               Double_t tofs=3., Int_t itshitpixel =AliHFEextraCuts::kBoth, 
				    //Double_t phimi=-1., Double_t phima=-1.,
	       Int_t icent=1,
	       Int_t assITS=2, Int_t assTPCcl=100,
               Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
               Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
               Double_t assITSpid = 3.0, Double_t assTOFpid = 0.0,
               Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE,
	       Int_t weightlevelback = -1, Int_t WhichWei = 0,
               Bool_t npeBeauty = kFALSE, Bool_t ipCharge = kFALSE, Bool_t isOpp = kFALSE, Int_t ipSys = 0)
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

  //printf("Argument passed to function to determine the centrality estimator %i\n", icent);
  if (icent == 2) TString cesti("V0M");
  else if (icent == 3) TString cesti("CL1");
  else if (icent == 4) TString cesti("ZNA");
  else TString cesti("V0A");
  //printf("Centrality estimator %s\n", cesti.Data());

  TString cweightsback("");
  if(weightlevelback>=0) {
    cweightsback += "Wa";
    if (WhichWei>0){
      cweightsback += WhichWei;
      //cweightsback += weightlevelback;
    }
  }

  if(beauty) {
     if(ipCharge && isOpp) TString cbeauty("BeautyIPopp");
     else if(ipCharge) TString cbeauty("BeautyIPcrg");
     else if(!ipCharge) TString cbeauty("Beauty");
     else TString cbeauty("BeautyWrong");
     cbeauty = cbeauty+Form("ip%dwei%d",ipSys,WhichWei);
  }
  else TString cbeauty("");

  TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCs%dTOFs%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%dce%s%s%s",
				   tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,itofs,assTPCcl,assTPCPIDcl,assITS,iassDCAr,
				   iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,cesti.Data(),cweightsback.Data(),cbeauty.Data()));

  printf("Add macro appendix %s\n", appendix.Data());

 if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors")) 
     gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigWeightFactors.C");
 if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpepPb")) 
     gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEnpepPb.C");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEnpepPb(useMC, isAOD, beauty, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					     tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, icent, etaIncMin, etaIncMax,
					     phimi*TMath::Pi()/9., phima*TMath::Pi()/9., 
					     assETAm, assETAp, assMinPt, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, 
					     assTPCSplus,assITSpid,assTOFpid,
					     useCat1Tracks, useCat2Tracks, weightlevelback, npeBeauty, ipCharge, isOpp, ipSys);

  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();

  if (useMC)	task->SetHasMCData(kTRUE);
  else		task->SetHasMCData(kFALSE);

  task->SelectCollisionCandidates(AliVEvent::kINT7);

  if(useMC&&(beauty || (weightlevelback>=0))) {
    ConfigWeightFactors(task,kFALSE,WhichWei);
    if(WhichWei==8) {
      ConfigWeightFactors(task,kTRUE,3);
      task->SetNonHFEsystematics(kTRUE);
    }
  }
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
