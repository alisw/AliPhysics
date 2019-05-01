AliAnalysisTask *AddTaskHFEnpepPb(Bool_t MCthere,
                                  Bool_t isAOD, 
                                  Bool_t kNPERef = kTRUE,
                                  Bool_t kNPERefTPConly = kFALSE,
                                  Bool_t kNPEkAny = kFALSE,
                                  Bool_t kNPEkCentrality = kFALSE,
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
  Double_t dEdxhmAOD[12] = {3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18,3.18};  // upper cut, 3 sigma
  Double_t dEdxhm1AOD[12] = {2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15};  // upper cut, 2 sigma
  Double_t dEdxhm2AOD[12] = {2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67,2.67};  // upper cut, 2.5 sigma
  Double_t dEdxhm3AOD[12] = {4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21,4.21};  // upper cut, 4 sigma
  Double_t tpcl0AOD[12]  = {-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94};  // 84% <-- hadron contamination
  Double_t tpcl1AOD[12]  = {-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68};  // 77.2%
  Double_t tpcl2AOD[12]  = {-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425};  // 69% <-- had cont
  Double_t tpcl3AOD[12]  = {-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17};  // 59.7%           <-- had cont
  Double_t tpcl4AOD[12]  = {-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043};  // 55%
  Double_t tpcl5AOD[12]  = {0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09};  // 49.9%                       <-- had cont
  Double_t tpcl6AOD[12]  = {0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216};  // 45%
  Double_t tpcl7AOD[12]  = {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35};  // 40%
  Double_t tpcl8AOD[12]  = {0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60};  // 30.7%

  // IP systematics
  Double_t ipParup0[3] = {0.006294,0.082263,-0.568599}; // +0.5 sigma from the reference 
  Double_t ipParup1[3] = {0.006665,0.083763,-0.571077}; // +0.7
  Double_t ipParup2[3] = {0.007214,0.086286,-0.575887}; // +1              ---> sys upper
  Double_t ipParup3[3] = {0.008144,0.090745,-0.585015}; // +1.5
  Double_t ipParup4[3] = {0.009010,0.094635,-0.590015}; // +2  
  Double_t ipParup5[3] = {0.010767,0.103133,-0.601536}; // +3 
  Double_t ipPardn0[3] = {0.005108,0.075929,-0.564700}; // -0.3
  Double_t ipPardn1[3] = {0.004505,0.073860,-0.551137}; // -0.5
  Double_t ipPardn2[3] = {0.004140,0.072280,-0.547989}; // -0.7
  Double_t ipPardn3[3] = {0.003617,0.069823,-0.542757}; // -1              ---> sys lower
  Double_t ipPardn4[3] = {0.002731,0.065968,-0.535332}; // -1.5
  Double_t ipPardn5[3] = {0.001809,0.061884,-0.523571}; // -2 
  Double_t ipPardn6[3] = {0.000905,0.057915,-0.512314}; // -2.5
  Double_t ipPardn7[3] = {0.000064,0.054140,-0.503857}; // -3

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
  Double_t EtaDalWei = 0.;
  

  if(kNPERef){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************

    // July 17: add by hand the variable to check the SPD status or not, only in few cases, to be removed!!!!
    // July 21: change reference to --> no ITS SA included!

    if(isAOD==1){ 
      // Old reference with ITS SA included
      //RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
      //		  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
      //		  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      
      // New reference (21.07) without ITS SA
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);

      if (MCthere){
      /*
	// Checks around the re-weighting, now done without ITS SA !!! (21.07)
	// vary the weight of the eta Dalitz contribution
	RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
                            dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
                            kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData, 1.2);
	RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
                            dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
                            kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData, 0.8);
	// ... with no weights
	RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			    dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			    kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);
      */
 	// this is an exceptional case, so I handle it by hand for the moment
	// It is to run with special weights (other than the default d3 to describe data, or ....
	//RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
	//		    dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
	//		    kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kb2weiData);
	//RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
	//		    dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
	//		    kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weib2);
      }
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
    // From July 21, without ITS SA
    // **************************************************************
    if(isAOD==1){
      // Reference: 50% PID cut. For the moment, old contamination file
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      /*// Reference: 50% PID cut. Here use the new contamination file
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 1, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // Reference: 50% PID cut. Here use the new contamination file
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 2, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // Reference: 50% PID cut. Here use the new contamination file
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 3, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */
      /* 1
      // SYSTEMATICS TPC CLUSTERS
      // for tracking
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 90, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 95, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 105, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 110, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 115, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // for PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, 60, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, 100, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */

      /* 2
      //SYSTEMATICS TPC PID
      // with the old contamination file
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // with the new contamination file
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 1, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 1, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 1, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl4AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 1, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 1, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl6AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 1, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */

      /*
      // SYSTEMATICS ITS
      // number of ITS hits
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 3, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 6, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // kAny, and vary the number of ITS hits
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kAny, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 2, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kAny, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 3, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kAny, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kAny, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // kFirst, and vary the number of ITS hits
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kFirst, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kFirst, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */
      /*
      // SYSTEMATICS ASSOCIATED
      // ITS PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, 2.0, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, 4.0, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // TPC PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm1, dEdxachm1, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm2, dEdxachm2, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // add TOF PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, 3.0, kTRUE, kFALSE, kWei, kd3weiData);
      // ITS number of hits
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, 3.0, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, 4.0, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, 5.0, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // TPC tracking clusters
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, 80, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, 100, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // TPC PID Clusters
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, 80, 70, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, 0., AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, 100, 80, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */
      }
  }
  
  if(kNPEkAny){
    // **************************************************************
    // 
    // task for kAny instead of kBoth, then vary the cut on ITS hits
    //
    // **************************************************************
    if(isAOD==1){
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kAny, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 2, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kAny, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 3, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kAny, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kAny, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // shortly add kFirst
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 4, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kFirst, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kFirst, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
    }
  }

  if(kNPEkCentrality){
    // **************************************************************
    // 
    // task for ZNA instead of V0A
    //
    // **************************************************************
    if(isAOD==1){
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEZNA, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
    }
  }
  
  if (kNPEpidsys){
    if(isAOD==1){
      // Inclusive: TPC pid cuts
      // vary lower cut
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // vary upper cut
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhm1AOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhm2AOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhm3AOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // Inclusive: TOF PID cuts
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 1.5, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 2.0, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 2.5, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, 4.0, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
    }
  }

  if(kNPEsystematics){
    // **************************************************************
    // 
    // Cut systematics 
    //
    // **************************************************************
    if(isAOD==1){

      /*
      // SYST 1
      // TPC clusters for tracking
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 90, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 95, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 105, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 110, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, 115, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty,120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // old checks
      //RegisterTaskNPEpPb( MCthere, isAOD, isBeauty,130, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
      //		  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
      //		  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);
      //RegisterTaskNPEpPb( MCthere, isAOD, isBeauty,140, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
      //		  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
      //		  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData);

      // TPC clusters for PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, 60, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl,100, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */
      /*
      // SYST 2
      // Variation of cuts on the associated electrons
      // Vary associated ITS PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, 2.0, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, 4.0, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // Vary associated TPC PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm1, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm2, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // add TOF PID
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, 3.0, kTRUE, kFALSE, kWei, kd3weiData);
      // vary cut on ITS hits
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, 3, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, 4, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, 5, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // Associated: vary cut on TPC tracking clusters
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, 80, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, 100, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // Associated: vary cut on TPC PID clusters
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, 80, 70, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, 100, 80, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */

      // SYST 3
      // Inclusive DCA cuts
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, 0.5, 1.0, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, 2.0, 5.0, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      // Associated DCA cuts
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  0.5, 1.0, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  2.0, 5.0, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      /*
      // Inclusive ITS hits
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 3, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      RegisterTaskNPEpPb( MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, 6, kDefDCAr, kDefDCAz, tpcl2AOD, 
			  dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl, 
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kd3weiData);
      */
    }
  }

  if(kBeauty){
    // **************************************************************
    //
    // Beauty task
    //
    // **************************************************************
    const int	kDefTPCclBeauty	= 100;  // before 110
    if(isAOD==1){
	RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
			   dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			   kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE, kTRUE);
        // last argument for  nonHFE systematics' container kTRUE = filled
    }
    else{
      // reference
      /*RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE, kFALSE);*/

      // reference : use this task if you want to run reference together with nonHFE systematics' container filled
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE, kTRUE);

      // modifying TPC cluster and TPC cluster at the same time
      /*RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 100, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);

      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 120, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);*/


      // hadron cont
/*
       RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 1, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
       RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 2, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
       RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 3, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
*/
/*
      //ip sys : done  1st train
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE, kFALSE, kFALSE, kFALSE, ipParup2);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE, kFALSE, kFALSE, kFALSE, ipPardn3); */
/*
      // tpc cluster : done 2nd train
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 90, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 100, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 105, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 115, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, 120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
*/

      // tpc pid cluster : done 3rd train
/*
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, 60, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
			  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, 100, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);

			  // dca xy : done 3rd train
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, 2, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, 0.5, kDefDCAz, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, 0.1, 0.5, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
*/
  
      // dca z : done 2nd train
/*      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, 4, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
      RegisterTaskNPEpPb( MCthere, isAOD, kBeauty, kDefTPCclBeauty, kDefTPCclPID, kDefITScl, kDefDCAr, 1, tpcl2AOD,
                          dEdxhmAOD, kDefTOFs, AliHFEextraCuts::kBoth, 0, kHFEV0A, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kd3weiData, 0, kFALSE);
			  */




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
	       Int_t SPDcheck = 0,
				    //Double_t phimi=-1., Double_t phima=-1.,
	       Int_t icent=1,
	       Int_t assITS=2, Int_t assTPCcl=100,
               Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
               Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
               Double_t assITSpid = 3.0, Double_t assTOFpid = 0.0,
               Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE,
	       Int_t weightlevelback = -1, Int_t WhichWei = 0, Double_t etadalwei = 0.,
               Bool_t npeBeauty = kFALSE, Bool_t nonHFEsys = kFALSE, Bool_t ipCharge = kFALSE, Bool_t isOpp = kFALSE, Double_t *ipPar=NULL)
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
  Int_t tpchi = 0;
  if(tpcdEdxcuthigh) tpchi = (Int_t) (tpcdEdxcuthigh[0]*1000.);
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
  Int_t ipSys = 0;
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
    Int_t wei_index = WhichWei + 10*weightlevelback;
    cweightsback += "Wa";
    if (WhichWei>0){
      //cweightsback += WhichWei;
      cweightsback += wei_index;
      if (etadalwei>0){
	cweightsback += etadalwei*10.;
      }
    }
  }

  if(beauty) {
     if(ipCharge && isOpp) TString cbeauty("BeautyIPopp");
     else if(ipCharge) TString cbeauty("BeautyIPcrg");
     else if(!ipCharge) TString cbeauty("Beauty");
     else TString cbeauty("BeautyWrong");
     if(ipPar) ipSys = (Int_t) (ipPar[0]*100000.);
     cbeauty = cbeauty+Form("ip%d",ipSys);
     if(nonHFEsys) cbeauty = cbeauty+"Syst";
  }
  else TString cbeauty("");

  TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%d%dDCAr%dz%dTPCs%dh%dTOFs%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%dce%s%s%s",
				   tpcCls,tpcClsPID,itsCls,ipixelany,SPDcheck,
				   idcaxy,idcaz,tpclow,tpchi,itofs,assTPCcl,assTPCPIDcl,assITS,iassDCAr,
				   iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,cesti.Data(),cweightsback.Data(),cbeauty.Data()));

  printf("Add macro appendix %s\n", appendix.Data());

 if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors")) 
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors.C");
 if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpepPb")) 
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigHFEnpepPb.C");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEnpepPb(useMC, isAOD, beauty, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
					     tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, SPDcheck, icent, etaIncMin, etaIncMax,
					     phimi*TMath::Pi()/9., phima*TMath::Pi()/9., 
					     assETAm, assETAp, assMinPt, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, 
					     assTPCSplus,assITSpid,assTOFpid,
					     useCat1Tracks, useCat2Tracks, weightlevelback, etadalwei,
					     npeBeauty, ipCharge, isOpp, ipPar);

  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();

  if (useMC)	task->SetHasMCData(kTRUE);
  else		task->SetHasMCData(kFALSE);

  task->SelectCollisionCandidates(AliVEvent::kINT7);

  if(useMC&&(beauty || (weightlevelback>=0))) {
    ConfigWeightFactors(task,nonHFEsys,WhichWei);
    task->SetNonHFEsystematics(nonHFEsys);
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


