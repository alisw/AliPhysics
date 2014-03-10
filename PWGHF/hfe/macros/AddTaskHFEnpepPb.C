AliAnalysisTask *AddTaskHFEnpepPb(Bool_t MCthere,
                                  Bool_t isAOD, 
                                  Bool_t kNPERef = kTRUE,
                                  Bool_t kNPERefTPConly = kTRUE,
                                  Bool_t kNPEkAny = kFALSE,
                                  Bool_t kNPEsystematics = kFALSE,
                                  Bool_t kNPEpidsys = kFALSE
){

  // Default settings (TOF-TPC pPb)
  const int	kDefTPCcl	= 110;
  const int	kDefTPCclPID	=  80;
  const int	kDefITScl	=   4;
  const double	kDefDCAr	=   1.;
  const double	kDefDCAz	=   2.;
  const double	kDefTOFs	=   3.;
  const double  kDefEtaIncMin = -0.8;
  const double  kDefEtaIncMax = 0.8;


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
  Double_t dEdxhmAOD[12] = {3.13,3.13,3.13,3.13,3.13,3.13,3.13,3.13,3.13,3.13,3.13,3.13};  // upper cut
  Double_t tpcl0AOD[12]  = {-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94,-0.94};  // 84%
  Double_t tpcl1AOD[12]  = {-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68,-0.68};  // 77.2%
  Double_t tpcl2AOD[12]  = {-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425,-0.425};  // 69%
  Double_t tpcl3AOD[12]  = {-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17};  // 59.7%
  Double_t tpcl4AOD[12]  = {-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043,-0.043};  // 55%
  Double_t tpcl5AOD[12]  = {0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09};  // 49.9%
  Double_t tpcl6AOD[12]  = {0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216,0.216};  // 45%
  Double_t tpcl7AOD[12]  = {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35};  // 40%
  Double_t tpcl8AOD[12]  = {0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60,0.60};  // 30.7%

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

  Double_t dEdxaclm[12], dEdxachm[12],dEdxaclm1[12], dEdxachm1[12],dEdxaclm2[12], dEdxachm2[12];
  for(int icent = 0; icent < 12; icent++){
    dEdxaclm[icent] = kassTPCSminus;
    dEdxachm[icent] = kassTPCSplus;
    dEdxaclm1[icent] = -2;
    dEdxachm1[icent] = 2;
    dEdxaclm2[icent] = -4;
    dEdxachm2[icent] = 4;
  }

  if(kNPERef){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
      if(isAOD==1){
	  RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, dEdxhmAOD, kDefTOFs,
			     AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl,
			     kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      }
      else{
	  RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,
			     AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl,
			     kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      }
    /*RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, 
     *                    AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl, 
     *                    kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE);
     */
  }

  if(kNPERefTPConly){
    // **************************************************************
    // 
    // Reference task for TPC-only on the inclusive leg
    //
      // **************************************************************
      if(isAOD==1){
	  RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2AOD, dEdxhmAOD, 0,
			     AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl,
			     kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      }
      else{
	  RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, 0,
			     AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl,
			     kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      }
    /*RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, 
     *                    AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl, 
     *                    kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE);
     */
  }

  if(kNPEkAny){
    // **************************************************************
    // 
    // task for kAny instead of kBoth
    //
    // **************************************************************
    RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, 
                        AliHFEextraCuts::kAny, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl, 
                        kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
    /*RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, 
     *                    AliHFEextraCuts::kAny, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl, 
     *                    kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE);
     */
  }

  if(kNPEsystematics){
    // **************************************************************
    // 
    // Cut systematics on the associated track
    //
    // **************************************************************
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 0.0, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, 0.0, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, 3, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, 80, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, 100, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, 40, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, 40, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, 80, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, 100, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, 0.5, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, 2, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, 1, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, 4, dEdxaclm, dEdxachm, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm1, dEdxachm1, kTRUE, kTRUE);
      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm2, dEdxachm2, kTRUE, kTRUE);

  }

  // PID systematics
  if(kNPEpidsys){
    const int kPIDvars = 4;
    const double pidvars[kPIDvars] = {2, 2.5, 3.5, 4};
    double dEdxtestlow[12], dEdxtesthigh[12];
    for(int itest = 0; itest < kPIDvars; itest++){
      for(int icent = 0; icent < 12; icent++){
        dEdxtestlow[icent] = -1. * pidvars[itest];
        dEdxtesthigh[icent] = pidvars[itest];
      }
//      RegisterTaskNPEpPb( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpc1, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, kHFEV0A, kDefEtaIncMin, kDefEtaIncMax, kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxtestlow, dEdxtesthigh);
    } 
  }

  return NULL;
}

//===============================================================================
AliAnalysisTask *RegisterTaskNPEpPb(Bool_t useMC, Bool_t isAOD, 
               Int_t tpcCls=120, Int_t tpcClsPID=80, 
               Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
               Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL, 
               Double_t tofs=3., Int_t itshitpixel =AliHFEextraCuts::kBoth, Int_t icent=1,
               Double_t etaIncMin = -0.8, Double_t etaIncMax = 0.8,
               Double_t assETAm=-0.8, Double_t assETAp=0.8, Int_t assITS=2, Int_t assTPCcl=100,
               Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
               Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
               Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE)
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

  //
  // Cuts on the associated leg
  //
  Int_t iassDCAr = (Int_t)(assDCAr*10);
  Int_t iassDCAz = (Int_t)(assDCAz*10);
  Int_t iassTPCSminus = assTPCSplus ? (Int_t)(assTPCSplus[0]*1000.) : 0;

  printf("Argument passed to function to determine the centrality estimator %i\n", icent);
  if (icent == 2) TString cesti("V0M");
  else if (icent == 3) TString cesti("CL1");
  else if (icent == 4) TString cesti("ZNA");
  else TString cesti("V0A");
  printf("Centrality estimator %s\n", cesti.Data());

  TString appendix(TString::Format("SPD%d_incTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%d_photTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%dce%s",ipixelany,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSminus,cesti.Data()));
  printf("Add macro appendix %s\n", appendix.Data());

 if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpepPb")) gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/pPb/ConfigHFEnpepPb.C");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisTaskHFE *task = ConfigHFEnpepPb(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, icent, etaIncMin, etaIncMax,
					     assETAm, assETAp, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, assTPCSplus, useCat1Tracks, useCat2Tracks);
  if(isAOD)
    task->SetAODAnalysis();
  else
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
