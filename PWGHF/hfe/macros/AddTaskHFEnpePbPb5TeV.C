AliAnalysisTask *AddTaskHFEnpePbPb5TeV(
                                   Bool_t MCthere = kFALSE,       // DATA
                                   //Bool_t MCthere = kTRUE,      // MC
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
  const int	kDefTPCcl	= 120;  // 100 (Andrea)
  const int	kDefTPCclPID	=  80;  // 90 (Andrea)
  const int kDefTPCclshared = 1.1;
  const int	kDefITScl	=   4;  // 5 (Andrea)
  const int kDefITSchi2percluster = -1; // cleanup removes badly matching tracks - effects high pt  (cut value = 36) ---> 36 default value for AOD
  const double	kDefDCAr	=   1.; // 2.4 (Andrea)
  const double	kDefDCAz	=   2.; // 3.2 (Andrea)
  const double	kDefTOFs	=   3.;
  const double	kDefITSs	=   2.; // -1,1 up to 1.5, -2,2 up to 3 (Andrea)
  const double  kDefEtaIncMin = -0.8;
  const double  kDefEtaIncMax = 0.8;
  const Bool_t   etacorrection   = kFALSE;
  const Bool_t   multicorrection = kFALSE;

  // --- TPC nsigma max and min cut ---
  Double_t dEdxhm[12] = {3.11,3.11,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};  
  Double_t tpcl1[12]  = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};                              // my 46% (mfaggin)
  Double_t tpcl2[12]  = {-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1};                  // my 50% (mfaggin)
  Double_t tpcl3[12]  = {0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18};                  // my 40% (mfaggin)
  Double_t tpcl4[12]  = {-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38};      // my 60% (mfaggin)

  // Default setting for the associated electron for the NonPhotonic Analysis
  const double	kassETAm = -0.8;                // -0.9 (Andrea)
  const double	kassETAp = 0.8;                 // 0.9 (Andrea)
  const int	kassITS		=   2;          // # cluster
  const int	kassTPCcl	= 60;           // 80 (Andrea)
  const int	kassTPCPIDcl	=  60;          // not used (Andrea) ---> directly in the filterbit
  const double	kassDCAr	=   1.0;        // not used (Andrea) ---> directly in the filterbit 2.4
  const double	kassDCAz	=   2.0;        // not used (Andrea) ---> directly in the filterbit 3.2
  const double	kassTPCSminus	=  -3.0;
  const double	kassTPCSplus	=   3.0;
  // --- Centrality selection ---
  const int     centrMin        = 0;
  const int     centrMax        = 10;

  Int_t kWei = -1;
  
        // if we want weights to be applied kWei MUST be 0!!!
  if (MCthere) kWei = 0; // default Pb-Pb
  enum {

    k11a10abisweiData = 6,  // LHC11a10abis weights 
    k11a10bplusweiData = 7,     // LHC11a10b_plus weights 
    k11a10bpluswei11a10abis   = 8,  // LHC11a10b_plus weights for LHC11a10abis
    k11a10bbisweiData = 9,  // LHC11a10bbis weights 

    k16g1 = 60         // 60: PbPb LHC16g1 minimum bias MC (mfaggin, 29/06/2017)
  };
  Int_t kWeightMC = k16g1;
  

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
                                                                
    // TPC low cut = 0 (tpcl1)
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // Additional tasks 

    // (sma, 6.6.2017) Vary the TOF cut to +-2 sigma)
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 2., 0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // (mfaggin, 15.6.2017) Vary the TOF cut to +-2.5 sigma)
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 2.5, 0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // (mfaggin, 15.6.2017) Vary the TOF cut to +-3.5 sigma)
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 3.5, 0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);

    // ----------------   
    if(MCthere) 
    {                                                                                                                  
        // TPC low cut = -0.1 (tpcl2) NO WEIGHTS
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
        // TPC low cut = -0.1 (tpcl2) WITH WEIGHTS
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    }
    else
    {
        // TPC low cut = -0.1 (tpcl2) NO WEIGHTS for data
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
    }
    // ----------------

    // TPC low cut = 0.16 (tpcl3)
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // TPC low cut = -0.38 (tpcl4)
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl4, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
/*
                // ---------------- TPC clusters (18.6.2017) ----------------
    // TPC cluster 110
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, 110, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // TPC cluster 115
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, 115, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // TPC cluster 125
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, 125, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // TPC cluster 130
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, 130, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);

                // ---------------- TPC clusters for PID (24.6.2017) ----------------
    // TPC cluster for PID 70
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // TPC cluster for PID 75
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, 75, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // TPC cluster for PID 85
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    // TPC cluster for PID 90
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 //kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
*/
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

   // TPC low cut = 0 (tpcl1)
/*
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm,  kDefTOFs,  kDefITSs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
*/
   // Additional tasks
/*
   // (sma, 6.6.2017) Vary the ITS cut to +-1
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs, 1., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
*/
   // (sma, 6.6.2017) Vary the ITS cut to +-3
/*
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs, 3., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
*/
   // TPC low cut = -0.1 (tpcl2)
   if(MCthere)
   {
        // WITH WEIGHTS   
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        // NO WEIGHTS
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
   }
   else
   {
        // NO WEIGHTS FOR DATA
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
   }
   // TPC low cut = 0.16 (tpcl3)
/*
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3, dEdxhm,  kDefTOFs,  kDefITSs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
*/
   // TPC low cut = -0.38 (tpcl4)
/*
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl4, dEdxhm,  kDefTOFs,  kDefITSs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
*/
   // TPC low cut = -0.1 (tpcl2) TOF+-2
/*
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, isBeauty, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  2.,  kDefITSs, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
*/

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
  // ------- to manage containers name with negative TPC low cut --------
  bool IsTPClowcutNegative = kFALSE;
  if(tpcdEdxcutlow) 
  {
        tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
        if(tpclow<0)
        {
                IsTPClowcutNegative = kTRUE;
                tpclow = 0 - tpclow;            // switched signed (ready to be used in the container name)
        }
  }
  // --------------------------------------------------------------------
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
    cweightsback += "_";
    cweightsback += wei;
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
  
  // ------- to manage containers name with negative TPC low cut --------
  TString appendix = "";                                                                   // letter 'm' added in this point (after TPCs)
  if(IsTPClowcutNegative)      appendix+=TString::Format("SPD%d_incTPCc%dTPCp%dITS%dDCAr%dz%dTPCsm%dTOFs%dITSs%dm%dt%d_photTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%d%s%s%s%s",ipixelany,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,iitss,imult,itofpos,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,cweightsback.Data(),cmvx.Data(),cbeauty.Data(),kfp.Data()); 
  else                         appendix+=TString::Format("SPD%d_incTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dITSs%dm%dt%d_photTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%d%s%s%s%s",ipixelany,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,iitss,imult,itofpos,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,cweightsback.Data(),cmvx.Data(),cbeauty.Data(),kfp.Data());
  //TString appendix(TString::Format("SPD%d_incTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dITSs%dm%dt%d_photTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%d%s%s%s%s",ipixelany,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,iitss,imult,itofpos,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,cweightsback.Data(),cmvx.Data(),cbeauty.Data(),kfp.Data()));      // old version
  // --------------------------------------------------------------------

  printf("Add macro appendix %s\n", appendix.Data());
  
  // --------------------------------------------------------------------
  // GRID version
  if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors_PbPb5TeV")) gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors_PbPb5TeV.C");
  if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpePbPb5TeV"))gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFEnpePbPb5TeV.C");
/*
  // GSI version
  // ----- my weights (mfaggin, 29/06/2017) -----
  //if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors")) gROOT->LoadMacro("$TRAIN_ROOT/util/hfe/configs/ConfigWeightFactors.C");
  if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors_PbPb5TeV")) gROOT->LoadMacro("$TRAIN_ROOT/util/hfe/configs/ConfigWeightFactors_PbPb5TeV.C");
  // --------------------------------------------
  if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpePbPb5TeV")) gROOT->LoadMacro("$TRAIN_ROOT/util/hfe/configs/ConfigHFEnpePbPb5TeV.C");
  // --------------------------------------------------------------------  
*/
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    AliAnalysisTaskHFE *task = ConfigHFEnpePbPb5TeV(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itss, itshitpixel, itschi2percluster, tpcsharedcluster, etacorr, multicorr, toflast, etaIncMin, etaIncMax,
					      assETAm, assETAp, assITS, assTPCcl, assTPCPIDcl,
					      assDCAr, assDCAz, assTPCSminus, assTPCSplus,
					      useCat1Tracks, useCat2Tracks, weightlevelback,usekfparticle);
  // old config file
  //AliAnalysisTaskHFE *task = ConfigHFEnpePbPb(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itss, itshitpixel, itschi2percluster, tpcsharedcluster, etacorr, multicorr, toflast, etaIncMin, etaIncMax,
  //					      assETAm, assETAp, assITS, assTPCcl, assTPCPIDcl,
  //					      assDCAr, assDCAz, assTPCSminus, assTPCSplus,
  //					      useCat1Tracks, useCat2Tracks, weightlevelback,usekfparticle);

  if(isAOD)
    task->SetAODAnalysis();
  else
    task->SetESDAnalysis();
  
  if (useMC)	task->SetHasMCData(kTRUE);
  else		task->SetHasMCData(kFALSE);

  //if(useMC&&(beauty || (weightlevelback>=0))) ConfigWeightFactors(task,kFALSE,wei);//2;For default PbPb
  if(useMC&&(beauty || (weightlevelback>=0))) ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei);//2;For default PbPb

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
