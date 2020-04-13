//////////////////////////////////////////////////////////////
//                                                          //
// filBIT x ---> at track level is tested the filterbit 2^x //
//                                                          //
//////////////////////////////////////////////////////////////
Int_t filBIT;

AliAnalysisTask *AddTaskHFEnpePbPb5TeV(
                                   //Bool_t MCthere = kFALSE,                     // DATA
                                   Bool_t MCthere = kTRUE,                    // MC
                                   Bool_t isAOD = kTRUE,
				   Bool_t kNPERef = kFALSE,                      // PID: TOF+TPC   ---> hardcoded kTRUE in AddTask_hfe_HFEnpePbPb5TeV.C
				   Bool_t kNPEkAny = kFALSE,
                                   Bool_t newCentralitySelection = kTRUE,       // kTRUE: new framework used; kFALSE: old framework used 
				   Bool_t kNPERefTPConly = kFALSE,
				   Bool_t kNPETOFITS = kTRUE,                   // PID: TOF+TPC+ITS
				   Bool_t kNPETOFlast = kFALSE,
				   Bool_t kNPEw      = kFALSE,
				   Bool_t kNPEkf  = kFALSE,
                                   Int_t  RunSyst = 0                            // switch statement useful for systematics (mfaggin, 06/07/2017)
                                   ,Int_t centrMin = 0          // min centrality
                                   ,Int_t centrMax = 10         // max centrality
                                   //,Int_t centrMin = 30          // min centrality
                                   //,Int_t centrMax = 50         // max centrality
                                   //,Int_t centrMin = 60          // min centrality
                                   //,Int_t centrMax = 80         // max centrality
                                   ,Bool_t kHadContSyst = kFALSE        // hadron contamination systematics
                                   ,Bool_t kPileUpIonutRejection = kFALSE
                                   ,Bool_t kCheckDCA = kFALSE    // additional check: set maximum value for DCA between inclusive and associated tracks (default: 3cm)
                                   ,Bool_t kPhiSystConsistentHadCont = kFALSE    // mfaggin, 13-Mar-2018
                                   ,Bool_t kRejectKinkParticles = kFALSE         // mfaggin, 12-Apr-2018
                                   ,Int_t chosen_filBIT = 2     // EXPONENT that then defines the filterbit (filterbit=2^chosenfilBIT) (mfaggin, 07-Jun-2018)
                                   ,Bool_t kTestDifferentFilBIT = kFALSE        // test a different filterbit (mfaggin, 07-Jun-2018)
                                   ,Bool_t kTestWeightSmallerCentBins = kFALSE   // test weights in smaller centrality bins (e.g.: 0-5% and 5-10% instead of 0-10%)
                                   ,Bool_t kTestNewWeights_MBfixedHIJING = kFALSE        // using weights for pi0, eta from MB production with fixed HIJING issues for pi0 decays
                                   ,Bool_t kTestNewWeights_MBfixedHIJING_asDefault = kTRUE        // using weights for pi0, eta from MB production with fixed HIJING issues for pi0 decays as default
                                   ,Bool_t kTestCookedWeightsFrom276 = kFALSE    // testing weights from charged pions at 5.02TeV cooked from charged pions at 2.76 TeV
                                   ,Bool_t kTestWeights_paperAfterCR1 = kTRUE   // testing weights from measured charged pions at 5 TeV (almost final, paper after CR1) (mfaggin, 13-Aug-2019)
                                   )		   
  
{
        filBIT = chosen_filBIT;
        // 
        // NB:
        //      - filBIT=2 : tested filterbit number 4  = kTrkITSConstraint 
        //      - filBIT=4 : tested filterbit number 16 = kTrkGlobalNoDCA
        //

        enum SystTipe
        {
                // ----- TOFonly cases -----
                 kTPCcutTOFcut  = 1     // syst. on TPC sigma and TOF sigma cut
                ,kTPCclst       = 2     // syst. on TPC cluster and TPC cluster for PID
                ,kITSclstANDspd = 3     // syst. on ITS cluster and SPD hits (*)
                ,kEta           = 4     // syst. on Eta interval
                ,kAssTPCclst    = 5     // syst. on associated particle TPC cluster and TPC cluster for PID 
                ,kAssITSclst    = 6     // syst. on associated particle ITS cluster 
                ,kAssMinpT      = 7     // syst. on associated particle minimum pT (mfaggin, 14th July 2017)
                ,kSPDAny        = 8     // (*) missing kAny case added
                // ----- ITS+TOF+TPC PID cases -----    (mfaggin, 25 January 2018)
                ,kTPCcut_ITS        = 9         // syst. on TPC PID cut
                ,kTOFcut_ITS        = 10        // syst. on TOF PID cut
                ,kTPCclst_ITS       = 11        // syst. on TPC cluster and TPC cluster for PID
                ,kITSclstANDspd_ITS = 12        // syst. on ITS cluster and SPD hits
                ,kEta_ITS           = 13        // syst. on Eta interval
                ,kAssTPCclst_ITS    = 14        // syst. on associated particle TPC cluster and TPC cluster for PID
                ,kAssITSclst_ITS    = 15        // syst. on associated particle ITS cluster
                ,kAssMinpT_ITS      = 16        // syst. on associated particle minimum pT
                ,kMix1              = 17        // syst. on many parameters simultaneous variation (PID + track cuts) 
                ,kMix2              = 18        // syst. on many parameters simultaneous variation (ass. track cuts) 
                // ----- syst. weights ----- (mfaggin, 06 March 2018)
                ,kWeightSyst        = 19        // syst. on weights used for the tagging efficiency
        };

        // mfaggin, 13-Mar-2018
        enum HadContFunc{
                 kPhi_0_14    = 1
                ,kPhi_14_326  = 2
                ,kPhi_326_2pi = 3
        };

  // Default settings (TOF-TPC PbPb)
  const int	kDefTPCcl	= 120;  // 100 (Andrea)
  const int	kDefTPCclPID	=  80;  // 90 (Andrea)
  const int kDefTPCclshared = 1.1;
  const int	kDefITScl	=   4;  // 5 (Andrea)
  const int kDefITSchi2percluster = -1; // cleanup removes badly matching tracks - effects high pt  (cut value = 36) ---> 36 default value for AOD
  const double	kDefDCAr	=   1.; // 2.4 (Andrea)
  const double	kDefDCAz	=   2.; // 3.2 (Andrea)
  const double	kDefTOFs	=   3.;
  //const double	kDefITSs	=   2.; // -1,1 up to 1.5, -2,2 up to 3 (Andrea)
  const double  kDefITSsMin = -4;       // mfaggin 15-Dec-2017
  const double  kDefITSsMax = 2;        // mfaggin 15-Dec-2017
  const double  kDefEtaIncMin = -0.8;
  const double  kDefEtaIncMax = 0.8;
  const Bool_t   etacorrection   = kFALSE;
  const Bool_t   multicorrection = kFALSE;

  // --- TPC nsigma max and min cut ---
  Double_t dEdxhm[12] = {3.11,3.11,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};  
  Double_t tpcl1[12]  = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};                              // my 46% (mfaggin)
  //Double_t tpcl2[12]  = {-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1};                  // my 50% (mfaggin) in 0-10%

  //Double_t tpcl3[12]  = {0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18};                  // my 40% (mfaggin) in 0-10%
  //Double_t tpcl4[12]  = {-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38,-0.38};      // my 60% (mfaggin) in 0-10%

  //Double_t tpcl3[12]  = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};                              // my 40% (mfaggin) in 0-10% (mfaggin, 06-Jul-2018)
  //Double_t tpcl4[12]  = {-0.42,-0.42,-0.42,-0.42,-0.42,-0.42,-0.42,-0.42,-0.42,-0.42,-0.42,-0.42};      // my 60% (mfaggin) in 0-10% (mfaggin, 06-Jul-2018)

  Double_t tpcl2[12];
  Double_t tpcl3[12];
  Double_t tpcl4[12];
  for(UInt_t i_tpcl=0; i_tpcl<12; i_tpcl++){

        //if(centrMin==0  && centrMax==10) tpcl2[i_tpcl]=-0.1;

        if(centrMin==0  && centrMax==10){
                tpcl2[i_tpcl]=-0.16;   // my 50% (mfaggin) in 0-10% (mfaggin, 06-Jul-2018)
                tpcl3[i_tpcl]= 0.1;
                tpcl4[i_tpcl]=-0.42;
        }
        if(centrMin==30 && centrMax==50){
                tpcl2[i_tpcl]= 0; 
                tpcl3[i_tpcl]= 0.23;
                tpcl4[i_tpcl]=-0.23;
        }
        if(centrMin==60 && centrMax==80){
                tpcl2[i_tpcl]= 0.2;     // (mfaggin, 20-Jun-2018)
                tpcl3[i_tpcl]= 0.43;
                tpcl4[i_tpcl]=-0.03;
        }
        if(centrMin==30 && centrMax==80){       // cases 30-50% and 60-80% merged in 1 file (mfaggin, 22-Jun-2018)
                printf("\nCase centrality 30-80%\n");
                tpcl2[i_tpcl]=0;                // NB: it MUST trigger a function present in the had_cont_functions file, even if it is not the right one for the specific centrality, otherwise even for the correct cases the contamination function is not passed to the task (the Bool_t returned by the ReadContaminationFunctions function is kFALSE if AT LEAST 1 function is not read)
                tpcl3[i_tpcl]=0;
                tpcl4[i_tpcl]=0;
                //if(i_tpcl==0 || i_tpcl==1){
                //        tpcl2[i_tpcl]=-0.16;
                //        tpcl3[i_tpcl]= 0.1;
                //        tpcl4[i_tpcl]=-0.42;
                //}
                if(i_tpcl==4 || i_tpcl==5){   
                        tpcl2[i_tpcl]= 0;
                        tpcl3[i_tpcl]= 0.23;
                        tpcl4[i_tpcl]=-0.23;
                }
                if(i_tpcl==7 || i_tpcl==8){    
                        tpcl2[i_tpcl]= 0.2;
                        tpcl3[i_tpcl]= 0.43;
                        tpcl4[i_tpcl]=-0.03;
                }
        }
        printf("i_tpcl=%d, tpcl2[%d]=%.2f\n",i_tpcl,i_tpcl,tpcl2[i_tpcl]);
  }

  // new case considered for ITS cut [-4,2] (mfaggin, 09 January 2018)
  Double_t dEdxhmALTERNATIVE[12] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0};                   // with tpcl2: about 47.5% (mfaggin)

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
  // --- min pT as. particle (mfaggin, 14th July 2017) ---
  // const double kassMinpT = 0.1;      // not used... hardcoded as default value in the RegisterTaskNPEPbPb function

  // --- Centrality selection ---
  //const int     centrMin        = 0;  // put as argument of the AddTask itself
  //const int     centrMax        = 10;

  Int_t kWei = -1;
  
        // if we want weights to be applied kWei MUST be 0!!!
  if (MCthere) kWei = 0; // default Pb-Pb
  enum {

    k11a10abisweiData = 6,  // LHC11a10abis weights 
    k11a10bplusweiData = 7,     // LHC11a10b_plus weights 
    k11a10bpluswei11a10abis   = 8,  // LHC11a10b_plus weights for LHC11a10abis
    k11a10bbisweiData = 9,  // LHC11a10bbis weights 

    k16g1 = 60,         // 60: PbPb LHC16g1 minimum bias MC (mfaggin, 29/06/2017)
    k16g1_2 = 61         // 61: PbPb LHC16g1 minimum bias MC using charged pion data spectra (mfaggin, 26/07/2017)

    // --- systematics on charged pion weights --- (mfaggin, 06-Mar-2018)
    ,k16g1_2_default = 62       // 62: equal to 61, but necessary to evaluate upper limit-lower limit systematics (small differences introduced due to daata spectra fits)
    ,k16g1_2_systUp  = 63       // 63: upper limit charged pion weights
    ,k16g1_2_systLow = 64       // 64: lower limit charged pion weights
    ,k16g1_2_systTiltUpDown = 65       // 65: tilting up-down charged pion weights
    ,k16g1_2_systTiltDownUp = 66       // 66: tilting down-up charged pion weights

    ,k16g1_3 = 67               // 67: PbPb LHC16g1 minimum bias MC using pi0 data spectra (mfaggin, 06-Mar-2018)

    ,k16g1_2_3050 = 70          // 70: PbPb LHC16g1 minimum bias MC using pi charged  data spectra for 30-50% centrality class (mfaggin, 07-Jun-2018)
    ,k16g1_2_6080 = 71          // 71: PbPb LHC16g1 minimum bias MC using pi charged data spectra for 60-80% centrality class (mfaggin, 20-Jun-2018)
    ,k16g1_2_all = 72           // 72: PbPb LHC16g1 minimum bias MC using pi charged data spectra for 0-10%, 30-50% and 60-80% centrality classes, all in one file (mfaggin, 22-Jun-2018)
    ,k16g1_2_smallercentbins = 73       // 73: PbPb LHC16g1 minimum bias MC using pi charged data spectra in smaller centrality bins (e.g.: 0-5% and 5-10% instead of 0-10%)

    ,k18e1_fixedHIJING = 80     // 80: PbPb LHC18e1 minimum bias MC with fixed HIJING issue on pi0 decay - pi charged data spectra used for weights (11/09/2018)

    // weights from charged pions at 5.02 TeV cooked from charged pions at 2.76 TeV with the ratio of charged particle spectrum at the two energy as follows:
    //  CHPION(5.02TeV) = CHPIONS(2.76TeV)*CHPART(5.02TeV)/CHPART(2.76TeV)
    ,k18e1_fixHIJING_weightsFrompiCharged276 = 81

    // weights from measured charged pions ate 5 TeV almost final (paper after CR1)
    ,k18e1_fixHIJING_weightsPion5TeV_afterCR1 = 82

  };
  //Int_t kWeightMC = k16g1;
  Int_t kWeightMC = 0;
  if(centrMin==0  && centrMax==10)      kWeightMC=k16g1_2;
  if(centrMin==30 && centrMax==50)      kWeightMC=k16g1_2_3050;
  if(centrMin==60 && centrMax==80)      kWeightMC=k16g1_2_6080;
  if(centrMin==30 && centrMax==80)      kWeightMC=k16g1_2_all;  // (mfaggin, 22-Jun-2018)

  // using weights for pi0, eta from MB production with fixed HIJING issues for pi0 decays as default
  if(kTestNewWeights_MBfixedHIJING_asDefault)   kWeightMC=k18e1_fixedHIJING;
  
  // check outputs
  printf("\n\n===== Centrality: %d,%d\n",centrMin,centrMax);
  printf("===== tpcl2[] (lower edge TPC cut) = {");
        for(UInt_t ilowtpc=0; ilowtpc<12; ilowtpc++){
                if(ilowtpc==11) printf("%.2f",tpcl2[ilowtpc]);
                else            printf("%.2f,",tpcl2[ilowtpc]);
        }
  printf("}\n");
  printf("===== kWeightMC (MC weight parameter): %d\n\n",kWeightMC);


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


  //const Bool_t isBeauty = kFALSE; // should be false to prevent inclusive analysis    ---> fixed to kFALSE into the RegisterTask function due to the large amount of arguments

  if(kNPERef){
    // **************************************************************
    // 
    // Reference task
    //
    // ************************************************************** 
        printf("\n#####\n");
        printf("##### kNPERef case");
        printf("\n#####\n");
    // ----------------   
    if(MCthere) 
    {                                                                                                                  
        // TPC low cut = -0.1 (tpcl2) NO WEIGHTS
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
        // TPC low cut = -0.1 (tpcl2) WITH WEIGHTS
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
    }
    else
    {
        // TPC low cut = -0.1 (tpcl2) NO WEIGHTS for data
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
    }
    // ----------------      
    // ****************************************************************
    // systematics
    // ****************************************************************
    // ----------------                     
    switch(RunSyst)
    {
        case kTPCcutTOFcut:     // *** syst. on TPC sigma and TOF sigma cut ***
                cout << endl << "###############################################################" << endl;
                cout         << "###   running systematics on TPC sigma and TOF sigma cut    ###" << endl;
                cout         << "###############################################################" << endl;
                // TPC low cut = 0 (tpcl1)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC low cut = 0.16 (tpcl3)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC low cut = -0.38 (tpcl4)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl4, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // (sma, 6.6.2017) Vary the TOF cut to +-2 sigma)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 2., 0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // (mfaggin, 15.6.2017) Vary the TOF cut to +-2.5 sigma)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 2.5, 0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // (mfaggin, 15.6.2017) Vary the TOF cut to +-3.5 sigma)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 3.5, 0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kTPCclst:          // *** syst. on TPC cluster and TPC cluster for PID ***
                cout << endl << "#######################################################################" << endl;
                cout         << "###   running systematics on TPC cluster and TPC cluster for PID    ###" << endl;
                cout         << "#######################################################################" << endl;
                // TPC cluster 110
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
110, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster 115
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
115, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster 125
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
125, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster 130
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
130, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster for PID 70
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster for PID 75
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, 75, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster for PID 85
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster for PID 90
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kITSclstANDspd:    // *** syst. on ITS cluster and SPD hits ***
                cout << endl << "############################################################" << endl;
                cout         << "###   running systematics on ITS cluster and SPD hits    ###" << endl;
                cout         << "############################################################" << endl;
                // ITS cluster 2
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, 2, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // ITS cluster 3
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, 3, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // ITS cluster 5
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // ITS cluster 6
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, 6, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // SPD hit kFirst
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kFirst, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // SPD hit kSecond
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kSecond, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // SPD hit kNone
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kNone, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // SPD hit kAny
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kAny, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                
        break;

        case kEta:
                cout << endl << "#########################################################" << endl;
                cout         << "###   running systematics on inclusive eta interval   ###" << endl;
                cout         << "#########################################################" << endl;
                // |Eta|<0.9
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, -0.9, 0.9,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // |Eta|<0.7
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, -0.7, 0.7,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // |Eta|<0.6
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, -0.6, 0.6,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // |Eta|<0.5
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, -0.5, 0.5,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kAssTPCclst:
                cout << endl << "###########################################################################################" << endl;
                cout         << "###   running systematics on associated particle TPC cluster and TPC cluster for PID    ###" << endl;
                cout         << "###########################################################################################" << endl;
                // associated particle TPC cluster 40
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 40, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster 50
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 50, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster 70
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 70, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster 80
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 80, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster for PID 40
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, 40, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster for PID 50
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, 50, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster for PID 70
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, 70, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster for PID 80
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, 80, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC); 
        break;

        case kAssITSclst:
                cout << endl << "##################################################################" << endl;
                cout         << "###   running systematics on associated particle ITS cluster   ###" << endl;
                cout         << "##################################################################" << endl;
                // associated particle ITS cluster 1
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 1, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle ITS cluster 3
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 3, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle ITS cluster 4
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 4, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle ITS cluster 5
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 5, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle ITS cluster 6
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 6, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kAssMinpT:        
                cout << endl << "#################################################################" << endl;
                cout         << "###   running systematics on associated particle minimum pT   ###" << endl;
                cout         << "#################################################################" << endl;
                // associated particle minimum pT 0.0
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.0);
                // associated particle minimum pT 0.05
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.05);
                // associated particle minimum pT 0.15
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.15);
                // associated particle minimum pT 0.20
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.20);
                // associated particle minimum pT 0.25
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.25);
                // associated particle minimum pT 0.30
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.30);
        break;

        case kSPDAny:
                cout << endl << "##############################################" << endl;
                cout         << "###   running AliHFEextraCuts::kAny case   ###" << endl;
                cout         << "##############################################" << endl;
                // SPD hit kAny
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kAny, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;
    }
  }

  if(kNPETOFlast){
    // **************************************************************
    // 
    // Apply TOF after TPC for mismatch background studies
    //
    // **************************************************************
        printf("\n#####\n");
        printf("##### kNPETOFlast case");
        printf("\n#####\n");
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kTRUE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
  }

/*
  if(kNPEkAny){
    // **************************************************************
    // 
    // task for kAny instead of kBoth
    //
    // **************************************************************
        printf("\n#####\n");
        printf("##### kNPEkAny case");
        printf("\n#####\n");
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kAny, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
  }
*/

  if(kNPEw && MCthere){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
        printf("\n#####\n");
        printf("##### kNPEw && MCthere case");
        printf("\n#####\n");
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE, 0,kWei);
  }

  if(kNPERefTPConly){
    // **************************************************************
    // 
    // Reference task
    //
    // **************************************************************
        printf("\n#####\n");
        printf("##### kNPERefTPConly case");
        printf("\n#####\n");
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, 0.,0.,0, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
  }


 if(kNPETOFITS){
   // **************************************************************
   // 
   // Reference task
   //
   // **************************************************************                                                                            
        printf("\n#####\n");
        printf("##### kNPETOFITS case");
        printf("\n#####\n");
   // TPC low cut = 0 (tpcl1)
/*
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
*/
   // Additional tasks
/*
   // (sma, 6.6.2017) Vary the ITS cut to +-1
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs, 1., 1, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
*/
   // (sma, 6.6.2017) Vary the ITS cut to +-3
/*
   RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs, 3., 3., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
*/
   // TPC low cut = -0.1 (tpcl2)
   if(MCthere)
   {
        // WITH WEIGHTS   
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        // NO WEIGHTS
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);

        // weights from measured charged pions ate 5 TeV almost final (paper after CR1)
        if(kTestWeights_paperAfterCR1){
                printf("\n#####\n");
                printf("##### kTestWeights_paperAfterCR1 case - weights from measured charged pions at 5 TeV, almost final (paper after CR1)");
                printf("\n#####\n");  
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k18e1_fixHIJING_weightsPion5TeV_afterCR1);
        }

        // testing weights from charged pions at 5.02TeV cooked from charged pions at 2.76 TeV
        if(kTestCookedWeightsFrom276){
                printf("\n#####\n");
                printf("##### kTestCookedWeightsFrom276 case - weights from charged pions cooked from 2.76TeV one tested");
                printf("\n#####\n");     
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k18e1_fixHIJING_weightsFrompiCharged276);
        }

        if(kTestNewWeights_MBfixedHIJING && !kTestNewWeights_MBfixedHIJING_asDefault){
                printf("\n#####\n");
                printf("##### kTestNewWeights_MBfixedHIJING case - fixed HIJING weights tested");
                printf("\n#####\n");      
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k18e1_fixedHIJING);  
        }

        if(kTestDifferentFilBIT){
                const Int_t original_filBIT = filBIT;
                if(original_filBIT==2)   filBIT=4;
                if(original_filBIT==4)   filBIT=2;
                // WITH WEIGHTS changed filBIT
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // NO WEIGHTS changed filBIT
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
        }
        filBIT = chosen_filBIT; // original filBIT restored

        if(kNPEkAny){
        // **************************************************************
        // 
        // task for kAny instead of kBoth
        //
        // **************************************************************
                printf("\n#####\n");
                printf("##### kNPEkAny case (in addition to kBoth one) for MC");
                printf("\n#####\n");
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kAny, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        }

        if(kTestWeightSmallerCentBins){
        // *****************************************************************************
        //
        // weights in smaller centrality bins (e.g.: 0-5% and 5-10% instead of 0-10%)
        //
        // *****************************************************************************
                printf("\n=====\n");
                printf("\n===== kTestWeightSmallerCentBins case for MC");
                printf("\n=====\n");
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k16g1_2_smallercentbins);

        }

   }
   else
   {
        // NO WEIGHTS FOR DATA
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
        if(kTestDifferentFilBIT){
                const Int_t original_filBIT = filBIT;
                if(original_filBIT==2)   filBIT=4;
                if(original_filBIT==4)   filBIT=2;
                // NO WEIGHTS changed filBIT
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
                }
        filBIT = chosen_filBIT; // original filBIT restored
                if(kPhiSystConsistentHadCont){
                        // hadron contamination function used: phi [0,1.4] (mfaggin, 13-Mar-2018)
                        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                        kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	                kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kFALSE,kFALSE,kPhi_0_14);
                        // hadron contamination function used: phi [1.4,3.26] (mfaggin, 13-Mar-2018)
                        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                        kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	                kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kFALSE,kFALSE,kPhi_14_326);
                        // hadron contamination function used: phi [3.26,2pi] (mfaggin, 13-Mar-2018)
                        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                        kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	                kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kFALSE,kFALSE,kPhi_326_2pi);
                }
        if(kNPEkAny){
        // **************************************************************
        // 
        // task for kAny instead of kBoth
        //
        // **************************************************************
                printf("\n#####\n");
                printf("##### kNPEkAny case (in addition to kBoth one) for DATA");
                printf("\n#####\n");
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kAny, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                }
   }

   switch(RunSyst){
        
        case kTPCcut_ITS:
                cout << endl << "#################################################" << endl;
                cout         << "###   running systematics on TPC sigma cut    ###" << endl;
                cout         << "#################################################" << endl;
                // WITH WEIGHTS - ITS cut [-4,2], TPC cut [0.18,3] (about 40%)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl3, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
                kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // WITH WEIGHTS - ITS cut [-4,2], TPC cut [-0.1,2] (about 47.5%)
                //RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                //kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhmALTERNATIVE,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // WITH WEIGHTS - ITS cut [-4,2], TPC cut [-0.38,3] (about 60%)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl4, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kTOFcut_ITS:
                cout << endl << "#################################################" << endl;
                cout         << "###   running systematics on TOF sigma cut    ###" << endl;
                cout         << "#################################################" << endl;   
                // Vary the TOF cut to +-2 sigma
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 2., kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // Vary the TOF cut to +-2.5 sigma
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 2.5, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // Vary the TOF cut to +-3.5 sigma
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, 3.5, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // Vary the ITS cut to +-3 sigma (symmetric cut)
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, -3, 3, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);           
        break;

        case kTPCclst_ITS:
                cout << endl << "#######################################################################" << endl;
                cout         << "###   running systematics on TPC cluster and TPC cluster for PID    ###" << endl;
                cout         << "#######################################################################" << endl;
                // TPC cluster 110
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                110, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster 115
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                115, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster 125
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                125, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster 130
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                130, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster for PID 70
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster for PID 75
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, 75, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster for PID 85
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TPC cluster for PID 90
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kITSclstANDspd_ITS:
                cout << endl << "############################################################" << endl;
                cout         << "###   running systematics on ITS cluster and SPD hits    ###" << endl;
                cout         << "############################################################" << endl;
                // ITS cluster 2
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, 2, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // ITS cluster 3
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, 3, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // ITS cluster 5
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, 5, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // ITS cluster 6
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, 6, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
		kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // SPD hit kFirst
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kFirst, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // SPD hit kSecond
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kSecond, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // SPD hit kNone
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kNone, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // SPD hit kAny
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs, kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kAny, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kEta_ITS:
                cout << endl << "#########################################################" << endl;
                cout         << "###   running systematics on inclusive eta interval   ###" << endl;
                cout         << "#########################################################" << endl;
                // |Eta|<0.7
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, -0.7, 0.7,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // |Eta|<0.6
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, -0.6, 0.6,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // |Eta|<0.5
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, -0.5, 0.5,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // 0<Eta<0.8
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, 0, 0.8,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // -0.8<Eta<0
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, -0.8, 0,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kAssTPCclst_ITS:
                cout << endl << "###########################################################################################" << endl;
                cout         << "###   running systematics on associated particle TPC cluster and TPC cluster for PID    ###" << endl;
                cout         << "###########################################################################################" << endl;
                // associated particle TPC cluster 40
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017 
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 40, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster 50
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 50, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster 70
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 70, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster 80
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 80, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster for PID 40
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, 40, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster for PID 50
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, 50, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster for PID 70
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, 70, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle TPC cluster for PID 80
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, 80, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC); 
        break;

        case kAssITSclst_ITS:
                // associated particle ITS cluster 1
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 1, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle ITS cluster 3
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 3, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle ITS cluster 4
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 4, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle ITS cluster 5
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 5, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // associated particle ITS cluster 6
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 6, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kAssMinpT_ITS:
                // associated particle minimum pT 0.0
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.0);
                // associated particle minimum pT 0.05
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.05);
                // associated particle minimum pT 0.15
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.15);
                // associated particle minimum pT 0.20
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.20);
                // associated particle minimum pT 0.25
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.25);
                // associated particle minimum pT 0.30
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.30);
        break;

        case kMix1:
                // TOF cut [-3,3], TPC cluster 110, TPC cluster for PID 70, ITS cluster 1  
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                110, 70, 1, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  3.,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-3,3], TPC cluster 125, TPC cluster for PID 85, ITS cluster 2
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                125, 85, 2, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  3.,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-3,3], TPC cluster 115, TPC cluster for PID 75, ITS cluster 3 
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                115, 75, 3, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  3.,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-3,3], TPC cluster 130, TPC cluster for PID 80, ITS cluster 5
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                130, 90, 5, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  3.,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-2.5,2.5], TPC cluster 110, TPC cluster for PID 70    
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                110, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  2.5,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-2.5,2.5], TPC cluster 115, TPC cluster for PID 75    
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                115, 75, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  2.5,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-2.5,2.5], TPC cluster 125, TPC cluster for PID 85    
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                125, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  2.5,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-2.5,2.5], TPC cluster 130, TPC cluster for PID 80    
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                130, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  2.5,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-3.5,3.5], TPC cluster 110, TPC cluster for PID 70    
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                110, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  3.5,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-3.5,3.5], TPC cluster 115, TPC cluster for PID 75    
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                115, 75, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  3.5,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-3.5,3.5], TPC cluster 125, TPC cluster for PID 85    
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                125, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  3.5,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
                // TOF cut [-3.5,3.5], TPC cluster 130, TPC cluster for PID 80    
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                130, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  3.5,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
        break;

        case kMix2:
                // ass. TPC cluster 40, ass. TPC cluster for PID 40, ass. ITS cluster 1, associated particle minimum pT 0.10
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 1, 40, 40, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.10);
                // ass. TPC cluster 50, ass. TPC cluster for PID 50 ass. ITS cluster 3, associated particle minimum pT 0.10
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 3, 50, 50, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.10);
                // ass. TPC cluster 70, ass. TPC cluster for PID 70 ass. ITS cluster 4, associated particle minimum pT 0.10
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 4, 70, 70, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.10);
                // ass. TPC cluster 80, ass. TPC cluster for PID 80 ass. ITS cluster 5, associated particle minimum pT 0.10
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, 5, 80, 80, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.10);
                // ass. TPC cluster 40, ass. TPC cluster for PID 40, associated particle minimum pT 0.15
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 40, 40, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.15);
                // ass. TPC cluster 50, ass. TPC cluster for PID 50, associated particle minimum pT 0.15
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 50, 50, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.15);
                // ass. TPC cluster 70, ass. TPC cluster for PID 70, associated particle minimum pT 0.15
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 70, 70, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.15);
                // ass. TPC cluster 80, ass. TPC cluster for PID 80, associated particle minimum pT 0.15
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 80, 80, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.15);
                // ass. TPC cluster 40, ass. TPC cluster for PID 40, associated particle minimum pT 0.20
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 40, 40, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.20);
                // ass. TPC cluster 50, ass. TPC cluster for PID 50, associated particle minimum pT 0.20
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 50, 50, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.20);
                // ass. TPC cluster 70, ass. TPC cluster for PID 70, associated particle minimum pT 0.20
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 70, 70, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.20);
                // ass. TPC cluster 80, ass. TPC cluster for PID 80, associated particle minimum pT 0.20
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, 80, 80, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.20);
        break;

        case kWeightSyst:
                cout << endl << "##################################################################" << endl;
                cout         << "###   running systematics on weights for tagging efficiency    ###" << endl;
                cout         << "##################################################################" << endl;
        // WITH WEIGHTS default ch. pions  
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k16g1_2_default);
        // WITH WEIGHTS upper limit ch. pions  
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k16g1_2_systUp);
        // WITH WEIGHTS lower limit ch. pions  
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k16g1_2_systLow);
        // WITH WEIGHTS tilting up-down ch. pions  
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k16g1_2_systTiltUpDown);
        // WITH WEIGHTS tilting down-up ch. pions  
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k16g1_2_systTiltDownUp);
        // WITH WEIGHTS pi0  
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,k16g1_3);
        break;

   }
 }

 
  if(kNPEkf){
    // **************************************************************
    // 
    // Use KF particle
    //
    // **************************************************************
        printf("\n#####\n");
        printf("##### kNPEkf case");
        printf("\n#####\n");
    RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl1, dEdxhm, kDefTOFs,0.,0., AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1,2,kFALSE,kFALSE,kFALSE,kTRUE);
  }



  if(kPileUpIonutRejection){
        // ****************************************************************************************
        // 
        // Pile-up rejection based on number of kTPCout tracks with VZERO multiplicity correlation
        //
        // ****************************************************************************************
        printf("\n#####\n");
        printf("##### kPileUpIonutRejection case");
        printf("\n#####\n");

        // pile-up rejection applied
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kTRUE);

  }
/*
 if(kCheckDCA){                                                                            
        printf("\n#####\n");
        printf("##### kCheckDCA case");
        printf("\n#####\n");

                // *** change DCA max between inclusive and associated tracks ***
                // max DCA 2cm
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kFALSE,2.0);
                // max DCA 1.5cm
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kFALSE,1.5);
                // max DCA 1cm
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kFALSE,1.0);
                // max DCA 0.5cm
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm, kDefTOFs,kDefITSsMin,kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
	        kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kFALSE,0.5);

 }*/

if(kHadContSyst){
        printf("\n#####\n");
        printf("##### kHadContSyst case");
        printf("\n#####\n");
                // WITH WEIGHTS   
                RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
                kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kFALSE, kDefEtaIncMin, kDefEtaIncMax,
                kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC,0.1,kFALSE//,3.0
                ,kHadContSyst);
}

if(kRejectKinkParticles){
        printf("\n#####\n");
        printf("##### kRejectKinkParticles case");
        printf("\n#####\n");
        // WITH WEIGHTS   
        RegisterTaskNPEPbPb( centrMin,centrMax,newCentralitySelection,MCthere, isAOD, //isBeauty, // mfaggin 15-Dec-2017
kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl2, dEdxhm,  kDefTOFs,  kDefITSsMin, kDefITSsMax, AliHFEextraCuts::kBoth, kDefITSchi2percluster, kDefTPCclshared, etacorrection, multicorrection, kTRUE, kDefEtaIncMin, kDefEtaIncMax,
			//kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,-1);
			 kassETAm, kassETAp, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kTRUE, kFALSE,kWei,kWeightMC);
}


  
  return NULL;

}

//===============================================================================

//===============================================================================
AliAnalysisTask *RegisterTaskNPEPbPb(
                                     Int_t centrMin = 0, Int_t centrMax = 100,
                                     Bool_t newCentralitySelection = kTRUE, // kTRUE: new framework used; kFALSE: old framework used 
                                     Bool_t useMC, Bool_t isAOD, //Bool_t beauty,
				     Int_t tpcCls=120, Int_t tpcClsPID=80, 
				     Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0, 
				     Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL, 
				     Double_t tofs=3., //Double_t itss=0., 
                                     Double_t itssMin=0., Double_t itssMax=0.,  // mfaggin 15-Dec-2017
                                     Int_t itshitpixel =AliHFEextraCuts::kBoth,
				     Double_t itschi2percluster = -1, Double_t tpcsharedcluster = 1.1,
				     Bool_t etacorr=kFALSE, Bool_t multicorr = kFALSE, Bool_t rejkinks = kFALSE,
				     Double_t etaIncMin = -0.8, Double_t etaIncMax = 0.8,
				     Double_t assETAm=-0.8, Double_t assETAp=0.8, Int_t assITS=2, Int_t assTPCcl=100,
				     Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
				     Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
				     Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE,
				     Int_t weightlevelback = -1,Int_t wei = 2,
                                     Double_t assMinpTvalue = 0.1,      // associated particle minimum pT (mfaggin, 14th July 2017)
                                     Bool_t applyPileUpRejection = kFALSE    // --- Pile-up rejection using correlation kTPCout-VZEROmult ---
                                     //,Double_t dcaMAX = 3.0     // DCA max between inclusive and associated tracks (mfaggin, 20-Feb-2018)
                                     ,Bool_t hadcontsyst = kFALSE
                                     ,Int_t phisystconsistenthadcont = 0        // mfaggin, 13-Mar-2018
				     //Bool_t releasemcvx = kFALSE,
				     //Bool_t ipCharge = kFALSE,
				     //Bool_t ipOpp = kFALSE
                                     //,
				     //Bool_t usekfparticle = kFALSE
                                     )
{
        // *** hardcoded value due to cint-limitation on number of parameters (mfaggin, 14th July 2017)
        Bool_t usekfparticle = kFALSE;
        Bool_t beauty = kFALSE;         // mfaggin 15-Dec-2017
	Bool_t releasemcvx = kFALSE;    // mfaggin 02-Feb-2018
	Bool_t ipCharge = kFALSE;       // mfaggin 02-Feb-2018
        Bool_t ipOpp = kFALSE;          // mfaggin 02-Feb-2018
        Bool_t toflast = kFALSE;        // mfaggin 12-Apr-2018

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
        tpclow = (Int_t) (tpcdEdxcutlow[7]*1000.);
        if(tpclow<0)
        {
                IsTPClowcutNegative = kTRUE;
                tpclow = 0 - tpclow;            // switched signed (ready to be used in the container name)
        }
  }
  Int_t tpchigh = (Int_t) (tpcdEdxcuthigh[4]*1000.);
  // --------------------------------------------------------------------
  Int_t itofs = (Int_t)(tofs*10.);
  Int_t iitssMin = (Int_t)(itssMin*10.);
        Int_t iitsMinNEGATIVE = 0;
        // ----- name with negative ITSPID low cut -----
        bool IsITSlowCutNegative = kFALSE;
        if(iitssMin<0){
                IsITSlowCutNegative = kTRUE;
                iitsMinNEGATIVE = 0-iitssMin;          // switched signed (ready to be used in the container name) 
        }
        // ---------------------------------------------
  Int_t iitssMax = (Int_t)(itssMax*10.);
  Int_t ipixelany = itshitpixel;
  Int_t imult = multicorr ? 1 : 0;
  Int_t itofpos = toflast ? 1 : 0;

  //
  // Cuts on the associated leg
  //
  Int_t iassDCAr = (Int_t)(assDCAr*10);
  Int_t iassDCAz = (Int_t)(assDCAz*10);
  Int_t iassMinpT = (Int_t) (assMinpTvalue*100);
  Int_t iassTPCSplus  = assTPCSplus ? (Int_t)(assTPCSplus[0]*1000) : 0;
  Int_t icat1 = useCat1Tracks ? 1 : 0;
  Int_t icat2 = useCat2Tracks ? 1 : 0;
  Int_t etaInclusiveMax = (Int_t)(etaIncMax*100);
  
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
  
// new version, mfaggin 15-Dec-2017
  //TString appendix(Form("SPD%d_incEta%dTPCc%dTPCp%dITS%dDCAr%dz%d",ipixelany,etaInclusiveMax,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz));

// new version, mfaggin 02-Aug-2018
  TString appendix(Form("SPD%d_incEta",ipixelany));
  int intEtaMin = static_cast<int> (etaIncMin*100);
  int intEtaMax = static_cast<int> (etaIncMax*100);
  if(intEtaMin != 0-intEtaMax){
        if(intEtaMin<0){
                intEtaMin = 0-intEtaMin;
                appendix+= TString::Format("m%d%dTPCc%dTPCp%dITS%dDCAr%dz%d",intEtaMin,intEtaMax,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz);
        }
        else    appendix+= TString::Format("%d%dTPCc%dTPCp%dITS%dDCAr%dz%d",intEtaMin,intEtaMax,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz);
  }
  else appendix+= TString::Format("%dTPCc%dTPCp%dITS%dDCAr%dz%d",etaInclusiveMax,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz);




  if(centrMin==30 && centrMax==80){
        printf("\n----------- tpclow %d\n",tpclow);
        if(tpclow>195 && tpclow<205)        appendix+="TPCPIDeff50";       // with multiple centralities the reference to the mean TPC value has no sense ---> TPC PID value used (NB: the min value of TPC checked is the one for 60-80% just for simplicity, namely tpcdEdxcutlow[7], line 1228)
        if(tpclow>425 && tpclow<435)           appendix+="TPCPIDeff40";
        if(tpclow>25 && tpclow<35)             appendix+="TPCPIDeff60";
  }
  else{
        if(IsTPClowcutNegative)       appendix+=TString::Format("TPCsMinm%d",tpclow);
        else                          appendix+=TString::Format("TPCsMin%d",tpclow);
        appendix+=TString::Format("Max%d",tpchigh);
  }

  appendix+=TString::Format("TOFs%d",itofs);
  if(IsITSlowCutNegative)       appendix+=TString::Format("ITSsMinm%d",iitsMinNEGATIVE);
  else                          appendix+=TString::Format("ITSsMin%d",iitssMin);
  appendix+=TString::Format("ITSsMax%dm%dt%d_photMinpT%dTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%d%s%s%s%s",iitssMax,imult,itofpos,iassMinpT,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,cweightsback.Data(),cmvx.Data(),cbeauty.Data(),kfp.Data());

  // DCA max between inclusive and associated tracks (mfaggin, 20-Feb-2018)
  //Int_t intmaxDCA = (Int_t) (dcaMAX*10);
  //if(dcaMAX<3.0)                appendix+=TString::Format("DCAmax%d",intmaxDCA);

  appendix += TString::Format("_filBIT%d",filBIT);
  if(applyPileUpRejection)      appendix+="_PileUpRejIonutCut";
  if(hadcontsyst)               appendix+="_hadContSyst";
  if(rejkinks)                  appendix+="_kinksReject";

  if(phisystconsistenthadcont==1)       appendix+="_hadcontphi014";      // mfaggin, 13-Mar-2018
  if(phisystconsistenthadcont==2)       appendix+="_hadcontphi14326";    // mfaggin, 13-Mar-2018
  if(phisystconsistenthadcont==3)       appendix+="_hadcontphi3262pi";   // mfaggin, 13-Mar-2018

/*
  // ------- to manage containers name with negative TPC low cut --------
  TString appendix = "";                                                                   // letter 'm' added in this point (after TPCs)
  if(IsTPClowcutNegative)      appendix+=TString::Format("SPD%d_incEta%dTPCc%dTPCp%dITS%dDCAr%dz%dTPCsm%dTOFs%dITSs%dm%dt%d_photMinpT%dTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%d%s%s%s%s",ipixelany,etaInclusiveMax,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,iitss,imult,itofpos,iassMinpT,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,cweightsback.Data(),cmvx.Data(),cbeauty.Data(),kfp.Data()); 
  else                         appendix+=TString::Format("SPD%d_incEta%dTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dITSs%dm%dt%d_photMinpT%dTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%d%s%s%s%s",ipixelany,etaInclusiveMax,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,iitss,imult,itofpos,iassMinpT,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,cweightsback.Data(),cmvx.Data(),cbeauty.Data(),kfp.Data());
*/
  //TString appendix(TString::Format("SPD%d_incTPCc%dTPCp%dITS%dDCAr%dz%dTPCs%dTOFs%dITSs%dm%dt%d_photTPCc%dTPCp%dITS%dDCAr%dDCAz%dTPCs%d%s%s%s%s",ipixelany,tpcCls,tpcClsPID,itsCls,idcaxy,idcaz,tpclow,itofs,iitss,imult,itofpos,assTPCcl,assTPCPIDcl,assITS,iassDCAr,iassDCAz,iassTPCSplus,cweightsback.Data(),cmvx.Data(),cbeauty.Data(),kfp.Data()));      // old version
  // --------------------------------------------------------------------

  printf("Add macro appendix %s\n", appendix.Data());
  
  // --------------------------------------------------------------------
  // GRID version
  if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors_PbPb5TeV")) gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors_PbPb5TeV.C");
  if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpePbPb5TeV"))gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFEnpePbPb5TeV.C");
  // GSI version
  // ----- my weights (mfaggin, 29/06/2017) -----
  //if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors_PbPb5TeV")) gROOT->LoadMacro("$TRAIN_ROOT/util/hfe/configs/ConfigWeightFactors_PbPb5TeV.C");
  // --------------------------------------------
  //if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpePbPb5TeV")) gROOT->LoadMacro("$TRAIN_ROOT/util/hfe/configs/ConfigHFEnpePbPb5TeV.C");
  // --------------------------------------------------------------------  

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    AliAnalysisTaskHFE *task = ConfigHFEnpePbPb5TeV(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, //itss, 
                                              itssMin, itssMax,         // mfaggin 15-Dec-2017
                                              itshitpixel, itschi2percluster, tpcsharedcluster, etacorr, multicorr, toflast, etaIncMin, etaIncMax,
					      assETAm, assETAp, assITS, assTPCcl, assTPCPIDcl,
					      assDCAr, assDCAz, assTPCSminus, assTPCSplus,
					      useCat1Tracks, useCat2Tracks, weightlevelback
                                              ,assMinpTvalue    // associated particle minimum pT syst. (mfaggin, 14th July 2017)
                                              //,dcaMAX           // DCA max between inclusive and associated tracks (mfaggin, 20-Feb-2018)
                                              ,hadcontsyst      // had. cont. systematics
                                              ,phisystconsistenthadcont // phi variation with consistent hadron contamination function (mfaggin, 13-Mar-2018)
                                              ,rejkinks         // kink particles rejection
                                              ,filBIT           // EXPONENT that then defines the filterbit (filterbit=2^chosenfilBIT)
                                              );
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
  if(useMC&&(beauty || (weightlevelback>=0))){
        if(wei==67)                 ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei,"nonHFEcorrect_PbPb5TeV_frompi0.root");
        else if(wei==70 || wei==71) ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei,"nonHFEcorrect_PbPb5TeV_fromchpions_3050.root");
        else if(wei==72)            ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei,"nonHFEcorrect_PbPb5TeV_fromchpions_allCentr.root");
        else if(wei==73)            ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei,"nonHFEcorrect_PbPb5TeV_fromchpions_smallerbins.root");
        else if(wei==80)            ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei,"nonHFEcorrect_PbPb5TeV_fromchpions_allCentr_newMB_fixedHIJING.root");
        else if(wei==81)            ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei,"nonHFEcorrect_PbPb5TeV_fromchpions_ScaledFrom276TeV_allCentr_newMB_fixedHIJING.root");
        else if(wei==82)            ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei,"nonHFEcorrect_PbPb5TeV_fromchpions_afterCR1_allCentr_newMB_fixedHIJING.root");
        else                        ConfigWeightFactors_PbPb5TeV(task,kFALSE,wei);
  }

  // ----- centrality selection -----
  task->SetCentralityCheck(newCentralitySelection,"V0M");
  task->SetCentralityInterval(centrMin,centrMax);               // all events outside the desired centrality interval are rejected
  // --------------------------------
  // ----- trigger selecton ---------
  if(!newCentralitySelection)   task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral); // old framework
  if(newCentralitySelection)    task->SelectCollisionCandidates(AliVEvent::kINT7);                                               // new framework
  // --------------------------------
  // --- Pile-up rejection using correlation kTPCout-VZEROmult ---
  task->SetApplyPileUpMultiplicitySelection(applyPileUpRejection);
  // -------------------------------------------------------------
  
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
