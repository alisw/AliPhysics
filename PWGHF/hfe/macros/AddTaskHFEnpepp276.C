AliAnalysisTask *AddTaskHFEnpepp276(Bool_t MCthere,
                                    Bool_t isAOD,
                                    Bool_t kNPERef = kTRUE,
                                    Bool_t kNPERefTPConly = kFALSE,
                                    Int_t RunSystematic = 0 // select systematic type
){

   enum SystematicType
   {
      kSystMultiple = 1, // Kink, ITS, chi2/TPCcluster, hadron contamination, SPD kAny
      kSystDCA = 2,
      kSystTPCcluster= 3,
      kSystTPCPID = 4,
      kSystTOFPID = 5,
      kSystAsociatedMultiple = 6, // ITS hits & standalone, TPC PID, TOF PID
      kSystAsociatedDCA = 7,
      kSystAsociatedTPCcluster = 8,
      kSystAsociatedMinpTWeights = 9
   };
   // Default settings (TOF-TPC pp)
   // ESD analysis of LHC11a, 2.76 TeV analysis

   const int	kDefTPCcl	= 110;
   const int	kDefTPCclPID	=  80;
   const int	kDefITScl	=   3;
   const double	kDefDCAr	=   1.;
   const double	kDefDCAz	=   2.;
   const double	kDefTOFs	=   3.;

   // TPC PID Cuts Inclusive leg:
   // General, if mean=0 and sigma=1:
   // Above 3 sigma we neglect 0.13%.
   // Cut in sigma (effective efficiency from cut to 3 sigma)
   // -1 (84%), -0.75 (77.2%), -0.5 (69%), -0.25 (59.7%), -0.129 (55%)October 4
   //  0 (49.9%), 0.122 (45%), 0.25 (40%), 0.5 (30.7%)

   Double_t dEdxhm[12] = {3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2};  // Above 3 sigma we neglect 0.13%
   Double_t dEdxhm20[12] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0};
   Double_t dEdxhm21[12] = {1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};
   Double_t dEdxhm22[12] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
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
   Double_t tpcl11[12]  = {-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5};  // 69%
   Double_t tpcl12[12]  = {-0.75,-0.75,-0.75,-0.75,-0.75,-0.75,-0.75,-0.75,-0.75,-0.75,-0.75,-0.75};  // 77.3%
   Double_t tpcl13[12]  = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};  // 84%
   Double_t tpcl14[12]  = {-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3};  // 90%
   Double_t tpcl15[12] = {-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2,-3.2};

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

   // For 2.76 TeV analysis, March 15, 2016:
   if (MCthere) kWei = 0;

   enum {
      kWeiLHC11b10a = 19,    // weights for published pp @ 2.76 TeV: LHC11b10a, Pythia min. bias
      kWeiLHC11b10b = 20,    // pass2                                LHC11b10b, Pythia + HF + pi0
      kWeiLHC12a9 = 21,      //                                      LHC12a9: Pythia + HF->e
      kWeiLHC12e6 = 22,      //                                      LHC12e6: Pythia min. bias
      kWeiLHC12f1a = 23,     // pass4                                LHC12f1a: Pythia min. bias
      kWeiLHC12f1b = 24,     //                                      LHC12f1b: Phojet min. bias
      kWeiLHC12f1a_f1b = 25, // MC closure: 12f1a to 12f1b
      kWeiLHC12f1b_f1a = 26, //             12f1b to 12f1a
      kWeiLHC12e6_tu = 27,   // pion spectrum tilted 'up'            LHC12e6: Pythia min. bias
      kWeiLHC12e6_td = 28,    // pion spectrum tilted 'down'          LHC12e6: Pythia min. bias
      kWeiLHC11b10b_12e6 = 29 // MC closure: 11b10b to 12e6
   };
   int kWeiData;
   // The re-weighting concerns the photonic sources. The tagging efficiency is best taken
   // from min bias MC for statistics reasons. Therefore the default is put for min bias MC.
   kWeiData = kWeiLHC12e6;
   //kWeiData = kWeiLHC12f1a;


   if(kNPERef){
      // **************************************************************
      //
      // Reference task
      //
      // **************************************************************

      if (!MCthere){
         // TPC lower cut at 0 sigma

         switch (RunSystematic) {
            case kSystMultiple:
               // Hadron Contamination
               // choice of the hadron contamination function - only Reference choice (TPC lower cut: -1 sigma cut, kBoth); default is Error function
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData,1);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData,2);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData,3);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData,4);

               break;
            case kSystDCA:
               break;
            case kSystTPCcluster:
               break;

            case kSystTPCPID:
               // variation of the lower TPC cut
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl5,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl10,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl11,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl12,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl14,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

               // vary the top TPC PID cut with -1 sigma cut
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm22, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm21, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm20, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
               break;

            case kSystTOFPID:
               break;
            case kSystAsociatedMultiple:
               break;
            case kSystAsociatedDCA:
               break;
            case kSystAsociatedTPCcluster:
               break;
            case kSystAsociatedMinpTWeights:
               break;

            default:
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
               break;
         }
      }

      // Reference (for MC with reference weight LHC12e6 for pass2)
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                        dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                        kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

      switch (RunSystematic) {
         case kSystMultiple:
            // include kink mothers
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 1, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            // ITS hits
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, 2, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, 4, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            // chi2/TPCcluster
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData, 0, 0.1, 3);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData, 0, 0.1, 5);

            // SPD selection
            // kFirst
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            // SPD kAny
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kAny, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            break;
         case kSystDCA:

            // DCA Var
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 2, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 0.5, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, 4, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, 1, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            break;

         case kSystTPCcluster:

            // TPC clusters
            RegisterTaskNPEpp( MCthere, isAOD, 100, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, 105, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, 115, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, 120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            // TPC PID custer
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, 75, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            break;

         case kSystTPCPID:
            break;

         case kSystTOFPID:

            // different TOF PID cuts
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, 2., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, 2.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, 3.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, 4., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, 5., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            break;

         case kSystAsociatedMultiple:
            // Number of ITS hits

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, 1, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, 3, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            // Try using again category 2 tracks (ITS SA)
            // and vary the number of hits from 2 to 4 (max!) because the ITS PID is required

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, 3, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, 4, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE, kWei, kWeiData);

            //TOF PID and TPC PID
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm1, dEdxachm1, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm2, dEdxachm2, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, 3.0, kTRUE, kFALSE, kWei, kWeiData);

            break;

         case kSystAsociatedDCA:
            // DCA associated
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              2, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              0.5, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              0.2, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, 4, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, 1, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, 0.5, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            break;

         case kSystAsociatedTPCcluster:
            //TPC cluster
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 50,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 60,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 70,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 80,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70,50,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 60,50,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 50,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            break;

         case kSystAsociatedMinpTWeights:
            //min pT associated
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData, 0, 0.0);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData, 0, 0.05);
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData, 0, 0.15);

            break;

         default:
            // SPD selection
            // kFirst
            RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                              dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

            break;
      }

      // different re-weighting (only for MC)
      if (MCthere){
         // with no weights for Tracking
         RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                           dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                           kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

         switch (RunSystematic) {
            case kSystMultiple:

               // include kink mothers
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 1, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);


               // ITS hits
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, 2, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, 4, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               // chi2/TPCcluster
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, -1, 0, 0, 0.1, 3);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, -1, 0, 0, 0.1, 5);
               // SPD selection
               // SPD kFirst
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);
               // SPD kAny
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kAny, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               break;

            case kSystDCA:

               // DCA Var
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 2, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 0.5, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, 4, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, 1, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               break;

            case kSystTPCcluster:

               // TPC clusters
               RegisterTaskNPEpp( MCthere, isAOD, 100, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, 105, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, 115, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, 120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               // TPC PID cluster
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, 75, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, 70, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               break;

            case kSystTPCPID:
               break;

            case kSystTOFPID:

               // different TOF PID cuts
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, 2., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, 2.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, 3.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, 4., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, 5., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               break;

               // Some systematics on the associated leg
            case kSystAsociatedMultiple:
               // Number of ITS hits
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, 1, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, 3, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               // Try using again category 2 tracks (ITS SA)
               // and vary the number of hits from 2 to 4 (max!) because the ITS PID is required

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, 3, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, 4, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kTRUE);

               //TOF PID and TPC PID
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm1, dEdxachm1, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm2, dEdxachm2, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, 3.0, kTRUE, kFALSE);

               break;

            case kSystAsociatedDCA:
               // DCA associated
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 2, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 0.5, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 0.2, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, 4, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, 1, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, 0.5, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               break;

            case kSystAsociatedTPCcluster:
               // TPC PID cluster
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 50,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 60,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 70,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 80,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               // TPC PID cluster
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70,50,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 60,50,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 50,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);

               break;

            case kSystAsociatedMinpTWeights:
               //tilted weights
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12e6_tu);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12e6_td);

               // min pT associated
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, -1, 0, 0, 0.0);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, -1, 0,  0, 0.05);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, -1, 0, 0, 0.15);

               // min pT associated + tilted weight
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12e6_tu, 0, 0.0);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12e6_tu, 0, 0.05);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12e6_tu, 0, 0.15);

               // min pT associated + tilted weight
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12e6_td, 0, 0.0);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12e6_td, 0, 0.05);

               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12e6_td, 0, 0.15);
               break;
            default:
               // SPD selection
               // SPD kFirst
               RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                 dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                 kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);
               break;
         }

         /*
          // MC closure on pass4 MC samples
          RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12f1a_f1b);
          RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiLHC12f1b_f1a);
          */

      }
   }
   //}

   if(kNPERefTPConly){
      // **************************************************************
      //
      // Reference task for TPC-only on the inclusive leg
      //
      // **************************************************************

      //TPC-only task: TPC PID cuts: 0 sigma
      RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0,
                        dEdxhm, 0., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                        kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE, kWei, kWeiData);

      if (MCthere){
         RegisterTaskNPEpp( MCthere, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0,
                           dEdxhm, 0., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                           kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid, kTRUE, kFALSE);
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
                                   Int_t weightlevelback = -1, Int_t WhichWei = 0,
                                   Int_t HadronContFunc=0, Double_t assMinPt = 0.1, Int_t Chi2perTPCcluster=4 )
{
   // Fixed values
   Double_t etaIncMin = -0.8; Double_t etaIncMax = 0.8;
   Double_t assETAm=-0.8; Double_t assETAp=0.8;
   // Test May 19,2016
   //Double_t etaIncMin = -0.9; Double_t etaIncMax = 0.9;
   //Double_t assETAm   = -0.9; Double_t assETAp=0.9;

   //
   // Cuts on the inclusive leg
   //
   Int_t idcaxy = (Int_t)(dcaxy*10.);
   Int_t idcaz = (Int_t)(dcaz*10.);
   Int_t tpclow = 0;
   if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
   Int_t tpchigh = 0;
   if(tpcdEdxcuthigh) tpchigh = (Int_t) (tpcdEdxcuthigh[0]*1000.);
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

   if (Chi2perTPCcluster!=4) {
      TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCs%dTOFs%dK%dChi%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%d%s",
                                       tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,itofs,iKink,Chi2perTPCcluster,assTPCcl,assTPCPIDcl,assITS,
                                       iassDCAr,iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,cweightsback.Data()));
   } else if (assMinPt!=0.1) {
      Int_t iMinPt = assMinPt*100;
      TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCs%dTOFs%dK%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%dminpT%d%s",
                                       tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,itofs,iKink,assTPCcl,assTPCPIDcl,assITS,
                                       iassDCAr,iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,iMinPt,cweightsback.Data()));
   } else if (tpchigh[1]!=3200) {
      TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCs%dTPCtop%dTOFs%dK%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%d%s",
                                       tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,tpchigh,itofs,iKink,assTPCcl,assTPCPIDcl,assITS,iassDCAr,
                                       iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,cweightsback.Data()));
   }else if(HadronContFunc!=0){
      TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCs%dTOFs%dK%dHadFunc%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%d%s",
                                       tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,itofs,iKink,HadronContFunc,assTPCcl,assTPCPIDcl,assITS,iassDCAr,
                                       iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,cweightsback.Data()));
   } else {
      TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCs%dTOFs%dK%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%d%s",
                                       tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,itofs,iKink,assTPCcl,assTPCPIDcl,assITS,iassDCAr,
                                       iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,cweightsback.Data()));
   }

   printf("Add macro appendix %s\n", appendix.Data());

   if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors"))
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors.C");
   if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpepp276"))
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHFEnpepp276.C");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
   //mgr->AddClassDebug("AliHFENonPhotonicElectron", 1);
   AliAnalysisTaskHFE *task = ConfigHFEnpepp276(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz, 
                                                tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, iKink,etaIncMin, etaIncMax,
                                                assETAm, assETAp, assMinPt, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus, 
                                                assTPCSplus,assITSpid,assTOFpid, useCat1Tracks, useCat2Tracks, weightlevelback, HadronContFunc,
                                                Chi2perTPCcluster);
   
   if(isAOD)
      task->SetAODAnalysis();
   else
      task->SetESDAnalysis();
   
   if (useMC)	task->SetHasMCData(kTRUE);
   else		task->SetHasMCData(kFALSE);
   
   task->SelectCollisionCandidates(AliVEvent::kMB);
   
   if(useMC && weightlevelback>=0) {
      ConfigWeightFactors(task,kFALSE,WhichWei,"nonHFEcorrect_pp276.root");
   }
   
   //create data containers
   TString containerName = mgr->GetCommonFileName();
   containerName += ":HFEtask";
   containerName += appendix.Data();
   printf("container name: %s\n", containerName.Data());
   
   
   task->ConnectOutput(1, mgr->CreateContainer(Form("HFE_Results_%s", appendix.Data()), TList::Class(),
                                               AliAnalysisManager::kOutputContainer, containerName.Data()));
   task->ConnectOutput(2, mgr->CreateContainer(Form("HFE_QA_%s", appendix.Data()), TList::Class(),
                                               AliAnalysisManager::kOutputContainer, containerName.Data()));
   
   mgr->ConnectInput(task,  0, cinput );
   
   mgr->AddTask(task);
   
   return NULL;
}
