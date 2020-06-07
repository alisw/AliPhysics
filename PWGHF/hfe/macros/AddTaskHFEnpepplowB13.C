#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include "$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHFEnpepplowB13.C"
#include "$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors.C"

AliAnalysisTask *RegisterTaskNpepplowB13(Bool_t useMC, Int_t RunSystematic = 0, Bool_t isAOD,
                                   Int_t tpcCls=120, Int_t tpcClsPID=80,
                                   Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0,
                                   Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL,
                                   Double_t tofs=3., Int_t itshitpixel =AliHFEextraCuts::kBoth,
                                   Int_t iKink = 0, Int_t assITS=2, Int_t assTPCcl=100,
                                   Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
                                   Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
                                   Double_t assITSpid = 3.0, Double_t assTOFpid = 0.0, Double_t assMinPt = 0.0,
                                   Double_t assETAm = -0.9,Double_t assETAp = 0.9,
                                   Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kFALSE,
                                   Int_t weightlevelback = -1, Int_t WhichWei = 0,Int_t HadronContFunc=0);
#endif


AliAnalysisTask *AddTaskHFEnpepplowB13(Bool_t MCthere,
                                 Bool_t isAOD = kTRUE,
                                 Bool_t kNPERef = kTRUE,
                                 Bool_t kNPERefTPConly = kFALSE,
                                 Int_t RunSystematic = 0 // select systematic type
){
    
    
   enum SystematicType
    {
        kSystMultiple = 1, // ITS, hadron contamination, SPD kAny, PrimaryVertexTyp
        kSystTPCcluster= 2,
        kSystPID = 3,
        //        kSystAsociatedMultiple = 6, // ITS hits & standalone, TPC PID, TOF PID
        kSystAsociatedDCA = 4,
        kSystAsociatedTPCcluster = 5,
        kSystAsociatedMinpTWeights = 6,
	kSystTrackPIDMixedCuts = 7,
	kSystAssoTrackPIDAssoMinPtMixedCuts = 8,
	kHC = 9,
	kSystAsociatedMinpTNew = 10,
	kSystAssoTrackPIDAssoMinPtMixedCutsNew = 11,
	kSystTrackPIDMixedCutsNew = 12,
	kWeights = 13
    };

    
    // Default settings (TOF-TPC pp)
    // ESD analysis of LHC15n, 5 TeV analysis
    
    const int	kDefTPCcl	= 100;
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
    
    // On ESD:
    Double_t dEdxhm[12] = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};  // Above 3 sigma we neglect 0.13%
    Double_t dEdxhm20[12] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0};
    Double_t dEdxhm25[12] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5};
    Double_t dEdxhm35[12] = {3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5};

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
    Double_t tpcl15[12] = {-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};
    Double_t tpcl16[12]  = {-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};  // 90%
    
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
    const double assETAm        =   -0.8;
    const double assETAp        =   0.8;
    const int	kasspTmin		=    0.0;
    
    
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
    
    // For 5 TeV analysis, Feb 20, 2017: no re-weighting at all, for the moment
    if (MCthere) kWei = 0;
    
    enum {
        // Keep some examples from the 2.76 TeV analysis. To be replaced.
        kWeiLHC11b10a = 19,    // weights for published pp @ 2.76 TeV: LHC11b10a, Pythia min. bias
        kWeiLHC11b10b = 20,    // pass2                                LHC11b10b, Pythia + HF + pi0
        kWeiLHC12a9 = 21,      //                                      LHC12a9: Pythia + HF->e
        kWeiLHC12e6 = 22,      //                                      LHC12e6: Pythia min. bias
        kWeiLHC18h1_tu = 59,      //                                      LHC18h1: tilted up weights
        kWeiLHC18h1_td = 60,      //                                      LHC18h1: tilted down weights
    };
    int kWeiData;
    // The re-weighting concerns the photonic sources. The tagging efficiency is best taken
    // from min bias MC for statistics reasons. Therefore the default is put for min bias MC.
    kWeiData = 58; //d12 for the low magnetic field case
    
    
    if(kNPERef){
        // **************************************************************
        //
        // Reference task
        //
        // **************************************************************
        RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE);
        
        RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData, 2);
      
        if (MCthere){
            RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                      dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                      kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
            
            switch (RunSystematic){
            
            case kWeights:
                         
            RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                      dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                      kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiLHC18h1_tu);
                                      
            RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                      dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                      kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiLHC18h1_td);
                                      
           break;
                                      
           }
        
                }
                
        
        if (!MCthere){
            switch (RunSystematic) {
                case kSystMultiple:
                    // Hadron Contamination
                    // choice of the hadron contamination function - only Reference choice (TPC lower cut: -1 sigma cut, kBoth); default is Error function
                    RegisterTaskNpepplowB13(  MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                      dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                      kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData,1);
                                      
                    RegisterTaskNpepplowB13(  MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                      dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                      kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData,3);
                    
                    break;
                case kHC:
                    RegisterTaskNpepplowB13(  MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                      dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                      kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData,1);
                             
                    RegisterTaskNpepplowB13(  MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                      dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                      kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData,3);
                    break;
                case kSystPID:
                    break;
            }
        }
        switch (RunSystematic) {
            case kSystMultiple:
                // ITS hits
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, 2, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, 4, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                // SPD selection
                // kFirst
                RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kFirst, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                // SPD kAny
                RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kAny, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                // DCA Var
                RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 2.4, 3.2, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 0.5, 1, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                break;
                
            case kSystTPCcluster:
                
                // TPC clusters
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 2.0, 3.0, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                                  
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 1.5, 3.5, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData); 
                                  
                 RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 1.8, 2.8, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData); 
                                  
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, 2.8, 4.2, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);                  
                                  
                //RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 90, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13, dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl, kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                /*RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 95, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 105, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 110, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 125, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 130, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                // TPC PID custer
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, 95, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, 100, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);*/
                
                break;
                
            case kSystPID:
                // variation of the lower TPC cut
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl11,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl16,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                // vary the top TPC PID cut with -1 sigma cut
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm20, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm25, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm35, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                // different TOF PID cuts
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, 2, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, 2.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, 3.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, 4, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                break;
                
                
                //            default:
                //
                //                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                //                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                //                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                //
                //                break;
       
   

	case kSystTrackPIDMixedCuts:

               RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl11,
                                  dEdxhm, 2, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 105, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0,
                                  dEdxhm, 2, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 110, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl16,
                                  dEdxhm, 2, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
               
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 130, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm20, 2.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, 90, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm25, 2.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, 85, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm35, 2.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                

                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, 95, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm20, 3.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, 100, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm35, 3.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 90, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm25, 3.5, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 95, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm35, 4, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 125, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl11,
                                  dEdxhm, 4, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 120, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl16,
                                  dEdxhm, 4, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                break;

       case kSystAsociatedDCA:
                // DCA associated cuts
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  2.4,3.2, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13(MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl,kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  0.5,1, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

            break;
        case kSystAsociatedTPCcluster:
                
                // Assosciated TPC clusters
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, 50,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, 70,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, 80,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
               
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, 90,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
            
                
                break;
                

case kSystAsociatedMinpTWeights:

		// AssociatedMinPt		

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.1,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
        
		break;



case kSystAssoTrackPIDAssoMinPtMixedCuts:

		// Associated Mixed Cuts 		

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.1,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70, 80,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 50,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.1,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70, 50,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 90,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                          
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, 50,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, 80,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                          
		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 80,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70, 90,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 90,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                          
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 90,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		break;


case kSystAsociatedMinpTNew:

		// AssociatedMinPt		

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.1,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.12,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.16,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.18,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.2,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.24,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.28,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.3,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		break;



case kSystAssoTrackPIDAssoMinPtMixedCutsNew:

		// Associated Mixed Cuts 		

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.12,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70, 80,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 50,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.12,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70, 50,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 90,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.12,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                          
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, 50,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.12,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, 80,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                          
		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.12,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 80,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 90, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.12,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 70, 90,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 90,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.12,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                          
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 80, 90,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl13,
                          dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, 50, 70,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,0.14,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);

		break;


              
            case kSystTrackPIDMixedCutsNew:
                
                // TPC clusters
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 90, 90, kDefITScl, 2.8, 4.2, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 95, 85, kDefITScl,  2.8, 4.2, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 100, 95, kDefITScl,  1.8, 2.8, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 105, 90, kDefITScl,  1.8, 2.8, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 110, 90, kDefITScl,  2.8, 4.2, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 105, 85, kDefITScl,  2.0, 3.0, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 110, 95, kDefITScl,  2.0, 3.0, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                // TPC PID custer
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 105, 90, kDefITScl,  1.5, 3.5, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 95, 95, kDefITScl,  1.5, 3.5, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 105, 95, kDefITScl,  2.0, 3.0, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, 95, 90, kDefITScl,  1.8, 2.8, tpcl13,
                                  dEdxhm, kDefTOFs, AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                                  kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
                
                break;
  

        }
        

     
        
    }
    if(kNPERefTPConly){
        // **************************************************************
        //
        // Reference task for TPC-only on the inclusive leg
        //
        // **************************************************************
        //TPC-only task: TPC PID cuts: 0 sigma
        if (MCthere){
            RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0,
                              dEdxhm, 0., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                              kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE, kWei, kWeiData);
        }
        
        RegisterTaskNpepplowB13( MCthere, RunSystematic, isAOD, kDefTPCcl, kDefTPCclPID, kDefITScl, kDefDCAr, kDefDCAz, tpcl0,
                          dEdxhm, 0., AliHFEextraCuts::kBoth, 0, kassITS, kassTPCcl, kassTPCPIDcl,
                          kassDCAr, kassDCAz, dEdxaclm, dEdxachm, kassITSpid, kassTOFpid,kasspTmin,assETAm,assETAp, kTRUE, kFALSE);
    }
    
    
    
    return NULL;
}

//===============================================================================
AliAnalysisTask *RegisterTaskNpepplowB13(Bool_t useMC, Int_t RunSystematic = 0, Bool_t isAOD,
                                   Int_t tpcCls=120, Int_t tpcClsPID=80,
                                   Int_t itsCls=4, Double_t dcaxy=1.0, Double_t dcaz=2.0,
                                   Double_t *tpcdEdxcutlow=NULL, Double_t *tpcdEdxcuthigh=NULL,
                                   Double_t tofs=3., Int_t itshitpixel =AliHFEextraCuts::kBoth,
                                   Int_t iKink = 0, Int_t assITS=2, Int_t assTPCcl=100,
                                   Int_t assTPCPIDcl=80, Double_t assDCAr=1.0, Double_t assDCAz=2.0,
                                   Double_t *assTPCSminus = NULL, Double_t *assTPCSplus=NULL,
                                   Double_t assITSpid = 3.0, Double_t assTOFpid = 0.0, Double_t assMinPt = 0.0,
                                   Double_t assETAm = -0.9,Double_t assETAp = 0.9,
                                   Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kFALSE,
                                   Int_t weightlevelback = -1, Int_t WhichWei = 0,Int_t HadronContFunc=0)
{
    // Fixed values
    Double_t etaIncMin = -0.8; Double_t etaIncMax = 0.8;
    //    Double_t assETAm=-0.9; Double_t assETAp=0.9;
    
    //
    // Cuts on the inclusive leg
    //
    Int_t iassMinPt = (Int_t)(assMinPt*100.);
    Int_t idcaxy = (Int_t)(dcaxy*10.);
    Int_t idcaz = (Int_t)(dcaz*10.);
    Int_t tpclow = 0;
    Int_t tpchigh = 0;
    if(tpcdEdxcutlow) tpclow = (Int_t) (tpcdEdxcutlow[0]*1000.);
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
    if (useCat1Tracks) phoTrack = 1; //eta_ass09 and ptmin 00
    if (useCat2Tracks) phoTrack = 2; //eta_ass08 and ptmin 01
    
    TString cweightsback("");
    if(weightlevelback>=0) {
        cweightsback += "Wa";
        if (WhichWei>0){
            cweightsback += WhichWei;
            //cweightsback += weightlevelback;
        }
    }
    
   

 TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCls%dTPCh%dTOFs%dK%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dAssMinpt%dtr%d%s%dSyst%d",
                                     tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,tpchigh,itofs,iKink,assTPCcl,assTPCPIDcl,assITS,iassDCAr,
                                     iassDCAz,iassTPCSminus,iassITS,iassTOF,iassMinPt,phoTrack,cweightsback.Data(),HadronContFunc, RunSystematic));
    
 
 /*   
    TString appendix(TString::Format("incTPCc%dp%dITS%dSPD%dDCAr%dz%dTPCls%dTPCh%dTOFs%dK%d_photTPCc%dp%dITS%dDCAr%dz%dTPCs%dITSs%dTOFs%dtr%d%s%d",
                                     tpcCls,tpcClsPID,itsCls,ipixelany,idcaxy,idcaz,tpclow,tpchigh,itofs,iKink,assTPCcl,assTPCPIDcl,assITS,iassDCAr,
                                     iassDCAz,iassTPCSminus,iassITS,iassTOF,phoTrack,cweightsback.Data(),HadronContFunc));
    */
    printf("Add macro appendix %s\n", appendix.Data());
    
    if(useMC&&!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigWeightFactors"))
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/ConfigWeightFactors.C");
    if(!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigHFEnpepplowB13"))
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHFEnpepplowB13.C");
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    //mgr->AddClassDebug("AliHFENonPhotonicElectron", 1);
    AliAnalysisTaskHFE *task = ConfigHFEnpepplowB13(useMC, isAOD, appendix, tpcCls, tpcClsPID, itsCls, dcaxy, dcaz,
                                              tpcdEdxcutlow, tpcdEdxcuthigh, tofs, 0, itshitpixel, iKink,etaIncMin, etaIncMax,
                                              assETAm, assETAp, assMinPt, assITS, assTPCcl, assTPCPIDcl, assDCAr, assDCAz, assTPCSminus,
                                              assTPCSplus,assITSpid,assTOFpid, useCat1Tracks, useCat2Tracks, weightlevelback,HadronContFunc);
    
    if(isAOD)
        task->SetAODAnalysis();
    else
        task->SetESDAnalysis();
    
    if (useMC)	task->SetHasMCData(kTRUE);
    else		task->SetHasMCData(kFALSE);
    
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    
    if(useMC && weightlevelback>=0) {
        //ConfigWeightFactors(task,kFALSE,WhichWei,"nonHFEcorrect_pp13_LowB_Final_Sep28.root");
        ConfigWeightFactors(task,kFALSE,WhichWei,"nonHFEcorrect_pp13_LowB_Final_April30_2020.root");
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
