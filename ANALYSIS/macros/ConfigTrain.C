{
   train_name      = "train_LHC09a4ESDanalysis";
   proof_cluster   = "alicecaf.cern.ch";
   useAFPAR        = 0;
   proof_dataset   = "/COMMON/COMMON/LHC09a4_run8100X#/esdTree";
   usePLUGIN       = 1;
   usePAR          = 0;
   useCPAR         = 0;
   root_version    = "v5-25-02";
   aliroot_version = "v4-18-08-AN";
   alien_datadir   = "/alice/sim/PDC_09/LHC09a4/";
   alien_outdir    = "/alice/sim/PDC_09/LHC09a4/analysis/ESD/pass1";
   maxMergeFiles   = 50;
   mergeExclude    = "AliAOD.root AliAOD.VertexingHF.root";
   nRunsPerMaster  = 10;
   nFilesPerJob    = 200;
   run_range[0]    = 81655;
   run_range[1]    = 81656;
   useDBG          = 0;
   useMC           = 1;
   useTAGS         = 0;
   useKFILTER      = 0;
   useTR           = 1;
   useCORRFW       = 0;
   useAODTAGS      = 0;
   saveTrain       = kFALSE;

   // Analysis modules
   iAODanalysis    = 0;
   iAODhandler     = 0;
   iESDfilter      = 0;
   iMUONcopyAOD    = 0; 
   iJETAN          = 0;  // AOD prod
   iPWG4partcorr   = 0;  // AOD prod
   iPWG4gammaconv  = 1;
   iPWG4omega3pi   = 1;
   iPWG3vertexing  = 0;  // AOD prod
   iPWG3hfe        = 1;
   iPWG2femto      = 1;
   iPWG2spectra    = 1;
     iPWG2protons      = 1;
     iPWG2checkcascade = 1;
     iPWG2perfcascade  = 1;
     iPWG2checkv0      = 1;
     iPWG2strange      = 1;
   iPWG2flow       = 1;
   iPWG2res        = 1;
   iPWG2kink       = 1;
     iPWG2kinkESDMC    = 1;
     iPWG2kinkres      = 1;
     iPWG2kinklikesign = 1;
   iPWG2unicor     = 1;
   iPWG2evchar     = 1;
   iPWG2forward    = 0;  // error

// Configuration fot the wagons
   configPWG2femto = $ALICE_ROOT/PWG2/FEMTOSCOPY/macros/Train/Train3/ConfigFemtoAnalysis.C;
}
