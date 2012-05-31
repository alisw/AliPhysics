{
   train_name      = "train_LHC09a3ESD";
   proof_cluster   = "alicecaf.cern.ch";
   useAFPAR        = 0;
   proof_dataset   = "/COMMON/COMMON/LHC09a4_run8100X#/esdTree";
   usePLUGIN       = 1;
   usePAR          = 0;
   useCPAR         = 0;
   root_version    = "v20091109";
   aliroot_version = "v4-18-10-AN";
   alien_datadir   = "/alice/sim/PDC_08b/LHC09a3";
   alien_outdir    = "/alice/sim/PDC_08b/LHC09a3/AOD";
   maxMergeFiles   = 50;
   mergeExclude    = "AliAOD.root AliAOD.VertexingHF.root AliAOD.Jets.root deltaAODPartCorr.root";
   nRunsPerMaster  = 10;
   nFilesPerJob    = 10;
   run_range[0]    = 90000;
   run_range[1]    = 90025;
   useDBG          = 0;
   useMC           = 1;
   useTAGS         = 1;
   useKFILTER      = 1;
   useTR           = 1;
   useCORRFW       = 1;
   useAODTAGS      = 0;
   saveTrain       = kFALSE;

   // Analysis modules
   iAODanalysis    = 0;
   iAODhandler     = 1;
   iESDfilter      = 1;
   iMUONcopyAOD    = 0;
   iJETAN          = 1;
   iPWG4partcorr   = 1;
   iPWG4gammaconv  = 1;
   iPWG4omega3pi   = 1;
   iPWG3vertexing  = 1;
   iPWG3hfe        = 1;
   iPWG3d2h        = 1;
   iPWG2femto      = 1;
   iPWG2spectra    = 1;
     iPWG2protons      = 1;
     iPWG2checkcascade = 1;
     iPWG2perfcascade  = 1;
     iPWG2checkv0      = 1;
     iPWG2strange      = 1;
   iPWG2flow       = 1;
   iPWG2res        = 1;
     iPWG2rsneff   = 1;
   iPWG2kink       = 1;
     iPWG2kinkESDMC    = 1;
     iPWG2kinklikesign = 1;
     iPWG2kinkKstarESD = 1;
     iPWG2kinkKstarMC  = 1;
     iPWG2kinkL1520ESD = 1;
     iPWG2kinkL1520MC  = 1;
     iPWG2kinkPhiESD   = 1;
     iPWG2kinkPhiMC    = 1;
   iPWG2unicor     = 1;
   iPWG2evchar     = 1;
   iPWG2forward    = 1;

// Configuration fot the wagons
   configPWG2femto = "$ALICE_ROOT/PWG2/FEMTOSCOPY/macros/Train/Train3/ConfigFemtoAnalysis.C";
   configPWG3d2h = "$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF_highmult.C";
}
