{
   train_name      = "train_LHC09a5";
   proof_cluster   = "alicecaf.cern.ch";
   useAFPAR        = 0;
   proof_dataset   = "/COMMON/COMMON/LHC09a4_run8100X#/esdTree";
   usePLUGIN       = 1;
   usePAR          = 0;
   useCPAR         = 0;
   root_version    = "v20090623";
   aliroot_version = "v4-17-04";
   alien_datadir   = "/alice/sim/PDC_09/LHC09a5/";
   alien_outdir    = "/alice/sim/PDC_09/LHC09a5/AOD";
   maxMergeFiles   = 50;
   mergeExclude    = "AliAOD.root AliAOD.VertexingHF.root AOD.tag.root";
   nRunsPerMaster  = 10;
   nFilesPerJob    = 100;
   run_range[0]    = 90000;
   run_range[1]    = 90040;
   useDBG          = 0;
   useMC           = 1;
   useTAGS         = 0;
   useKFILTER      = 1;
   useTR           = 1;
   useCORRFW       = 0;
   useAODTAGS      = 1;
   saveTrain       = kFALSE;

   // Analysis modules
   iAODanalysis    = 0;
   iAODhandler     = 1;
   iESDfilter      = 1;
   iMUONcopyAOD    = 0;
   iJETAN          = 1;
   iPWG4partcorr   = 1;
   iPWG4gammaconv  = 0;
   iPWG2femto      = 1;
   iPWG2spectra    = 1;
   iPWG2flow       = 1;
   iPWG2res        = 0;
   iPWG2kink       = 1;
   iPWG2unicor     = 1;
   iPWG2evchar     = 1;
}
