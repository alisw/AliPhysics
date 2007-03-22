void SETUP()
{
  // Load the ANALYSIS library
   gSystem->Load("libANALYSISRL");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IANALYSISRL");
   gROOT->ProcessLine(".include ANALYSISRL");

   // Set our location, so that other packages can find us
   gSystem->Setenv("ANALYSISRL_INCLUDE", "ANALYSISRL");
}
