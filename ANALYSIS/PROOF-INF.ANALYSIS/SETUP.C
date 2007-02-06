void SETUP()
{
  // Load the ANALYSIS library
   gSystem->Load("libANALYSIS");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IANALYSIS");
   gROOT->ProcessLine(".include ANALYSIS");

   // Set our location, so that other packages can find us
   gSystem->Setenv("ANALYSIS_INCLUDE", "ANALYSIS");
}
