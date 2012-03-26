void SETUP()
{
  // Load the ANALYSIS library
   gSystem->Load("libANALYSIS");

   // Set the include paths
   gROOT->ProcessLine(".include ANALYSIS");

   // Set our location, so that other packages can find us
   gSystem->Setenv("ANALYSIS_INCLUDE", "ANALYSIS");

   // Set our lib coordinates, so that other packages can link to us
   TString lib = TString::Format("-L%s -lANALYSIS", gSystem->WorkingDirectory());
   gSystem->Setenv("ANALYSIS_LIBS", lib.Data());
}
