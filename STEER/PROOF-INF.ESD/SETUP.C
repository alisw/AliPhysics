void SETUP()
{
   // Load the ESD library
   gSystem->Load("libANALYSIS_NEW");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IANALYSIS_NEW");
   gROOT->ProcessLine(".include ANALYSIS_NEW");

   // Set our location, so that other packages can find us
   gSystem->Setenv("ANALYSIS_NEW_INCLUDE", "ANALYSIS_NEW");
}
