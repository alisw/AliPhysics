void SETUP()
{

   // Load the ESD library
   gSystem->Load("libAnalysisCheck");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IAnalysisCheck");
   gROOT->ProcessLine(".include PWG4");

   // Set our location, so that other packages can find us
   gSystem->Setenv("AnalysisCheck_INCLUDE", "AnalysisCheck");
}
