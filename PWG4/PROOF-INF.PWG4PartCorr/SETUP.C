void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWG4PartCorr");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4PartCorr");
   gROOT->ProcessLine(".include PWG4PartCorr/PartCorr");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4PartCorr_INCLUDE", "PWG4PartCorr/PartCorr");
}
