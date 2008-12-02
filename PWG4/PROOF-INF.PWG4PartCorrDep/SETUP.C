void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWG4PartCorrDep");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4PartCorrDep");
   gROOT->ProcessLine(".include PWG4PartCorrDep/PartCorrDep");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4PartCorrDep_INCLUDE", "PWG4PartCorrDep/PartCorrDep");
}
