void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWG4Gamma");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4");
   gROOT->ProcessLine(".include PWG4");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4Gamma_INCLUDE", "PWG4");
}
