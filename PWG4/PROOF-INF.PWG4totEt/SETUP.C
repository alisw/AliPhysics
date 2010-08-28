void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWG4totEt");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4totEt");
   gROOT->ProcessLine(".include PWG4totEt/totEt");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4totEt_INCLUDE", "PWG4totEt/totEt");
}
