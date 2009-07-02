void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWG4omega3pi");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4omega3pi");
   gROOT->ProcessLine(".include PWG4omega3pi/omega3pi");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4omega3pi_INCLUDE", "PWG4omega3pi/omega3pi");
}
