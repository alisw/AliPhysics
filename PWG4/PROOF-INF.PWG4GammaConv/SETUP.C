void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWG4GammaConv");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4GammaConv");
   gROOT->ProcessLine(".include PWG4GammaConv/GammaConv");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4GammaConv_INCLUDE", "PWG4GammaConv/GammaConv");
}
