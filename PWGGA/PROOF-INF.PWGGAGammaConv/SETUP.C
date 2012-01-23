void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWGGAGammaConv");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGGAGammaConv");
   gROOT->ProcessLine(".include PWGGAGammaConv/GammaConv");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGGAGammaConv_INCLUDE", "PWGGAGammaConv/GammaConv");
}
