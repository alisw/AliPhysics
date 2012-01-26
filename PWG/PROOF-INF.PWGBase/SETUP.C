void SETUP()
{

   // Load the JET-Tasks library
   gSystem->Load("libPWGBase");

   // Set the Include paths
   gSystem->AddIncludePath("-IPWGBase");
   gROOT->ProcessLine(".include PWG4Base/Base");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGBase_INCLUDE", "PWGBase/Base");
}
