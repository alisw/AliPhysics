void SETUP()
{

   // Load the JET-Tasks library
   gSystem->Load("libPWGTools");

   // Set the Include paths
   gSystem->AddIncludePath("-IPWGTools");
   gROOT->ProcessLine(".include PWGTools/Tools");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGTools_INCLUDE", "PWGTools/Tools");
}
