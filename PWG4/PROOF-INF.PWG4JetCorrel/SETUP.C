void SETUP()
{
   gSystem->Load("libPWG4JetCorrel");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4JetCorrel");
   gROOT->ProcessLine(".include PWG4JetCorrel/JetCorrel");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4JetCorrel_INCLUDE", "PWG4JetCorrel/JetCorrel");
}
