void SETUP()
{

   // Load the JET-Tasks library
   gSystem->Load("libPWG4JetTasks");


   gSystem->AddIncludePath("-I$ROOTSYS/include -IJETAN");
   gROOT->ProcessLine(".include JETAN");

   // Set the Include paths
   gSystem->AddIncludePath("-IPWG4");
   gROOT->ProcessLine(".include PWG4");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4JetTasks_INCLUDE", "PWG4");
}
