void SETUP()
{

   // Load the TRD PWG library
   gSystem->Load("libPWGTRD");

   // Set the Include paths
   gSystem->AddIncludePath("-IPWG/TRD");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGTRD_INCLUDE", "PWG/TRD");
}
