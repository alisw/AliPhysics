void SETUP()
{
   // Load the ESD library
   gSystem->Load("libAOD");

   // Set the include paths
   gROOT->ProcessLine(".include AOD");

   // Set our location, so that other packages can find us
   gSystem->Setenv("AOD_INCLUDE", "AOD");
}
