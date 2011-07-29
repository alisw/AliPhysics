void SETUP()
{
   // Load some ROOT libraries

   gSystem->Load("libVMC");
   gSystem->Load("libNet");
   gSystem->Load("libTree");

   // Load libSTEERBase, AOD depends on it
   gSystem->Load("libSTEERBase");

   // Load the AOD library
   gSystem->Load("libAOD");

   // Set the include paths
   gROOT->ProcessLine(".include AOD");

   // Set our location, so that other packages can find us
   gSystem->Setenv("AOD_INCLUDE", "AOD/AOD");
}
