void SETUP()
{
   // Load some ROOT libraries

   gSystem->Load("libVMC");
   gSystem->Load("libNet");
   gSystem->Load("libTree");

   // Load libSTEERBase, ESD depends on it
   gSystem->Load("libSTEERBase");

   // Load the ESD library
   gSystem->Load("libESD");

   // Set the include paths
   gROOT->ProcessLine(".include ESD");

   // Set our location, so that other packages can find us
   gSystem->Setenv("ESD_INCLUDE", "ESD/ESD");
}
