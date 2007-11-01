void SETUP()
{
   // Load the STEERBase library
   gSystem->Load("libVMC");
   gSystem->Load("libSTEERBase");

   // Set the include paths
   gROOT->ProcessLine(".include STEERBase");

   // Set our location, so that other packages can find us
   gSystem->Setenv("STEERBase_INCLUDE", "STEERBase");
}
