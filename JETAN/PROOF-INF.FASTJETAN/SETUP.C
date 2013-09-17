void SETUP()
{
   // Load some ROOT libraries
   gSystem->Load("libFASTJETANdev");

   // Set the Inlucde paths
   gROOT->ProcessLine(".include FASTJETANdev/DEV");

   // Set our location, so that other packages can find us
   gSystem->Setenv("FASTJETANDEV_INCLUDE", "FASTJETANdev/DEV");
}
