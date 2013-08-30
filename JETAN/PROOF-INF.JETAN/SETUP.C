void SETUP()
{
   // Load some ROOT libraries
   gSystem->Load("libJETANdev");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IJETANdev");
   gROOT->ProcessLine(".include JETANdev/DEV");

   // Set our location, so that other packages can find us
   gSystem->Setenv("JETANDEV_INCLUDE", "JETANdev/DEV/");
}
