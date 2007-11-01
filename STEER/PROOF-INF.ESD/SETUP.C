void SETUP()
{
   // Load some ROOT libraries
   gSystem->Load("libEG");
   gSystem->Load("libGeom");
   gSystem->Load("libProof");

   // Load the ESD library
   gSystem->Load("libESD");

   // Set the include paths
   gROOT->ProcessLine(".include ESD");

   // Set our location, so that other packages can find us
   gSystem->Setenv("ESD_INCLUDE", "ESD");
}
