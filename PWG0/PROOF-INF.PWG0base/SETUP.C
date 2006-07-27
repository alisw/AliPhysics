void SETUP()
{
   // Load some ROOT libraries
   gSystem->Load("libEG");
   gSystem->Load("libGeom");

   // Load the ESD library
   gSystem->Load("libESD");

   gSystem->Load("libPWG0base");

   // Set the include paths
   gROOT->ProcessLine(".include PWG0base");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG0base_INCLUDE", "PWG0base");
}
