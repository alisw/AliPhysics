void SETUP()
{
   // Load some ROOT libraries
   gSystem->Load("libFASTJETAN");

   // Set the Inlucde paths
   gROOT->ProcessLine(".include FASTJETAN");

   // Set our location, so that other packages can find us
   gSystem->Setenv("FASTJETAN_INCLUDE", "FASTJETAN");
}
