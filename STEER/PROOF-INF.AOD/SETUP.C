void SETUP()
{
   // Load the ESD library
   gSystem->Load("libAOD");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IAOD");
   gROOT->ProcessLine(".include AOD");

   // Set our location, so that other packages can find us
   gSystem->Setenv("AOD_INCLUDE", "AOD");
}
