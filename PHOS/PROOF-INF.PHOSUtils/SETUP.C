void SETUP()
{

   // Load library
   gSystem->Load("libPPHOSUtils");

   // Set the Include paths
//   gSystem->SetIncludePath("-I$ROOTSYS/include -IPHOS");
   gROOT->ProcessLine(".include PHOS");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PHOSUtils_INCLUDE", "PHOS");
}
