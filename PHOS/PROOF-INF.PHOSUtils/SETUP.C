void SETUP()
{

   // Load library
   gSystem->Load("libPHOSUtils");

   // Set the Include paths
//   gSystem->SetIncludePath("-I$ROOTSYS/include -IPHOS");
   gROOT->ProcessLine(".include PHOSUtils");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PHOSUtils_INCLUDE", "PHOSUtils");
}
