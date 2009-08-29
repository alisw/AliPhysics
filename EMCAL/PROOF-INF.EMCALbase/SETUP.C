void SETUP()
{

   // Load library
   gSystem->Load("libEMCALbase");

   // Set the Include paths
//   gSystem->SetIncludePath("-I$ROOTSYS/include -IEMCAL");
   gROOT->ProcessLine(".include EMCALbase");

   // Set our location, so that other packages can find us
   gSystem->Setenv("EMCALbase_INCLUDE", "EMCALbase");
}
