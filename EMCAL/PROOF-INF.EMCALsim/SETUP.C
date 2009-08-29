void SETUP()
{

   // Load library
   gSystem->Load("libEMCALsim");

   // Set the Include paths
//   gSystem->SetIncludePath("-I$ROOTSYS/include -IEMCAL");
   gROOT->ProcessLine(".include EMCALsim");

   // Set our location, so that other packages can find us
   gSystem->Setenv("EMCALsim_INCLUDE", "EMCALsim");
}
