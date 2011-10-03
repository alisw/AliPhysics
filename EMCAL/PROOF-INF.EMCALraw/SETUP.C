void SETUP()
{

   // Load library
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libEMCALraw.so");

   // Set the Include paths
//   gSystem->SetIncludePath("-I$ROOTSYS/include -IEMCAL");
   gROOT->ProcessLine(".include EMCALraw");

   // Set our location, so that other packages can find us
   gSystem->Setenv("EMCALsim_INCLUDE", "EMCALraw");
}
