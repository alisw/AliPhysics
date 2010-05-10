void SETUP()
{

   // Load library
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libEMCALUtils.so");

   // Set the Include paths
//   gSystem->SetIncludePath("-I$ROOTSYS/include -IEMCAL");
   gROOT->ProcessLine(".include EMCALUtils");

   // Set our location, so that other packages can find us
   gSystem->Setenv("EMCALUtils_INCLUDE", "EMCALUtils");
}
