void SETUP()
{

   // Load library
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libEMCALbase.so");

   // Set the Include paths
//   gSystem->SetIncludePath("-I$ROOTSYS/include -IEMCAL");
   gROOT->ProcessLine(".include EMCALbase");

   // Set our location, so that other packages can find us
   gSystem->Setenv("EMCALbase_INCLUDE", "EMCALbase");
}
