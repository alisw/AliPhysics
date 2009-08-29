void SETUP()
{

   // Load library
   gSystem->Load("libEMCALrec");

   // Set the Include paths
//   gSystem->SetIncludePath("-I$ROOTSYS/include -IEMCAL");
   gROOT->ProcessLine(".include EMCALrec");

   // Set our location, so that other packages can find us
   gSystem->Setenv("EMCALrec_INCLUDE", "EMCALrec");
}
