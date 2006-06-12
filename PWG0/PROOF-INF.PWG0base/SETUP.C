void SETUP()
{
   gSystem->Load("libPWG0base");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG0base");
   gROOT->ProcessLine(".include PWG0base");
}
