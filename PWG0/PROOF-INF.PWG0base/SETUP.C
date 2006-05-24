void SETUP()
{
   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG0base");
   gROOT->ProcessLine(".include PWG0base");
}
