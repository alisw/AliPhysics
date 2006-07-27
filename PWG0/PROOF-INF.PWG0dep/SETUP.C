void SETUP()
{
   // we assume PWG0base (and thus ESD) already loaded

   // this package depends on STEER
   gSystem->Load("libVMC");
   gSystem->Load("libMinuit");
   gSystem->Load("libSTEER");

   gSystem->Load("libPWG0dep");

   // Set the Include paths
   gROOT->ProcessLine(".include PWG0dep");
}
