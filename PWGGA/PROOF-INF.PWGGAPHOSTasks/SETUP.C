// $Id$

void SETUP()
{
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libPWGGAPHOSTasks.so");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGGAPHOSTasks");
   gROOT->ProcessLine(".include PWGGAPHOSTasks/PHOSTasks");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGGAPHOSTasks_INCLUDE", "PWGGAPHOSTasks/PHOSTasks");
}
