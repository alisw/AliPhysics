// $Id$

void SETUP()
{
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libPWG4UserTasks.so");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4UserTasks");
   gROOT->ProcessLine(".include PWG4UserTasks/UserTasks");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4UserTasks_INCLUDE", "PWG4UserTasks/UserTasks");
}
