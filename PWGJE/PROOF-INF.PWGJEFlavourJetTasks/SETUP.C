// $Id: SETUP.C 58098 2012-08-06 04:23:06Z loizides $

void SETUP()
{
  // Load the library
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->Load(ocwd+"/libPWGJEFlavourJetTasks.so");

  // Set the Include paths
  gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGJEFlavourJetTasks");
  gROOT->ProcessLine(".include PWGJEFlavourTasks/FlavourJetTasks");
  
  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGJEFlavourJetTasks_INCLUDE", "PWGJEFlavourTasks/FlavourJetTasks");
}
