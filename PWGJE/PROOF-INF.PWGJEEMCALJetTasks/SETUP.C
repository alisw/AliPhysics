// $Id$

void SETUP()
{
  // Load the library
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->Load(ocwd+"/libPWGJEEMCALJetTasks");

  // Set the Include paths
  gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGJEEMCALJetTasks");
  gROOT->ProcessLine(".include PWGJEEMCALTasks/EMCALJetTasks");
  
  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGJEEMCALJetTasks_INCLUDE", "PWGJEEMCALTasks/EMCALJetTasks");
}
