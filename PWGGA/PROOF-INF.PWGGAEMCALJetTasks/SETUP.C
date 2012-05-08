// $Id$

void SETUP()
{
  // Load the library
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->Load(ocwd+"/libPWGGAEMCALJetTasks.so");

  // Set the Include paths
  gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGGAEMCALJetTasks");
  gROOT->ProcessLine(".include PWGGAEMCALTasks/EMCALJetTasks");
  
  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGGAEMCALJetTasks_INCLUDE", "PWGGAEMCALTasks/EMCALJetTasks");
}
