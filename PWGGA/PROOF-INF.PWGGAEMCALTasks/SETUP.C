// $Id$

void SETUP()
{
  // Load the library
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->Load(ocwd+"/libPWGGAEMCALTasks.so");

  // Set the Include paths
  gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGGAEMCALTasks");
  gROOT->ProcessLine(".include PWGGAEMCALTasks/EMCALTasks");
  
  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGGAEMCALTasks_INCLUDE", "PWGGAEMCALTasks/EMCALTasks");
}
