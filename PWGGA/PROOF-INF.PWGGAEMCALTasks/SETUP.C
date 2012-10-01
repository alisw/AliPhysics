// $Id$

void SETUP()
{
  // Load the library
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->Load(ocwd+"/libPWGGAEMCALTasks.so");

  // Set the Include paths
 gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGEMCAL -IPWGGAEMCALTasks");
 gROOT->ProcessLine(".include PWGGAEMCAL/EMCAL");
 gROOT->ProcessLine(".include PWGGAEMCALTasks/EMCALTasks");
  
  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGGAEMCALTasks_INCLUDE", "PWGGAEMCALTasks/EMCALTasks");
}
