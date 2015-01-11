// $Id: SETUP.C 54472 2012-02-07 12:21:23Z loizides $

void SETUP()
{
  // Load the library
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->Load(ocwd+"/libPWGEMCAL");

  // Set the Include paths
  gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGEMCAL");
  gROOT->ProcessLine(".include PWGEMCAL/EMCAL");
  
  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGEMCAL_INCLUDE", "PWGEMCAL/EMCAL");
}
