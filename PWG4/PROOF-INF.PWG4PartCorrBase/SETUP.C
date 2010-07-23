void SETUP()
{

   // Load the ESD library
   //gSystem->Load("libPWG4PartCorrBase");
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libPWG4PartCorrBase.so");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4PartCorrBase");
   gROOT->ProcessLine(".include PWG4PartCorrBase/PartCorrBase");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4PartCorrBase_INCLUDE", "PWG4PartCorrBase/PartCorrBase");
}
