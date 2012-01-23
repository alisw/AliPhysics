void SETUP()
{

   // Load the ESD library
   //gSystem->Load("libPWGGACaloTrackCorrBase");
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libPWGGACaloTrackCorrBase.so");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGGACaloTrackCorrBase");
   gROOT->ProcessLine(".include PWGGACaloTrackCorrBase/CaloTrackCorrBase");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGGACaloTrackCorrBase_INCLUDE", "PWGGACaloTrackCorrBase/CaloTrackCorrBase");
}
