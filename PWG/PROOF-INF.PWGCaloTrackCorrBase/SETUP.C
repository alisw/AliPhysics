void SETUP()
{

   // Load the ESD library
   //gSystem->Load("libPWGCaloTrackCorrBase");
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libPWGCaloTrackCorrBase");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGCaloTrackCorrBase");
   gROOT->ProcessLine(".include PWGCaloTrackCorrBase/CaloTrackCorrBase");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGCaloTrackCorrBase_INCLUDE", "PWGCaloTrackCorrBase/CaloTrackCorrBase");
}
