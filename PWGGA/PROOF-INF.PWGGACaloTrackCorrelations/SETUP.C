void SETUP()
{

   // Load the ESD library
   //gSystem->Load("libPWGGACaloTrackCorrelations");
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libPWGGACaloTrackCorrelations.so");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGCaloTrackCorrBase  -IPWGGACaloTrackCorrelations");
   gROOT->ProcessLine(".include PWGCaloTrackCorrBase/CaloTrackCorrBase");
   gROOT->ProcessLine(".include PWGGACaloTrackCorrelations/CaloTrackCorrelations");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGGACaloTrackCorrelations_INCLUDE", "PWGGACaloTrackCorrelations/CaloTrackCorrelations");
}
