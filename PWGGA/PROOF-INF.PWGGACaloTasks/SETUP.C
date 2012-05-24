void SETUP()
{

   // Load the ESD library
   //gSystem->Load("libPWGGACaloTasks");
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libPWGGACaloTasks.so");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGGACaloTasks");
   gROOT->ProcessLine(".include PWGGACaloTasks/CaloTasks");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGGACaloTasks_INCLUDE", "PWGGACaloTasks/CaloTasks");
}
