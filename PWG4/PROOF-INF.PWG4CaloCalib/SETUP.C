void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWG4CaloCalib");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWG4CaloCalib");
   gROOT->ProcessLine(".include PWG4CaloCalib/CaloCalib");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG4CaloCalib_INCLUDE", "PWG4CaloCalib/CaloCalib");
}
