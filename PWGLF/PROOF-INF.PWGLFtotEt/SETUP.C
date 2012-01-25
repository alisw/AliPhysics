void SETUP()
{

   // Load the ESD library
   gSystem->Load("libPWGLFtotEt");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGLFtotEt");
   gROOT->ProcessLine(".include PWGLFtotEt/totEt");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGLFtotEt_INCLUDE", "PWGLFtotEt/totEt");
}
