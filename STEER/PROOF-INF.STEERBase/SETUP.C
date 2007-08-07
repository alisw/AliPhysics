void SETUP()
{
   // Load the ESD library
   gSystem->Load("libSTEERBase");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -ISTEERBase");
   gROOT->ProcessLine(".include STEERBase");

   // Set our location, so that other packages can find us
   gSystem->Setenv("STEERBase_INCLUDE", "STEERBase");
}
