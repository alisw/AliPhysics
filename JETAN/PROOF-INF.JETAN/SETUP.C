void SETUP()
{
   // Load some ROOT libraries
   gSystem->Load("libEG");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC")
   gSystem->Load("libESD");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IJETAN");
   gROOT->ProcessLine(".include JETAN");

   // Set our location, so that other packages can find us
   gSystem->Setenv("JETAN_INCLUDE", "JETAN");
}
