void SETUP()
{
   // Load some ROOT libraries
   gSystem->Load("libJETANdev");

   // Set the Inlucde paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IJETAN");
   gROOT->ProcessLine(".include JETAN");

   // Set our location, so that other packages can find us
   gSystem->Setenv("JETAN_INCLUDE", "JETAN");
}
