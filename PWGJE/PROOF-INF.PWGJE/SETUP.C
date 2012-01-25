void SETUP()
{

   // Load the JET-Tasks library
   gSystem->Load("libPWGJE");


   gSystem->AddIncludePath("-I$ROOTSYS/include -IJETAN");
   gROOT->ProcessLine(".include JETAN");

   // Set the Include paths
   gSystem->AddIncludePath("-IPWGJE");
   gROOT->ProcessLine(".include PWGJE/");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGJE_INCLUDE", "PWGJE");
}
