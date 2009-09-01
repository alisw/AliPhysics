void SETUP()
{
  // Load the ANALYSIS library
   gSystem->Load("libTENDER");

   // Set the include paths
   gROOT->ProcessLine(".include TENDER/Tender");

   // Set our location, so that other packages can find us
   gSystem->Setenv("Tender_INCLUDE", "TENDER/Tender");
}
