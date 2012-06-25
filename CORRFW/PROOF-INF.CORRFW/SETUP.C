void SETUP()
{
  // Load the CF library
   gSystem->Load("libCORRFW");

   // Set the include paths
   gROOT->ProcessLine(".include CORRFW");

   // Set our location, so that other packages can find us
   gSystem->Setenv("CORRFW_INCLUDE", "CORRFW");

   // Set our lib coordinates, so that other packages can link to us
   TString lib = TString::Format("-L%s -lCORRFW", gSystem->WorkingDirectory());
   gSystem->Setenv("CORRFW_LIBS", lib.Data());
}
