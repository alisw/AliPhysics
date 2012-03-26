void SETUP()
{
  // Load the ANALYSIS library
   gSystem->Load("libANALYSISalice");

   // Set the include paths
   gROOT->ProcessLine(".include ANALYSISalice");

   // Set our location, so that other packages can find us
   gSystem->Setenv("ANALYSISalice_INCLUDE", "ANALYSISalice");

   // Set our lib coordinates, so that other packages can link to us
   TString lib = TString::Format("-L%s -lANALYSISalice", gSystem->WorkingDirectory());
   gSystem->Setenv("ANALYSISalice_LIBS", lib.Data());
}
