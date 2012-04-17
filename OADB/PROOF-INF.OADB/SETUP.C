Int_t SETUP()
{
  // load library
  gSystem->Load("libOADB");

  // set include path
  gROOT->ProcessLine(".include OADB");

  // Set our location, so that other packages can find us
  gSystem->Setenv("OADB_INCLUDE", "OADB");
  
  // path for the root files
  const char* oadbPath = gSystem->Getenv("OADB_PATH");
  if (!oadbPath || oadbPath[0] == '\0')
    gSystem->Setenv("OADB_PATH", "OADB");

  // Set our lib coordinates, so that other packages can link to us
  TString lib = TString::Format("-L%s -lOADB", gSystem->WorkingDirectory());
  gSystem->Setenv("OADB_LIBS", lib.Data());

  // We're happy
  return 0;

}
