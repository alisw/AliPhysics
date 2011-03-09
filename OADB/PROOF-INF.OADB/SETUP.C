Int_t SETUP()
{
  // load library
  gSystem->Load("libOADB");

  // set include path
  gROOT->ProcessLine(".include OADB");

  // Set our location, so that other packages can find us
  gSystem->Setenv("OADB_INCLUDE", "OADB");
  
  // path for the root files
  gSystem->Setenv("OADB_PATH", "OADB");

  // We're happy
  return 0;

}
