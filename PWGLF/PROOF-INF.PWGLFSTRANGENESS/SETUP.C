void SETUP() {
  CheckLoadLibrary("libPWGLFSTRANGENESS");

  // Set the include paths
  gROOT->ProcessLine(".include PWGLFSTRANGENESS/STRANGENESS");
  gROOT->ProcessLine(".include PWGLFSTRANGENESS/STRANGENESS/LambdaK0PbPb");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGLFSTRANGENESS_INCLUDE", "PWGLFSTRANGENESS");
}

Int_t CheckLoadLibrary(const char* library) {
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
