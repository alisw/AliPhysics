void SETUP()
{
   // Load some ROOT libraries
   CheckLoadLibrary("libEG");
   CheckLoadLibrary("libGeom");

   // Load the ESD library
   CheckLoadLibrary("libESD");

   CheckLoadLibrary("libPWG0base");

   // Set the include paths
   gROOT->ProcessLine(".include PWG0base");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG0base_INCLUDE", "PWG0base");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
