void SETUP()
{
   // Load some ROOT libraries
   CheckLoadLibrary("libEG");
   CheckLoadLibrary("libGeom");

   // Load the ESD library
  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libESD");
  CheckLoadLibrary("libPWGUDbase");


   // Set the include paths
   gROOT->ProcessLine(".include PWGUDbase/base");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGUDbase_INCLUDE", "PWGUDbase/base");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
