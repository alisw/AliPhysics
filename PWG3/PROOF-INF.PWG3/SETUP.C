void SETUP()
{
    // Load some ROOT libraries
    CheckLoadLibrary("libTree");
    CheckLoadLibrary("libGeom");
    CheckLoadLibrary("libVMC");

    // Load the ESD library
    CheckLoadLibrary("libANALYSIS");
    CheckLoadLibrary("libSTEERBase");
    CheckLoadLibrary("libESD");
    CheckLoadLibrary("libAOD");
    CheckLoadLibrary("libANALYSISalice");
    CheckLoadLibrary("libPWG3");


   // Set the include paths
   gROOT->ProcessLine(".include PWG3");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG3_INCLUDE", "PWG3");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
