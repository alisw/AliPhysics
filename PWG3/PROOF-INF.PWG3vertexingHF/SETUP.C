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
    CheckLoadLibrary("libCORRFW");
    CheckLoadLibrary("libANALYSISalice");
    CheckLoadLibrary("libPWG3vertexingHF");

   // Set the include paths
   gROOT->ProcessLine(".include PWG3vertexingHF");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG3base_INCLUDE", "PWG3vertexingHF");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
