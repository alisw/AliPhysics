void SETUP()
{
    // Load some ROOT libraries
    CheckLoadLibrary("libTree");
    CheckLoadLibrary("libGeom");
    CheckLoadLibrary("libVMC");
    CheckLoadLibrary("libMinuit");

    // Load the ESD library
    CheckLoadLibrary("libANALYSIS");
    CheckLoadLibrary("libSTEERBase");
    CheckLoadLibrary("libESD");
    CheckLoadLibrary("libAOD");
    CheckLoadLibrary("libCORRFW");
    CheckLoadLibrary("libANALYSISalice");
    CheckLoadLibrary("libPWGHFbase");
    CheckLoadLibrary("libPWGHFvertexingHF");
    CheckLoadLibrary("libPWGHFhfe");
    CheckLoadLibrary("libPWGHFcorrelationHF");

   // Set the include paths
   gROOT->ProcessLine(".include PWGHFcorrelationHF/correlationHF");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGHFcorrelationHF_INCLUDE", "PWGHFcorrelationHF/correlationHF");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(library, "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
