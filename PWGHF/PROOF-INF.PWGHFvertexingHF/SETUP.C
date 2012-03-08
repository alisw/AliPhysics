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
    CheckLoadLibrary("libPWGflowBase");
    CheckLoadLibrary("libPWGflowTasks");
    CheckLoadLibrary("libPWGHFbase");
    CheckLoadLibrary("libPWGHFvertexingHF");

   // Set the include paths
   gROOT->ProcessLine(".include PWGHFvertexingHF/vertexingHF");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGHFvertexingHF_INCLUDE", "PWGHFvertexingHF/vertexingHF");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
