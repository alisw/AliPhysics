


void SETUP()
{
  // Load some ROOT libraries
  CheckLoadLibrary("libCore");
  CheckLoadLibrary("libTree");
  CheckLoadLibrary("libGeom");
  CheckLoadLibrary("libVMC");
  CheckLoadLibrary("libMinuit");
  
  // Load the AliROOT library
  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libSTEERBase");
  CheckLoadLibrary("libESD");
  CheckLoadLibrary("libAOD");
  CheckLoadLibrary("libANALYSISalice");
  CheckLoadLibrary("libCORRFW");
  CheckLoadLibrary("libPWG3dielectron");
  
  // Set the include paths
  gROOT->ProcessLine(".include PWG3dielectron/dielectron");
    
  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG3dielectron_INCLUDE", "PWG3dielectron");
}


Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library
  
  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
