void SETUP()
{
    // Load some ROOT libraries
    CheckLoadLibrary("libTree");
    CheckLoadLibrary("libGeom");
    CheckLoadLibrary("libVMC");
    CheckLoadLibrary("libMinuit");
    CheckLoadLibrary("libProofPlayer");

    // Load the ESD libraries
    CheckLoadLibrary("libANALYSIS");
    CheckLoadLibrary("libSTEERBase");
    CheckLoadLibrary("libESD");
    CheckLoadLibrary("libAOD");

  // Load the STEER libraries
    CheckLoadLibrary("libRAWDatabase");
    CheckLoadLibrary("libRAWDatarec");
    CheckLoadLibrary("libCDB");
    CheckLoadLibrary("libSTEER");
    CheckLoadLibrary("libPhysics");
    CheckLoadLibrary("libANALYSISalice");

    // Load the MUON libraries
    CheckLoadLibrary("libMUONcore");
    CheckLoadLibrary("libMUONmapping");
    CheckLoadLibrary("libMUONraw");
    CheckLoadLibrary("libMUONcalib");
    CheckLoadLibrary("libMUONgeometry");
    CheckLoadLibrary("libMUONtrigger");
    CheckLoadLibrary("libMUONbase");
    CheckLoadLibrary("libMUONrec");

    CheckLoadLibrary("libPWGmuondep");

   // Set the include paths
   gROOT->ProcessLine(".include PWGmuondep/muondep");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGmuondep_INCLUDE", "PWGmuondep/muondep");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
