void SETUP()
{
   // we assume PWG0base (and thus ESD) already loaded
   CheckLoadLibrary("libMinuit");

   // this package depends on STEER
   CheckLoadLibrary("libVMC");
   CheckLoadLibrary("libMinuit");
   CheckLoadLibrary("libSTEER");

   // more packages to access the alice event header
   CheckLoadLibrary("libEVGEN");
   CheckLoadLibrary("libFASTSIM");
   CheckLoadLibrary("libmicrocern");
   CheckLoadLibrary("libpdf");    // both are needed because old aliroot needs pdf, new one lhapdf
   CheckLoadLibrary("liblhapdf"); // 
   CheckLoadLibrary("libpythia6");
   CheckLoadLibrary("libEGPythia6");
   CheckLoadLibrary("libAliPythia6");

   CheckLoadLibrary("libPWG0dep");

   // Set the Include paths
   gROOT->ProcessLine(".include PWG0dep");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
