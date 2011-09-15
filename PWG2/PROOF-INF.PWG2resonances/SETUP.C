void SETUP()
{

   gSystem->SetDynamicPath(Form("%s:%s", gSystem->pwd(), gSystem->GetDynamicPath()));
   CheckLoadLibrary("libPWG2resonances");

   gROOT->ProcessLine(".include PWG2resonances");
   gROOT->ProcessLine(".include PWG2resonances/RESONANCES");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWG2resonances_INCLUDE", "PWG2resonances/RESONANCES");
}

Int_t CheckLoadLibrary(const char* library)
{
   // checks if a library is already loaded, if not loads the library

   if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
      return 1;

   return gSystem->Load(library);
}
