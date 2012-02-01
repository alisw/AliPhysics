void SETUP()
{

   gSystem->SetDynamicPath(Form("%s:%s", gSystem->pwd(), gSystem->GetDynamicPath()));
   CheckLoadLibrary("libPWGLFresonances");

   gROOT->ProcessLine(".include PWGLFresonances");
   gROOT->ProcessLine(".include PWGLFresonances/RESONANCES");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGLFresonances_INCLUDE", "PWGLFresonances/RESONANCES");
}

Int_t CheckLoadLibrary(const char* library)
{
   // checks if a library is already loaded, if not loads the library

   if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
      return 1;

   return gSystem->Load(library);
}
