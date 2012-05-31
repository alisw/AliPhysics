void SETUP()
{
   // Load some ROOT libraries

   gSystem->Load("libVMC");
   gSystem->Load("libNet");
   gSystem->Load("libTree");

   // Load libSTEERBase, AOD depends on it
   gSystem->Load("libSTEERBase");

   // Load the AOD library
   TString dypath = TString::Format("%s:%s", gSystem->WorkingDirectory(), gSystem->GetDynamicPath());
   gSystem->SetDynamicPath(dypath);
   gSystem->Load("libAOD");

   // Set the include paths
   gROOT->ProcessLine(".include AOD/AOD");

   // Set our location, so that other packages can find us
   gSystem->Setenv("AOD_INCLUDE", "AOD/AOD");

   // Set our lib coordinates, so that other packages can link to us
   TString lib = TString::Format("-L%s -lAOD", gSystem->WorkingDirectory());
   gSystem->Setenv("AOD_LIBS", lib.Data());
}

