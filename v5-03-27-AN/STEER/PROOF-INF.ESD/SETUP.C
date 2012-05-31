void SETUP()
{
   // Load some ROOT libraries

   gSystem->Load("libVMC");
   gSystem->Load("libNet");
   gSystem->Load("libTree");

   // Load libSTEERBase, ESD depends on it
   gSystem->Load("libSTEERBase");

   // Load the ESD library
   TString dypath = TString::Format("%s:%s", gSystem->WorkingDirectory(), gSystem->GetDynamicPath());
   gSystem->SetDynamicPath(dypath);
   gSystem->Load("libESD");

   // Set the include paths
   gROOT->ProcessLine(".include ESD/ESD");

   // Set our location, so that other packages can find us
   gSystem->Setenv("ESD_INCLUDE", "ESD/ESD");

   // Set our lib coordinates, so that other packages can link to us
   TString lib = TString::Format("-L%s -lESD", gSystem->WorkingDirectory());
   gSystem->Setenv("ESD_LIBS", lib.Data());
}

