void SETUP()
{
   // Load some Root libraries needed by STEERBase
   gSystem->Load("libVMC");
   gSystem->Load("libNet");
   gSystem->Load("libTree");
   gSystem->Load("libPhysics");
   gSystem->Load("libMinuit");

   // Load the STEERBase library
   TString dypath = TString::Format("%s:%s", gSystem->WorkingDirectory(), gSystem->GetDynamicPath());
   gSystem->SetDynamicPath(dypath);
   gSystem->Load("libSTEERBase");

   // Set the include paths
   gROOT->ProcessLine(".include STEERBase/STEERBase");

   // Set our location, so that other packages can find us
   gSystem->Setenv("STEERBase_INCLUDE", "STEERBase/STEERBase");

   // Set our lib coordinates, so that other packages can link to us
   TString lib = TString::Format("-L%s -lSTEERBase", gSystem->WorkingDirectory());
   gSystem->Setenv("STEERBase_LIBS", lib.Data());
}
