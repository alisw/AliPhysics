void SETUP()
{
   // Load some ROOT libraries

   gSystem->Load("libVMC");
   gSystem->Load("libNet");
   gSystem->Load("libTree");

   // Load libSTEERBase, CDB depends on it
   gSystem->Load("libSTEERBase");

   // Load the CDB library
   TString dypath = TString::Format("%s:%s", gSystem->WorkingDirectory(), gSystem->GetDynamicPath());
   gSystem->SetDynamicPath(dypath);
   gSystem->Load("libCDB");

   // Set the include paths
   gROOT->ProcessLine(".include CDB/CDB");

   // Set our location, so that other packages can find us
   gSystem->Setenv("CDB_INCLUDE", "CDB/CDB");

   // Set our lib coordinates, so that other packages can link to us
   TString lib = TString::Format("-L%s -lCDB", gSystem->WorkingDirectory());
   gSystem->Setenv("CDB_LIBS", lib.Data());
}
