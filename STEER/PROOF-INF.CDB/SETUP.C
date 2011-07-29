void SETUP()
{
   // Load some ROOT libraries

   gSystem->Load("libVMC");
   gSystem->Load("libNet");
   gSystem->Load("libTree");

   // Load libSTEERBase, CDB depends on it
   gSystem->Load("libSTEERBase");

   // Load the CDB library
   gSystem->Load("libCDB");

   // Set the include paths
   gROOT->ProcessLine(".include CDB");

   // Set our location, so that other packages can find us
   gSystem->Setenv("CDB_INCLUDE", "CDB/CDB");
}
