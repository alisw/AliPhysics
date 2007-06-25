void SETUP()
{
  
  // load all libraries
  //gProof->Exec("gROOT->Macro(\"$ALICE_ROOT/macros/loadlibsREC.C\")");
   // Load the ESD library
   gSystem->Load("libTPCcalib");

   // Set the Inlucde paths
   //gSystem->SetIncludePath("-I$ROOTSYS/include -ITPCcalib");
   //gROOT->ProcessLine(".include TPCcalib");

   // Set our location, so that other packages can find us
   //gSystem->Setenv("TPCcalib_INCLUDE", "TPCcalib");
}
