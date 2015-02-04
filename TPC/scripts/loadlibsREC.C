/// \file loadlibsREC.C

void loadlibsREC ()
{
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER  -I$ALICE_ROOT/TPC -I$ALICE_ROOT/TPC/CalibMacros -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD ");
  gSystem->Load("libCore");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libGui");

  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libNet");
  gSystem->Load("libTree");
  gSystem->Load("libProof");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSIS");
  
  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCrec");
  gSystem->Load("libTPCcalib");
  
  gSystem->Load("libSTAT");
  
  gSystem->Load("libThread");
}
