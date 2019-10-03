
{
  printf("\nWELCOME to ALICE\n\n"); 
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER\
   -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW -I$ALICE_ROOT/kdtree/include");
  //load library neccesary for TPC alone 
  gSystem->Load("$(ROOTSYS)/lib/libEG");
  gSystem->Load("$(ROOTSYS)/lib/libGeom");  
  gSystem->Load("$(ROOTSYS)/lib/libVMC");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTEER");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPCrec");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPCsim");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTRDrec");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTRDsim");
  //gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTRD");
  //gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTRACKING");
  
  
}

