{
  printf("\nWELCOME to ALICE\n\n");
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/install/include -I$ALICE_ROOT/STEER\
   -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW  -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/TPCBase  -I$ALICE_ROOT/TPC/TPCRec -I$ALICE_ROOT/TPC/TPCCalib -I$ALICE_PHYSICS/../src/PWGPP/TPC/");
  //load library neccesary for TPC alone
  //gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libANALYSIS");
  //gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libANALYSISalice");
  //gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPCcalib");


}

