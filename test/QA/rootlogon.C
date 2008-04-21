{
  printf("\nWELCOME to ALICE QA\n\n"); 
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TRD -I$ALICE_ROOT/RAW") ; 
  //load library neccesary for TPC alone 
//  gSystem->Load("$(ROOTSYS)/lib/libEG");
//  gSystem->Load("$(ROOTSYS)/lib/libGeom");  
//  gSystem->Load("$(ROOTSYS)/lib/libVMC");
  //gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTRACKING");
  
  
}

