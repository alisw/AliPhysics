
{
  printf("\nWELCOME to the magic world of Level3\n\n"); 

  //load library neccesary for TPC alone 
  gSystem->Load("$(ROOTSYS)/lib/libPhysics");
  gSystem->Load("$(ROOTSYS)/lib/libEG");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTEER");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libCONTAINERS");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPC");
  gSystem->Load("/usr/local/franken/lib/MLUC/lib/linux-i386/libMLUC.so");
  gSystem->Load("$(ALICE)/mylibs/libAliL3");
  
  gStyle->SetStatBorderSize(1);
  gStyle->SetTitleBorderSize(0);
  printf("TPC libraries loaded\n");
  printf("Level3 libraries loaded\n");
}
