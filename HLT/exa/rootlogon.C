
{
  printf("\nWELCOME to the magic world of Level3\n\n"); 

  gSystem->Load("$(ROOTSYS)/lib/libPhysics");
  gSystem->Load("$(ROOTSYS)/lib/libEG");

  if(1)
    {
      if(getenv("MLUCDIR")) {
        gSystem->Load("$(ALICE_ROOT)/lib/libSTEER");
        gSystem->Load("$(ALICE_ROOT)/lib/libCONTAINERS");
        gSystem->Load("$(ALICE_ROOT)/lib/libTPC");
      } else {
        gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTEER");
        gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libCONTAINERS");
        gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPC");
      }
      cout<<"TPC libraries loaded"<<endl;
    }
  if(0)
    {
      if(getenv("MLUCDIR")) {
        gSystem->Load("$(MLUCDIR/lib/linux-i386/libMLUC");
        gSystem->Load("$(TOPDIR)/lib_ALIROOT/libAliL3");
        gSystem->Load("$(TOPDIR)/lib_ALIROOT/libAliL3Hough");
        gSystem->Load("$(TOPDIR)/lib_ALIROOT/libAliL3Comp");
      } else {
        gSystem->Load("$(LEVEL3)/kip/MLUC/lib/linux-i386/libMLUC.so");
        gSystem->Load("$(LEVEL3)/lib_$(USERNAME)/libAliL3");
        gSystem->Load("$(LEVEL3)/lib_$(USERNAME)/libAliL3Hough");
        gSystem->Load("$(LEVEL3)/lib_$(USERNAME)/libAliL3Comp");
      }
      cout<<"HLT libraries loaded"<<endl;
    }

}
