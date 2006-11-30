
{
  printf("\nWELCOME to the magic world of Level3\n\n"); 

  gSystem->Load("$(ROOTSYS)/lib/libPhysics");
  gSystem->Load("$(ROOTSYS)/lib/libEG");

  if(0)
    {
      gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTEER");
      gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libCONTAINERS");
      gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPC");
      cout<<"TPC libraries loaded"<<endl;
    }
  if(1)
    {
      gSystem->Load("$(LEVEL3)/kip/MLUC/lib/linux-i386/libMLUC.so");
      gSystem->Load("$(LEVEL3)/lib_$(USERNAME)/libAliHLT");
      gSystem->Load("$(LEVEL3)/lib_$(USERNAME)/libAliHLTHough");
      gSystem->Load("$(LEVEL3)/lib_$(USERNAME)/libAliHLTComp");
      cout<<"HLT libraries loaded"<<endl;
    }
  gROOT->LoadMacro("$(HOME)/alirootcode/XFunct.C");
  gStyle->SetStatBorderSize(1);
  gStyle->SetTitleBorderSize(0);

}
