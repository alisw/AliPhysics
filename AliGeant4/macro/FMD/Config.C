void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libFMD");

  AliFMD* FMD = 0;
  switch (version) {
    case 0: FMD =  new AliFMDv0("FMD", "normal FMD");     break;
    case 1: FMD  = new AliFMDv1("FMD", "FMDv1 detector"); break;
  }  

//=================== FMD parameters ============================
}  
