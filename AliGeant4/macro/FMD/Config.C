void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libFMD");

  AliFMD* FMD = 0;
  switch (version) {
    case 0: FMD =  new AliFMDv0("FMD", "FMDv0 detector"); break;
    case 1: FMD  = new AliFMDv1("FMD","normal FMD");        break;
  }  

//=================== FMD parameters ============================
}  
