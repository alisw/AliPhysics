void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libPMD");

  AliPMD *PMD = 0;
  switch (version) {
    case 0: PMD  = new AliPMDv0("PMD","normal PMD");      break;
    case 1: PMD  = new AliPMDv1("PMD", "PMDv1 detector"); break;
    case 2: PMD  = new AliPMDv2("PMD", "PMDv2 detector"); break;
  }  

//=================== PMD parameters ============================
}
