void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libEMCAL");

  AliEMCAL* EMCAL = 0;
  switch (version) {
    case 0: EMCAL  = new AliEMCALv0("EMCAL", "EMCALv0"); break;
    case 1: EMCAL  = new AliEMCALv1("EMCAL", "EMCALArch1a"); break;
  }  

//=================== EMCAL parameters ============================
}
