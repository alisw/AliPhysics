void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libCRT");

  AliCRT* CRT = 0;
  switch (version) {
    case 0: CRT  = new AliCRTv0("CRT", "normal ACORDE"); break;
  }  

//=================== EMCAL parameters ============================
}
