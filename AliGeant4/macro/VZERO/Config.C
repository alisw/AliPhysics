void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libVZERO");

  AliVZERO* VZERO = 0;
  switch (version) {
    case 0: VZERO  = new AliVZEROv0("VZERO", "detector VZEROv0"); break;
    case 1: break;
    case 2: VZERO  = new AliVZEROv2("VZERO", "normal VZERO"); break;
  }  

//=================== VZERO parameters ============================
}
