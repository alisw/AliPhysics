void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libRICH");

  AliRICH *RICH = 0;
  switch (version) {
    case 0: RICH  = new AliRICHv0("RICH","normal RICH"); break;
    case 1: RICH  = new AliRICHv1("RICH","normal RICH"); break;
    case 2: RICH  = new AliRICHv2("RICH","normal RICH"); break;
  }  
}

