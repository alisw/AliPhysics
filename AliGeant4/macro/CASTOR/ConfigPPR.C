void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libCASTOR");

  AliCASTOR* CASTOR = 0;
  switch (version) {
    case 1: CASTOR  = new AliCASTORv1("CASTOR","normal CASTOR"); break;
  }  

//=================== CASTOR parameters ============================
}
