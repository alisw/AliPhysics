void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTRUCT");

  AliBODY* BODY = 0;
  switch (version) {
    case 0: BODY = new AliBODY("BODY","Alice envelop"); break;
  }  

//=================== Alice BODY parameters =============================
}  
