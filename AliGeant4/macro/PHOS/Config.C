void Config(Int_t version)
{  
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libPHOS");

  AliPHOS* PHOS = 0;
  switch (version) {
    case 0: PHOS  = new AliPHOSv0("PHOS","GPS2");                break;
    case 1: PHOS  = new AliPHOSv1("PHOS","GPS2");                break;
  }  

//=================== PHOS parameters ===========================
}  
