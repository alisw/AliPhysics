void Config(Int_t version)
{  
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libPHOS");

  AliPHOS* PHOS = 0;
  switch (version) {
    case 0: PHOS  = new AliPHOSv0("PHOS","GPS2");                break;
    case 1: PHOS  = new AliPHOSv1("PHOS","GPS2");                break;
    case 2: PHOS  = new AliPHOSv2("PHOS","GPS2");                break;
    case 3: PHOS  = new AliPHOSv3("PHOS","GPS2");                break;
    case 4: PHOS  = new AliPHOSv4("PHOS","GPS2");                break;
    case 5: PHOS  = new AliPHOSvFast("PHOS","GPS2");                break;
  }  

//=================== PHOS parameters ===========================
}  
