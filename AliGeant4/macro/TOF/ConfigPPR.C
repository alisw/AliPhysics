void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTOF");

  AliTOF *TOF = 0;
  switch (version) {
    case 0: TOF  = new AliTOFv0("TOF", "TOFv0 detector"); break;
    case 1: TOF  = new AliTOFv1("TOF","normal TOF");      break;
    case 2: TOF  = new AliTOFv2("TOF", "TOFv2 detector"); break;
    case 3: TOF  = new AliTOFv3("TOF", "TOFv3 detector"); break;
    case 4: TOF  = new AliTOFv4("TOF", "TOFv4 detector"); break;
  }   

//=================== TOF parameters ============================
}  
