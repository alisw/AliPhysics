void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libminicern");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libMUON");

  AliMUON *MUON = 0;
  switch (version) {
    case 0: MUON  = new AliMUONv0("MUON","normal MUON"); break;
    case 1: MUON  = new AliMUONv1("MUON","normal MUON"); break;
  }  
}
