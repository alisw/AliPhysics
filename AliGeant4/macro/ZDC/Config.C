void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libZDC");

  AliZDC* ZDC = 0;
  switch (version) {
    case 1: ZDC  = new AliZDCv1("ZDC","normal ZDC"); break;
  }  

//=================== ZDC parameters ============================
}  
