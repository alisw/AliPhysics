void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTART");

  AliSTART* START = 0;
  switch (version) {
    case 0: START  = new AliSTARTv0("START","START Detector"); break;
    case 1: START  = new AliSTARTv1("START","START Detector"); break;
  }  

//=================== START parameters ============================

}
