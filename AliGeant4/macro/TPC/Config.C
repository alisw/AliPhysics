void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPC");

  AliTPC* TPC = 0;
  switch (version) {
    case 0: TPC  = new AliTPCv0("TPC","Default"); break;
    case 1: TPC  = new AliTPCv1("TPC","Default"); break;
    case 2: TPC  = new AliTPCv2("TPC","Default"); break;
    case 3: TPC  = new AliTPCv3("TPC","Default"); break;
  }   
}
