void Config(Int_t version)
{
  gSystem->Load("libmicrocern");
  gSystem->Load("libMUON");

  AliMUON *MUON = 0;
  switch (version) {
    case 0: MUON  = new AliMUONv0("MUON","normal MUON"); break;
    case 1: MUON  = new AliMUONv1("MUON","default"); break;
  }  
}
