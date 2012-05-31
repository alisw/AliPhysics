{

  TFile *f = TFile::Open("$ALICE_ROOT/OCDB/EMCAL/Align/Data/Run0_999999999_v0_s0.root");
  AliCDBEntry *entry = (AliCDBEntry*)f->Get("AliCDBEntry");
  TClonesArray *aligndata = dynamic_cast<TClonesArray*>  (entry->GetObject());
  aligndata->Print("");
}
