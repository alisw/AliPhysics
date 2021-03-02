void AddObjects(AliOADBContainer& from, AliOADBContainer& to, TString passOverWrite);

void AddOADBFileTOetaMaps(const TString inputFile, const TString fileToUpdate = "$ALICE_PHYSICS_SRC/OADB/COMMON/PID/data/TPCetaMaps.root", TString passOverWrite = "")
{
  TFile* fin = TFile::Open(inputFile);
  TIter nextKey(fin->GetListOfKeys());

  TFile* fupdate = TFile::Open(fileToUpdate, "update");
  TList* keysUpdata = fupdate->GetListOfKeys();

  for ( auto kin : *fin->GetListOfKeys() ) {
    AliOADBContainer* contin = (AliOADBContainer*)fin->Get(kin->GetName());

    if ( keysUpdata->FindObject(kin->GetName()) ) {
      AliOADBContainer* cout = (AliOADBContainer*)fupdate->Get(kin->GetName());
      AddObjects(*contin, *cout, passOverWrite);
      cout->Write(0, TObject::kOverwrite);
    }
    else {
      contin->Write();
    }
  }

  fupdate->Close();
}

void AddObjects(AliOADBContainer& from, AliOADBContainer& to, TString passOverWrite)
{
  for (int i = 0; i < from.GetNumberOfEntries(); ++i) {
    const int lower = from.LowerLimit(i);
    const int upper = from.UpperLimit(i);
    const char* pass = passOverWrite.IsNull() ? from.GetPassNameByIndex(i)->GetName() : passOverWrite.Data();
    TObject* o = from.GetObjectByIndex(i);

    to.AppendObject(o, lower, upper, pass);
  }
}
