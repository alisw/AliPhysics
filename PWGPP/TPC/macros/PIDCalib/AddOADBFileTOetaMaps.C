void AddObjects(AliOADBContainer& from, AliOADBContainer& to, TString passOverWrite, Int_t firstRun, Int_t lastRun);

void AddOADBFileTOetaMaps(const TString inputFile, const TString fileToUpdate = "$ALICE_PHYSICS_SRC/OADB/COMMON/PID/data/TPCetaMaps.root", TString passOverWrite = "", Int_t firstRun = -1, Int_t lastRun = -1)
{
  TFile* fin = TFile::Open(inputFile);
  TIter nextKey(fin->GetListOfKeys());

  TFile* fupdate = TFile::Open(fileToUpdate, "update");
  TList* keysUpdata = fupdate->GetListOfKeys();

  for (auto kin : *fin->GetListOfKeys()) {
    AliOADBContainer* contin = (AliOADBContainer*)fin->Get(kin->GetName());

    if (keysUpdata->FindObject(kin->GetName())) {
      AliOADBContainer* cout = (AliOADBContainer*)fupdate->Get(kin->GetName());
      AddObjects(*contin, *cout, passOverWrite, firstRun, lastRun);
      cout->Write(0, TObject::kOverwrite);
    } else {
      contin->Write();
    }
  }

  fupdate->Close();
}

void AddObjects(AliOADBContainer& from, AliOADBContainer& to, TString passOverWrite, Int_t firstRun, Int_t lastRun)
{
  for (int i = 0; i < from.GetNumberOfEntries(); ++i) {
    int lower = from.LowerLimit(i);
    int upper = from.UpperLimit(i);
    const char* pass = passOverWrite.IsNull() ? from.GetPassNameByIndex(i)->GetName() : passOverWrite.Data();
    TObject* o = from.GetObjectByIndex(i);
    if (firstRun > 0) {
      lower = firstRun;
    }
    if (lastRun > 0) {
      upper = lastRun;
    }

    const int index = to.GetIndexForRun((lower + upper) / 2);
    if (index < 0) {
      printf("Adding object for run range %d - %d\n", lower, upper);
      to.AppendObject(o, lower, upper, pass);
    } else {
      printf("Updating object for run range %d - %d\n", lower, upper);
      to.UpdateObject(index, o, lower, upper, pass);
    }
  }
}
