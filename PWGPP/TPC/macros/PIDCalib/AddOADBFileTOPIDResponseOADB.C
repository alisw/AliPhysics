void AddObjects(AliOADBContainer& from, AliOADBContainer& to);

void AddOADBFileTOPIDResponseOADB(const TString inputOADB, const TString pidResponseOADBin = "$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCPIDResponseOADB.root", const TString pidResponseOADBout = "$ALICE_PHYSICS_SRC/OADB/COMMON/PID/data/TPCPIDResponseOADB.root")
{
  AliOADBContainer cin("TPCSplines");
  cin.InitFromFile(inputOADB, "TPCSplines");

  AliOADBContainer cout("TPCSplines");
  cout.InitFromFile(pidResponseOADBin, "TPCSplines");

  AddObjects(cin, cout);

  cout.WriteToFile(pidResponseOADBout);
}

void AddObjects(AliOADBContainer& from, AliOADBContainer& to)
{
  for (int i = 0; i < from.GetNumberOfEntries(); ++i) {
    const int lower = from.LowerLimit(i);
    const int upper = from.UpperLimit(i);
    const char* pass = from.GetPassNameByIndex(i)->GetName();
    TObject* o = cin.GetObjectByIndex(i);

    to.AppendObject(o, lower, upper, pass);
  }
}
