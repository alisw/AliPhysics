void dumpMapInfo(AliOADBContainer& c);
void recursiveDump(TObject* o, int depth = 0);

void dumpEtaMaps(const TString fileName)
{
  //AliLog::SetClassDebugLevel("TH2D", 0);
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  auto f = TFile::Open(fileName);

  std::vector<std::string> containerNames;
  for (const auto k : *f->GetListOfKeys()) {
    containerNames.emplace_back(k->GetName());
  }
  std::sort(containerNames.begin(), containerNames.end());

  for (const auto& name : containerNames) {
    auto c = (AliOADBContainer*)f->Get(name.data());
    printf("==================================================================\n");
    printf("Container: %s\n", name.data());
    dumpMapInfo(*c);
    printf("\n");
  }
}

void dumpMapInfo(AliOADBContainer& c)
{
  if (c.GetDefaultList()) {
    c.GetDefaultList()->Print();
  }

  for (int i = 0; i < c.GetNumberOfEntries(); ++i) {
    const int lower = c.LowerLimit(i);
    const int upper = c.UpperLimit(i);
    const char* pass = c.GetPassNameByIndex(i)->GetName();
    auto o = c.GetObjectByIndex(i);

    printf("%d - %d (%s): %s %s\n", lower, upper, pass, o->GetName(), o->GetTitle());
    recursiveDump(o, 1);
    printf("\n");
  }
}

void recursiveDump(TObject* o, int depth)
{
  const TString format = TString::Format("%%%ds %%s %%s", depth * 2);
  printf(format.Data(), "", o->GetName(), o->GetTitle());
  if (o->InheritsFrom(TCollection::Class())) {
    printf("\n");
    for (auto k : *(TCollection*)o) {
      recursiveDump(k, depth + 1);
    }
  }
  else {
    printf(" MD5 = %s\n", AliTPCPIDResponse::GetChecksum(o).Data());
  }
}
