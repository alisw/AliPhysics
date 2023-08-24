void dumpPIDResponse(const TString name)
{
  AliOADBContainer c("TPCSplines");
  c.InitFromFile(name, "TPCSplines");

  AliTPCPIDResponse tpcResp;

  for (int i = 0; i < c.GetNumberOfEntries(); ++i) {
    const int lower = c.LowerLimit(i);
    const int upper = c.UpperLimit(i);
    const char* pass = c.GetPassNameByIndex(i)->GetName();
    const TObjArray* arr = (TObjArray*)c.GetObjectByIndex(i);

    printf("==================================================================\n");
    printf("%d - %d (%s): %s\n", lower, upper, pass, arr->GetName());
    arr->Print();
    tpcResp.SetSplinesFromArray((TObjArray*)arr->FindObject("Splines"));
  }
}
