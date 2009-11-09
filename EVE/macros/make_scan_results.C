// Label is used to store scanning result for tracks and trackelts.
//
// BIT(1) stores the original selection.
// BIT(0) stores the user selection (set to same value as b1 at init).
//
// This allows to check all possible combinations.

struct XXX
{
  const char   *bname;
  const char   *oname;
  TBranch      *branch;
  AliESDVertex *vert;
};

XXX vvv[3] = {
  { "VT",   "PrimVertTracks" },
  { "VTPC", "PrimVertTPC" },
  { "VSPD", "PrimVertSPD" }
};

void make_scan_results()
{
  TFile *f = TFile::Open("scan_results.root", "UPDATE");

  f->Delete("SR;*");

  T = new TTree("SR", "Scanning results");

  TClonesArray* ts = new TClonesArray("AliESDtrack", 32);
  TBranch * tb = T->Branch("T", &ts);
  delete ts;

  AliMultiplicity *ms = 0;
  TBranch *mb = T->Branch("M", &ms);

  for (Int_t v = 0; v < 3; ++v)
  {
    vvv[v].vert   = 0;
    vvv[v].branch = T->Branch(vvv[v].bname, &vvv[v].vert);
  }

  for (Int_t i=0; i<=9999; ++i)
  {
    TString name;

    name.Form("Tracks_%04d", i);
    ts = (TClonesArray*) f->Get(name);
    if (ts == 0)
      continue;

    name.Form("Tracklets_%04d", i);
    ms = (AliMultiplicity*) f->Get(name);
    if (ms == 0)
      Error("make_scan_results", "'%s' not found.", name.Data());

    tb->SetAddress(&ts);
    mb->SetAddress(&ms);

    for (Int_t v = 0; v < 3; ++v)
    {
      name.Form("%s_%04d", vvv[v].oname, i);
      vvv[v].vert = (AliESDVertex*) f->Get(name);
      if (vvv[v].vert == 0)
        Error("make_scan_results", "'%s' not found.", name.Data());
      vvv[v].branch->SetAddress(&vvv[v].vert);
    }

    T->Fill();

    delete ts;
    delete ms;
    for (Int_t v = 0; v < 3; ++v) delete vvv[v].vert;
  }

  T->Write();

  f->Close();
  delete f;
}
