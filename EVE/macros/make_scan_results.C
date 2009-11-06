// Label is used to store scanning result for tracks and trackelts.
//
// BIT(1) stores the original selection.
// BIT(0) stores the user selection (set to same value as b1 at init).
//
// This allows to check all possible combinations.


void make_scan_results()
{
  TFile *f = TFile::Open("scan_results.root", "UPDATE");

  f->Delete("SR;*");

  T = new TTree("SR", "Scanning results");

  TClonesArray* clones = new TClonesArray("AliESDtrack", 32);
  TBranch * tb = T->Branch("T", &clones);

  AliMultiplicity *mult = 0;
  TBranch *mb = T->Branch("M", &mult);


  for (Int_t i=0; i<=9999; ++i)
  {
    TString name;

    name.Form("Tracks_%04d", i);
    TClonesArray* ts = (TClonesArray*) f->Get(name);

    name.Form("Tracklets_%04d", i);
    AliMultiplicity* ms =  (AliMultiplicity*) f->Get(name);

    if (ts && ms)
    {
      tb->SetAddress(&ts);
      mb->SetAddress(&ms);
      T->Fill();
    }
    else if ((ts && !ms) || (!ts && ms))
    {
      Error("make_scan_results", "Only one of tracks/tracklets exists for index %d.", i);
    }

  }

  T->Write();

  f->Close();
  delete f;
}
