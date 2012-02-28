void UpdateEventplaneOADB(TString oadbname, const char* updatename, Int_t runBegin = 0, Int_t runEnd = 0)
{

  Bool_t newcont = kFALSE;
  
  if (!runBegin || !runEnd) newcont = kTRUE;
  
  TFile* oadb = 0; 
  if (oadbname.Length() > 0.) oadb = TFile::Open(oadbname, "READ");
  TFile* in   = TFile::Open(updatename);

  AliOADBContainer* cont = 0;
  if (oadb) cont = (AliOADBContainer*)oadb->Get("FMDphidist");
  if (!cont) {
    if (newcont) cont = new AliOADBContainer("FMDphidist");
    else         Fatal("Something is wrong. There is no container, but you supplied a runrange...\n");
  }

  TList* list   = (TList*)in->Get("FMDEventplaneSums");
  TList* eplist = (TList*)list->FindObject("fmdEventplaneFinder");
  TH1D*  hist   = (TH1D*) eplist->FindObject("hPhiDist");

  if (!newcont) {
    hist->SetName(Form("%d-%d", runBegin, runEnd));
    cont->AppendObject(hist, runBegin, runEnd);
  }
  else {
    hist->SetName("Default");
    cont->AddDefaultObject(hist);
  }

  TFile* out = TFile::Open("new_fmdEPoadb.root", "RECREATE");
  out->Close();

  cont->WriteToFile("new_fmdEPoadb.root");

  Printf("Wrote new OADB object to file new_%s, please check that everything is OK and replace the old file", oadbname);

  out = TFile::Open("new_fmdEPoadb.root");
  new TBrowser();

}
