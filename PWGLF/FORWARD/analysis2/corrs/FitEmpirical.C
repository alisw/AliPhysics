void FitEmpirical(Bool_t update=false)
{
  const char* path="$ALICE_PHYSICS_SOURCE/OADB/PWGLF/FORWARD/EMPIRICAL";
  Int_t       runNo=138190;
  const char* which="default";
  
  TFile* file = TFile::Open(Form("%s/empirical_%09d.root",path,runNo),
			    (update ? "UPDATE" : "READ"));
  if (!file) return 0;
  Printf("File opened: %s", file->GetPath());

  TObject* obj = file->Get(Form("Forward/%s", which));
  if (!obj) return 0;

  if (!obj->IsA()->InheritsFrom(TH1::Class())) return 0;

  TH1* h = static_cast<TH1*>(obj);
  h->SetDirectory(0);
  h->Draw();

  TF1* fc0 = new TF1("fc0", "pol3", -3.5, -1.7);
  fc0->SetLineColor(kRed+1);
  fc0->SetParNames("C10","C11","C12","C13");
  h->Fit(fc0, "+R");

  TF1* fa0 = new TF1("fa0", "pol3", 1.7, 2.075);
  fa0->SetLineColor(kGreen+3);
  fa0->SetParNames("A00","A01","A02","A03","A04");
  h->Fit(fa0, "+R");

  TF1* fa1 = new TF1("fa1", "pol2", fa0->GetXmax(), 2.3);
  fa1->SetLineColor(kGreen+2);
  fa1->SetParNames("A10","A11","A12","A13","A14");
  h->Fit(fa1, "+R");

  TF1* fa2 = new TF1("fa2", "pol4", fa1->GetXmax(), 3.7);
  fa2->SetLineColor(kGreen+1);
  fa2->SetParNames("A20","A21","A22","A23","A24");
  h->Fit(fa2, "+R");

  TF1* fa3 = new TF1("fa3", "pol4", fa2->GetXmax(), 5);
  fa3->SetLineColor(kBlue+1);
  fa3->SetParNames("A30","A31","A32","A33","A34");
  h->Fit(fa3, "+R");

  TArrayD   pars(30);
  TObjArray names(30);
  TString   formula;
  Int_t     j    = 0;
  TF1*      fs[] = { fc0, fa0, fa1, fa2, fa3, 0 };
  TF1**     fp   = fs;
  while (*fp) {
    TF1* ft = *fp;
    fp++;
    // if (ft == fa1) continue;
    Printf("Function %p: %s [%f,%f] %d", ft, ft->GetTitle(),
	   ft->GetXmin(), ft->GetXmax(), ft->GetNpar());
    Double_t    x1   = ft->GetXmin(); // if (x1 < -3) x1 = -3.5;
    Double_t    x2   = ft->GetXmax();
    Bool_t      incl = (x2 < 0) || x2 > 4.5;
    // if (ft == fa0) { x2 = 2.3; fa0->SetRange(x1,x2); }
    // if (ft == fa2) { x1 = 2.1; fa2->SetRange(x1,x2); }
    const char* expr = ft->GetTitle();
    if (!formula.IsNull()) formula.Append("+");
    formula.Append(Form("(%5.2f<=x&&x<%s%5.2f)*%s(%2d)",
			x1, (incl ? "=" : ""), x2, expr, j));
    for (Int_t i = 0; i < ft->GetNpar(); i++, j++) {
      pars[j] = ft->GetParameter(i);
      names.AddAt(new TObjString(ft->GetParName(i)), j);
    }
  }

  
  TF1* f = new TF1("param", formula.Data(), -4, 6);
  f->SetLineColor(kGray+1);
  for (Int_t i = 0; i < j; i++) {
    f->SetParName(i, names.At(i)->GetName());
    f->SetParameter(i, pars[i]);
  }
  // h->Fit(f,"+ME");
  f->Draw("Same");
  f->SetNpx(1000);
	 
  file->cd("Forward");
  gDirectory->ls();
  if (update) f->Write();
  gDirectory->ls();
  
		   
}
