//____________________________________________________________________
TObject*
GetO(TDirectory* dir, const char* name, TClass* cls=0)
{
  if (!dir) {
    Warning("GetO", "No directory passed");
    return 0;
  }

  TObject* o = dir->Get(name);
  if (!o) {
    Warning("GetO", "object %s not found in %s",
	    name, dir->GetPath());
    return 0;
  }
  if (!cls) return o;
  if (!o->IsA()->InheritsFrom(cls)) {
    Warning("GetO", "Object %s in %s is not a %s, but a %s",
	    name, dir->GetPath(), cls->GetName(), o->ClassName());
    return 0;
  }
  return o;
}
//____________________________________________________________________
TDirectory* GetD(TDirectory* dir, const char* name)
{
  return static_cast<TDirectory*>(GetO(dir,name,TDirectory::Class()));
}

//____________________________________________________________________
TH1* GetH1(TDirectory* dir, const char* name)
{
  return static_cast<TH1*>(GetO(dir,name,TH1::Class()));
}

//____________________________________________________________________
TH1* CompareOne(TDirectory* newDir,
		TDirectory* oldDir,
		Double_t    c1,
		Double_t    c2,
		Double_t&   min,
		Double_t&   max,
		TLegend*    l)
{
  TString name;
  name.Form("cent%03dd%02d_%03dd%02d",
	    Int_t(c1), Int_t(c1*100)%100,
	    Int_t(c2), Int_t(c2*100)%100);
  TDirectory* newSubDir = GetD(newDir, name);
  TDirectory* oldSubDir = GetD(oldDir, name);
  if (!newSubDir || !oldSubDir) return 0;

  TH1* newRes = GetH1(newSubDir, "dndeta");
  TH1* oldRes = GetH1(oldSubDir, "dndeta");
  if (!newRes || !oldRes) return 0;

  TH1* ratio = static_cast<TH1*>(newRes->Clone(name));
  ratio->SetDirectory(0);
  ratio->SetTitle(Form("%5.1f - %5.1f%%", c1, c2));
  ratio->SetYTitle("New / Old");
  ratio->Divide(oldRes);

  if (l) {
    TLegendEntry* e = l->AddEntry("", Form("%3.0f - %3.0f%%", c1, c2), "f");
    e->SetFillStyle(1001);
    e->SetFillColor(ratio->GetMarkerColor());
  }
  min = TMath::Min(min, ratio->GetMinimum());
  max = TMath::Max(max, ratio->GetMaximum());
  return ratio;
}

//____________________________________________________________________
void
CompareResults(const char* newName="NewTaskNewPost.root",
	       const char* oldName="OldTaskNewCorrect.root",
	       const char* newTitle="New",
	       const char* oldTitle="Old")
{
  TFile* newFile = TFile::Open(newName,"READ");
  TFile* oldFile = TFile::Open(oldName,"READ");
  if (!newFile || !oldFile) return;

  TH1* newCent = GetH1(newFile, "realCent");
  TH1* oldCent = GetH1(oldFile, "realCent");
  if (!newCent || !oldCent) return;

  TString  t; t.Form("#it{R}=#frac{%s}{%s}", oldTitle, newTitle);
  TCanvas* c     = new TCanvas("c", t, 1200, 800);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.20);
  TLegend* l     = new TLegend(1-c->GetRightMargin(),
			       c->GetBottomMargin(),
			       1, 1-c->GetTopMargin(),
			       t);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  THStack* stack = new THStack("ratios","");
			       
  Double_t min = +1e6;
  Double_t max = -1e6;
  TH1* one = 0;
  for (Int_t i = newCent->GetNbinsX(); i--;) {
    Double_t c1 = newCent->GetXaxis()->GetBinLowEdge(i+1);
    Double_t c2 = newCent->GetXaxis()->GetBinUpEdge(i+1);
    Info("", "c1=%f c2=%f", c1, c2);
    TH1*     r  = CompareOne(newFile, oldFile, c1, c2, min, max, l);    
    if (!r) continue;
    if (!one) {
      one = static_cast<TH1*>(r->Clone("one"));
      one->SetDirectory(0);
      one->Reset();
      for (Int_t j = 1; j <= one->GetNbinsX(); j++) {
	one->SetBinContent(j,1);
	one->SetBinError  (j,0);
      }
    }
    // r->Add(one, i-1);
    // r->Scale(TMath::Power(10,i));
    stack->Add(r);
  }
  stack->Draw("nostack");
  stack->SetMinimum(0.95*min);
  stack->SetMaximum(1.05*max);
  stack->GetHistogram()->SetXTitle("#eta");
  stack->GetHistogram()->SetYTitle("#it{R}");
  l->Draw();
  c->Modified();
  c->Update();
  c->cd();
  c->SaveAs(Form("%sover%s.png", oldTitle, newTitle));
}

  
  
