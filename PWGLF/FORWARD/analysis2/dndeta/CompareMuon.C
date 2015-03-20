//====================================================================
TCanvas* MakeCanvas(const char* name, const char* title,
		    const char* img=0)
{
  TImage* dimg = 0;
  if (img) {
    TImage* dimg = TImage::Open(img);
    if (!dimg) {
      Error("Failed to open image file: %s", img);
      return 0;
    }
  }
  UInt_t w = (dimg ? dimg->GetWidth() : 1000);
  UInt_t h = (dimg ? dimg->GetHeight() : 900);
  Printf("Canvas size: (%dx%d)", w, h);
  TCanvas* c = new TCanvas(name, title, w, h);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);

  if (dimg) {
    c->SetTopMargin(0);
    c->SetRightMargin(0);
    c->SetLeftMargin(0);
    c->SetRightMargin(0);
    c->cd();
    // dimg->Scale(c->GetWw(), c->GetWh());
    dimg->Draw("x");
    // dimg->SetDrawOption("X");

    TPad* p = new TPad("overlay", "Overlay", 0, 0, 1, 1, 0, 0, 0);
    p->SetFillStyle(0);
    p->SetFillColor(0);
    p->SetLineColor(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetNumber(10);
    c->cd();
    p->Draw();
  }

  c->Modified();
  c->Update();
  c->cd();

  return c;
}

TPair*
GetMine(const TString& trigger,
	Bool_t rebinned=true,
	Bool_t empirical=true)
{
  gROOT->SetMacroPath(Form("%s:%s", gROOT->GetMacroPath(),
			   "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta"));

  if (!gROOT->GetClass("Drawer")) gROOT->LoadMacro("Drawer.C++");

  TString  system("pPb");
  TString  system2("Pbp");
  UShort_t sNN  = 5023;
  
  TString t1(trigger); t1.ReplaceAll("X", "A");
  TString t2(trigger); t2.ReplaceAll("X", "C");
  THStack* s1 = Drawer::GetStack(0, system,  sNN, t1, rebinned, empirical);
  THStack* s2 = Drawer::GetStack(0, system2, sNN, t2, rebinned, empirical);

  if (!s1 || !s2) return 0; 

  THStack* res = Drawer::Symmetrice(s1, s2);
  TPair*   ret = new TPair(res, 0);

  TMultiGraph* other = Drawer::GetOther(system, sNN, t1);
  if (!other) 
    Warning("", "No other data for %s,%d,%s", 
	    system.Data(),sNN,trigger.Data());
  else
    ret->SetValue(other);

  return ret;
}

void CompareMuonDensity(const char* trig)
{
  gStyle->SetFrameLineStyle(0);

  const char* png
    = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta/dNdeta_muon.png";
  TCanvas* c = MakeCanvas("CompareMuonDensity", "Muon",
			  gSystem->ExpandPathName(png));

  TPair* mine = GetMine(trig);
  if (!mine) {
    Warning("", "Didn't get mine");
    return;
  }
  mine->Print();
	    
  TVirtualPad* p = c->GetPad(10);
  if (!p) {
    Error("CompareMuon", "No overlay pad");
    return;
  }
  p->SetTopMargin(0.08);
  p->SetRightMargin(0.06);
  p->SetBottomMargin(0.16);
  p->SetLeftMargin(0.18);
  // p->SetFillStyle(1001);
  // p->SetFillColor(kRed);
  p->Modified();
  p->Update();
  
  p->cd();
  THStack* data = static_cast<THStack*>(mine->Key());
  data->SetTitle(Form("1/N dN_{ch}/d#eta (%s) and 1/#deltac dN_{#mu}/d#eta",
		      trig));
  data->Draw("nostack y+");
  data->GetHistogram()->GetXaxis()->SetRangeUser(-4.5,4.5);
  data->GetHistogram()->GetXaxis()->SetLabelSize(0);  
  data->GetHistogram()->GetYaxis()->SetNdivisions(305);
  data->GetHistogram()->GetXaxis()->SetNdivisions(305);
  // data->SetMaximum(30);
  if (mine->Value()) {
    mine->Value()->Draw("p");
  }
  TFrame* f =
    static_cast<TFrame*>(p->GetListOfPrimitives()->FindObject("TFrame"));
  if (!f) {
    Warning("", "No frame");
  }
  else {
    f->SetFillStyle(0);
  }

  c->Modified();
  c->Update();
  c->cd();

  c->Print(Form("plots/CompareMuonDensity%s.png",trig));
}


void CompareMuonRatio(const char* trig)
{
  gStyle->SetFrameLineStyle(0);
  
  const char* png
    = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta/dNdeta_muon_ratio.png";
  TCanvas* c = MakeCanvas("CompareMuonRatio", "Muon",
			  gSystem->ExpandPathName(png));

  TPair* mine = GetMine(trig);
  if (!mine) {
    Warning("", "Didn't get mine");
    return;
  }
  mine->Print();
	    
  TVirtualPad* p = c->GetPad(10);
  if (!p) {
    Error("CompareMuon", "No overlay pad");
    return;
  }
  p->SetTopMargin(0.08);
  p->SetRightMargin(0.06);
  p->SetBottomMargin(0.16);
  p->SetLeftMargin(0.18);
  // p->SetFillStyle(1001);
  // p->SetFillColor(kRed);
  p->Modified();
  p->Update();
  
  p->cd();

  THStack* data   = static_cast<THStack*>(mine->Key());
  THStack* ratios = new THStack(*data);
  TList*   dlist  = ratios->GetHists();
  TH1*     dperif = static_cast<TH1*>(dlist->Last());
  dlist->Remove(dperif);
  TIter  dnext(dlist);
  TH1*   dh = 0;
  while ((dh = static_cast<TH1*>(dnext()))) {
    dh->Divide(dperif);
    dh->SetMaximum(37.3);
  }
  ratios->SetTitle(Form("R_{cp} for 1/N dN_{ch}/d#eta (%s) "
			"and 1/#deltac dN_{#mu}/d#eta", trig));
  ratios->Draw("nostack y+");
  ratios->GetHistogram()->GetXaxis()->SetRangeUser(-4.5,4.5);
  ratios->GetHistogram()->GetYaxis()->SetNdivisions(305);
  ratios->GetHistogram()->GetXaxis()->SetLabelSize(0);  
  ratios->GetHistogram()->GetXaxis()->SetNdivisions(305);
  // ratios->GetHistogram()->SetMaximum(39);
  
  if (mine->Value()) {
    TMultiGraph*       odata   = static_cast<TMultiGraph*>(mine->Value());
    TMultiGraph*       oratios = new TMultiGraph;
    TList*             olist   = odata->GetListOfGraphs();
    olist->Sort();
    TGraphAsymmErrors* operif  = static_cast<TGraphAsymmErrors*>(olist->Last());
    olist->Remove(operif);

    TIter              onext(olist);
    TGraphAsymmErrors* og = 0;
    while ((og = static_cast<TGraphAsymmErrors*>(onext())))
      oratios->Add(Drawer::GOverG(og, operif, 0, 0), "p");

    oratios->Draw("p");
  }
  TFrame* f =
    static_cast<TFrame*>(p->GetListOfPrimitives()->FindObject("TFrame"));
  if (!f) {
    Warning("", "No frame");
  }
  else {
    f->SetFillStyle(0);
  }

  c->Modified();
  c->Update();
  c->cd();

  c->Print(Form("plots/CompareMuonRatio%s.png",trig));
}
  
void
CompareMuon(const char* trig="CENTZNX")
{
  CompareMuonDensity(trig);
  CompareMuonRatio(trig);
}
