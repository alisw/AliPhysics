Int_t 
DetIndex(Int_t d, Char_t r)
{
  return (d == 1 ? 0 : (d-1) * 2 - 1) + (r == 'I' ? 0 : 1);
}
Int_t
RingColor(Int_t d, Char_t r)
{
  return (d == 1 ? kRed : d == 2 ? kGreen : kBlue) + (r == 'I' ? 1 : 3);
}

void CompareMatDep() {
  gStyle->SetOptTitle(kFALSE);

  Int_t spread = 10;
  TMultiGraph* mg = new TMultiGraph;

  for(Int_t d = 1; d<=3;d++) {
    Int_t nr = (d == 1 ? 1 : 2);
    for(Int_t q = 0; q < nr; q++) {
      Char_t        r = (q == 0 ? 'I' : 'O');
      TGraphErrors* g = new TGraphErrors();
      g->SetName(Form("fmd%d%c", d, r));
      g->SetTitle(Form("FMD%d%c (+%d)", d, r, spread*(DetIndex(d,r)+1)));
      g->SetLineColor(RingColor(d, r));
      g->SetMarkerColor(RingColor(d, r));
      g->SetFillColor(RingColor(d, r));
      g->SetFillStyle(3002);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(1.3);
      mg->Add(g);
    }
  }
  
  TGraphErrors* g = new TGraphErrors();
  g->SetName("all");
  g->SetTitle("All");
  g->SetLineColor(kBlack);
  g->SetMarkerColor(kBlack);
  g->SetFillColor(kGray);
  g->SetFillStyle(3002);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.3);
  mg->Add(g);

  const char* base = "/mnt/samsung/oldhome/ALICE/FMDanalysis/BackgroundCorrection/BackgroundStudy/TestOfMaterialBudget/";

  Int_t  samples[] = { -20, -15, -10, -5, 0, 5, 10 , 15, 20, 999 };
  Int_t* sample    = samples;
  Int_t  i         = 0;
  while (*sample < 999) { 
    TFile* fin = 0;
    if(*sample < 0) 
      fin = TFile::Open(Form("%s/minus%dpercent/fmd_dNdeta.root",
			     base,-*sample));
    else if(*sample > 0) 
      fin = TFile::Open(Form("%s/plus%dpercent/fmd_dNdeta.root",base,*sample));
    if(*sample == 0) 
      fin = TFile::Open(Form("%s/refdist/fmd_dNdeta.root", base));
    if (!fin) { 
      Error("CompareMatDep", "Failed to open input file for %d", *sample);
      sample++;
      i++;
      continue;
    }
    Float_t twa = 0;
    Float_t tsw = 0;

    for(Int_t d = 1; d<=3;d++) {
      Int_t nr = (d == 1 ? 1 : 2);
      for(Int_t q = 0; q < nr; q++) {
	Char_t        r = (q == 0 ? 'I' : 'O');
	TH1F* hRingMult = (TH1F*)fin->Get(Form("hRingMult_FMD%d%c",
					      d,r));
	Float_t wav    = 0;
	Float_t sumofw = 0;

	// Weighted average 
	for(Int_t n = 1;n<=hRingMult->GetNbinsX();n++) {
	  if(hRingMult->GetBinContent(n) > -1 && 
	     hRingMult->GetBinError(n) > 0 ) {
	    Float_t weight = 1/TMath::Power(hRingMult->GetBinError(n),2);
	    wav            += weight*hRingMult->GetBinContent(n);
	    sumofw         += weight;
	    twa            += weight*hRingMult->GetBinContent(n);
	    tsw             += weight;
	  }
	}
	
	TGraphErrors* g = 
	  static_cast<TGraphErrors*>(mg->GetListOfGraphs()->At(DetIndex(d, r)));
	g->SetPoint(i, *sample, 100*wav/sumofw+spread*(DetIndex(d,r)+1));
	g->SetPointError(i, 0, 100*1/TMath::Sqrt(sumofw)); 
	Info("CompareMatDep", "Setting point %d @ %d of %s to %f +/- %f", 
	     i, *sample, g->GetName(), 100*wav/sumofw, 
	     100*1/TMath::Sqrt(sumofw));
      } // for q 
    } // for d

    TGraphErrors* g = 
      static_cast<TGraphErrors*>(mg->GetListOfGraphs()->At(5));
    g->SetPoint(i, *sample, 100*twa/tsw);
    g->SetPointError(i, 0, 100*1/TMath::Sqrt(tsw)); 
    Info("CompareMatDep", "Setting point %d @ %d of %s to %f +/- %f", 
	 i, *sample, g->GetName(), 100*twa/tsw, 
	 100*1/TMath::Sqrt(tsw));

    sample++;
    i++;
  }
  
  TCanvas* c  = new TCanvas("c", "C", 800, 800);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.13);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetGridx();
  c->SetGridy();

  mg->Draw("a pl");
  mg->GetHistogram()->SetXTitle("#Delta#rho [%]");
  mg->GetHistogram()->SetYTitle("Average relative difference [%]");
  mg->GetHistogram()->GetXaxis()->SetTitleFont(132);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(132);
  mg->GetHistogram()->GetXaxis()->SetNdivisions(10);
  mg->GetHistogram()->GetYaxis()->SetTitleFont(132);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(132);
  mg->GetHistogram()->GetYaxis()->SetNdivisions(1010);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
  // c->Clear();
  // mg->Draw("a e3");

  TLegend* leg = gPad->BuildLegend(0.15,0.7,0.4,0.945);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->SetFillStyle(0);

  TIter next(mg->GetListOfGraphs());
  TObject* o = 0;

  c->Print("materialDependence.png");
}
