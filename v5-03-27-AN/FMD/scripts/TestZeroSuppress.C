TGraph*
makeInput(Int_t n, TArrayI& vals, TF1* f)
{
  vals.Set(n);
  TGraph* data = new TGraph(n);
  data->SetName("input");
  data->SetTitle("Input");
  data->SetMarkerColor(kRed);
  data->SetLineColor(kRed);
  data->SetMarkerStyle(20);
  
  for (Int_t i = 0; i < data->GetN(); i++) { 
    Double_t v = f->Eval(i);
    data->SetPoint(i, i, v);
    vals[i]   = Int_t(v);
  }
  return data;
}

TGraph* makeFlat(const char* name, Int_t n, Float_t val, TArrayF& vals)
{
  vals.Set(n);
  vals.Reset(val);
  TGraph* g = new TGraph(n);
  g->SetName(name);
  g->SetTitle(name);
  g->SetMarkerStyle(21);
  
  for (Int_t i = 0; i < g->GetN(); i++) 
    g->SetPoint(i, i, val);
  return g;
}
  
void
runTest(TF1* f, AliFMDRawWriter* w, Float_t pv, Float_t nv)
{
  static Int_t num = 1;
  
  Int_t    n         = Int_t(f->GetXmax());
  UShort_t threshold = 1;
  TArrayI  vals;
  TArrayF  peds;
  TArrayF  noise;
  TArrayF  dummy;
  TGraph*  input      = makeInput(n, vals, f);
  TGraph*  gPeds      = makeFlat("peds",      n, pv,        peds);
  TGraph*  gNoise     = makeFlat("noise",     n, nv,        noise);
  TGraph*  gThreshold = makeFlat("threshold", n, threshold, dummy);
  gThreshold->SetLineStyle(2);
  gPeds->SetLineColor(kGreen);
  gNoise->SetLineColor(kMagenta);

  w->ZeroSuppress(vals.fArray, n, peds.fArray, noise.fArray, threshold);

  TGraph* output = new TGraph(n);
  output->SetName("output");
  output->SetLineColor(kBlue);
  output->SetMarkerColor(kBlue);
  output->SetMarkerStyle(21);
  for (Int_t i = 0; i < output->GetN(); i++) 
    output->SetPoint(i, i, vals[i]);

  TCanvas* c = new TCanvas(Form("c%02d", num), 
			   Form("Zero suppression test %d", num));
  c->SetFillColor(kWhite);
  input->Draw("APL");
  input->GetHistogram()->SetMinimum(0);
  gPeds->Draw("L same");
  gThreshold->Draw("L same");
  output->Draw("PL same");
  f->SetLineStyle(3);
  f->Draw("same");
  num++;
}

  

void
TestZeroSuppress()
{
  AliLog::SetModuleDebugLevel("FMD", 1);
  AliFMDRawWriter* w = new AliFMDRawWriter(0);

  TF1* simple = new TF1("simple", "[0] + [1] * (x >= [2] && x <= [3])", 0, 16);
  simple->SetParameters(4, 10, 2.5, 6.5);

  runTest(simple, w, 3.5, 0);

  w.SetPedSubtract();
  runTest(simple, w, 3.5, 0);

  w.SetPedSubtract(kFALSE);
  simple->SetParameters(4, 10, 3.5, 4.5);
  runTest(simple, w, 3.5, 0);

  TF1* two = new TF1("two", 
		     "[0]+[1]*((x>=[2]&&x<=[3])||(x>=[4]&&x<=[5]))", 0, 16);
  two->SetParameters(4, 10, 2.5, 4.5, 9.5, 12.5);
  runTest(two, w, 3.5, 0);

}

  
