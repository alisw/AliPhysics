void
Example(Bool_t fit=true)
{
  // Load the class - if not already done 
  if (!gROOT->GetClass("GraphSysErr"))
    gROOT->LoadMacro("GraphSysErr.C+g");

  // Adjust size along X of common errors 
  gStyle->SetErrorX(.2);
  // Adjust size of hat, cap, ... - depends on canvas size!
  gStyle->SetEndErrorSize(10);

  // Make our object 
  GraphSysErr* gse = new GraphSysErr("foo", "Gaussian");
  // Draw data with-out ticks 
  gse->SetDataOption(GraphSysErr::kNoTick);
  gse->SetXTitle("X");
  gse->SetYTitle("Y");
  // Set some key/value pairs 
  gse->SetKey("laboratory", "The Center");
  gse->SetKey("accelerator", "Mega Collider");
  gse->SetKey("detector", "Huge Experiment");
  gse->SetKey("author", "Christensen");
  gse->SetKey("reference","Jour.All.Things A1,999");
  gse->SetKey("doi","9999-9999-9999-9999");
  gse->SetKey("abstract", "The data");
  gse->SetKey("location", "In the paper");
  gse->SetKey("reackey", "graviton -> tachyons");
  gse->SetKey("obskey", "GUT");
  // Adding qualifiers 
  gse->AddQualifier("question", "Life, universe, and everything");
  
  // Two sources of common errors one relative, one absolue
  UInt_t cm1 = gse->DefineCommon("Common 0.05", false, .05);
  UInt_t cm2 = gse->DefineCommon("Common 10%", true, .1);
  
  // Two sources of point-to-point errors, one relative, one absolute
  UInt_t pp1 = gse->DeclarePoint2Point("Point-to-Point 0.1-0.2", true);
  UInt_t pp2 = gse->DeclarePoint2Point("Point-to-Point 5-10%", false);
  
  // Set options on summed errors (in case of option COMBINED)
  gse->SetSumLineColor(kRed+2);
  gse->SetSumLineWidth(2);
  gse->SetSumTitle("All errors");
  gse->SetSumOption(GraphSysErr::kHat);
  
  // Set attributes of common errors 
  gse->SetSysFillColor(cm1, kRed+2);
  gse->SetSysFillStyle(cm1, 3001);
  gse->SetSysLineColor(cm1, kRed+2);
  gse->SetSysFillColor(cm2, kCyan+2);
  gse->SetSysFillStyle(cm2, 3001);
  gse->SetSysOption(cm1, GraphSysErr::kBox);
  gse->SetSysOption(cm2, GraphSysErr::kRect);
  
  // Set attributes of other errors 
  gse->SetSysLineColor(pp1, kBlue+2);
  gse->SetSysLineWidth(pp1, 2);
  gse->SetSysLineColor(pp2, kGreen+2);
  gse->SetSysLineWidth(pp2, 3);
  gse->SetSysOption(pp1, GraphSysErr::kBar);
  gse->SetSysOption(pp2, GraphSysErr::kHat);
  
  // Fill a histogram with a Guassian random deviate 
  TH1* h = new TH1F("h", "h", 30, -3, 3);
  h->Sumw2();
  h->SetDirectory(0);
  h->FillRandom("gaus",1000);
  h->Scale(1./1000, "width");
  
  // Fill in the data points 
  for (Int_t i = 0; i < h->GetNbinsX(); i++) { 
    Int_t    bin = i+1;
    Double_t x   = h->GetXaxis()->GetBinCenter(bin);
    Double_t y   = h->GetBinContent(bin);
    Double_t sta = h->GetBinError(bin);
    Double_t w   = h->GetXaxis()->GetBinWidth(bin);
    
    // Set data 
    gse->SetPoint(i, x, y);
    gse->SetPointError(i, w/2, w/2);
    gse->SetStatError(i, sta);
    
    // Set point-to-point errors 
    gse->SetSysError(pp1, i, 0., gRandom->Uniform(0.1, 0.2));
    gse->SetSysError(pp2, i, 0., 0., 
		     gRandom->Uniform(0.05, 0.1),
		     gRandom->Uniform(0.05, 0.1));
  } 
  // Remove temporary histogram
  delete h;

  // Build our canvas 
  TCanvas* c = new TCanvas("c","c", 1400, 1000);
  c->SetFillColor(0);
  c->SetFillStyle(0);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);
  
  // Draw or fit (and draw) a Guassian to the data
  const char* option = "STACK stat axis quad split max west";
  if (!fit) 
    gse->Draw(option);
  else 
    gse->Fit("gaus", "SQ", option, -3, 3);

  // Make a legend 
  TLegend* l = c->BuildLegend(0.7,0.7,0.97,0.97);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  // update the canvas and print
  c->Modified();
  c->Update();
  c->cd();
  c->Print("Example.png");
}

  
