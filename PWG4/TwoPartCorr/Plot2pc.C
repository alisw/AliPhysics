
// Simple macro to plot AutoCorr output histos

void Plot2pc(const Int_t trial = 6, Bool_t printPDF = false)
{
  const char* histFile = Form("output/try%i/try%i.root", trial, trial);

  TObjArray* c1List = new TObjArray();
  TObjArray* c2List = new TObjArray();
  TObjArray* TH1List = new TObjArray();
  TObjArray* TH2List = new TObjArray();

  TFile* inFile = new TFile(histFile, "read");

  // TODO: Move this to another macro or a histo manager
  int nMultBinsRead = hMultBins->GetNbinsX();
  const int nMultBins = nMultBinsRead;
  int nPtBinsRead = hPtBins->GetNbinsX();
  const int nPtBins = nPtBinsRead;

  TString sMultBin[nMultBins];
  for (int i=0; i<nMultBins; i++) {
    double m1 = hMultBins->GetBinLowEdge(i+1);
    double m2 = m1 + hMultBins->GetBinWidth(i+1);
    sMultBin[i] = Form("%.0f - %.0f tracks", m1, m2);
  }
  TString sPtBin[nPtBins];
  for (int i=0; i<nPtBins; i++) {
    double pt1 = hPtBins->GetBinLowEdge(i+1);
    double pt2 = pt1 + hPtBins->GetBinWidth(i+1);
    sPtBin[i] = Form("%.2f - %.2f GeV/c", pt1, pt2);
  }

  gStyle->SetOptStat(0); // Turn on indiv. stats w/ TH1::SetStats(1)


  TIter next(inFile->GetListOfKeys()); 
  TKey *key;
  while ((key=(TKey*)next())) {
    TString className(key->GetClassName());
    TString keyName(key->GetName());
    printf("%10s %20s\n", className.Data(), keyName.Data());

    // Get TH1Fs
    if (className.Contains("TH1F")) {
      TH1F* h1 = (TH1F*)inFile->Get(keyName.Data());
      TH1List->Add(h1);
      Int_t i = TH1List->IndexOf(h1);
      c1List->Add(new TCanvas(Form("c1_%i",i), 
			      Form("c1_%s",keyName.Data()), 1));
    }

    // Get TH2Fs
    if (className.Contains("TH2F")) {
      TH2F* h2 = (TH2F*)inFile->Get(keyName.Data());
      TH2List->Add(h2);
      Int_t i = TH2List->IndexOf(h2);
      c2List->Add(new TCanvas(Form("c2_%i",i), 
			      Form("c2_%s",keyName.Data()), 1));
    }


  }

  // Give reasonable titles
  for (int i=0; i<TH2List->GetEntries(); i++) {
    for (int iM=0; iM<nMultBins; iM++) {
      for (int iPt=0; iPt<nPtBins; iPt++) {
	if (keyName.EndsWith(Form("_%i_%i", iM, iPt))) {
	  TH2F* h = (TH2F*)TH2List->At(i);
	  h->SetTitle(Form("%s, %s;#Delta#phi;#Delta#eta;", 
			   sMultBin[iM].Data(), sPtBin[iPt].Data()));
	}
      }
    }
  }
  
  cuts->Print();

  Int_t nPrinted = 0; // Number of canvases printed

  // Draw TH1s  
  for (int i=0; i<c1List->GetEntries(); i++) {
    TCanvas* c = (TCanvas*)c1List->At(i); 
    TH1F* h1 = (TH1F*)TH1List->At(i);
    TString s(h1->GetName());
    TString options("");
    Bool_t skipPrint = false;

    if (s.Contains("Bins")) { 
      skipPrint = true;
      continue;
    }
    if (s.Contains("hPt") || s.Contains("hMass")) 
      c->SetLogy();
    if (s.Contains("Sel") || s.Contains("hMassBg")) {
      h1->SetLineColor(kBlue);
      c = (TCanvas*)c1List->At(i-1); 
      options = "same";
      skipPrint = true;
    }
    if (s.Contains("hMult"))
      h1->SetStats(1);


    c->cd();
    h1->Draw(options.Data());


    if (printPDF) {
      c->Print(Form("plots/pngs/%s_try%i.png", c->GetName(), trial));                          
      if (nPrinted==0)
	// "[" opens ps file but doesn't print c
	c->Print(Form("plots/try%i.ps%s", trial, "["));
      
      c->Print(Form("plots/try%i.ps%s", trial, ""));
      nPrinted++;
    }
  }

  // Draw TH2s
  for (int i=0; i<c2List->GetEntries(); i++) {
    TCanvas* c = (TCanvas*)c2List->At(i);
    TH2F* h2 = (TH2F*)TH2List->At(i);
    TString s(h2->GetName());
    TString options("");
    Bool_t skipPrint = false;

    if (s.Contains("Sig") || 
	s.Contains("Bkg") || 
	s.Contains("SB")  ||
	s.Contains("hDR")) 
      options = "surf1";
    if (s=="hMixEff")
      options = "colz";


    c->cd();
    h2->Draw(options.Data());

    if (printPDF) {
      c->Print(Form("plots/pngs/%s_try%i.png",         
		       c->GetName(), trial));                          

      c->Print(Form("plots/try%i.ps%s",         
		    trial, i==c2List->LastIndex()? "]" : ""));

      cout<<Form("plots/try%i.ps%s",         
		 trial, i==c2List->GetEntries()-1 ? "]" : "")<<endl;
      nPrinted++;
    }
  }
  
  TString base = Form("plots/try%i", trial); 
  TString cmd = "ps2pdf " + base + ".ps " + base + ".pdf";             
  gSystem->Exec(cmd.Data()); 

  
  return;
}  
 
