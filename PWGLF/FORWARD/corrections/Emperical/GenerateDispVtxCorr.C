/** 
 * Open a file 
 * 
 * @param filename Name of file 
 * 
 * @return Return value 
 */
TFile* OpenFile(const char* filename)
{
  TFile* file = TFile::Open(filename, "READ");
  if (!file) 
    Error("OpenFile", "Failed to open the file %s for reading", filename);
  return file;
}
/** 
 * Get an object from parent list 
 * 
 * @param l    Parent list 
 * @param name Name of object to find 
 * 
 * @return 0 in case of problems, object pointer otherwise 
 */
TObject* GetObject(const TList* l, const char* name)
{
  if (!l) {
    Warning("GetObject", "No parent list given");
    return 0;
  }
  TObject* o = l->FindObject(name);
  if (!o) { 
    Warning("GetObject", "Object named %s not found in list %s", 
	    name, l->GetName());
    return 0;
  }
  return o;
}

/** 
 * Get a list from another list 
 * 
 * @param l    Parent list 
 * @param name Name of list to find 
 * 
 * @return 0 in case of problems, object pointer otherwise 
 */
TList* GetList(const TList* l, const char* name)
{
  TObject* o = GetObject(l, name);
  if (!o) return 0;
  return static_cast<TList*>(o);
}
/** 
 * Get a histogram from a list 
 * 
 * @param l    Parent list 
 * @param name Name of histogram to find 
 * 
 * @return 0 in case of problems, object pointer otherwise 
 */
TH1* GetHist(const TList* l, const char* name)
{
  TObject* o = GetObject(l, name);
  if (!o) return 0;
  return static_cast<TH1*>(o);
} 	  
/** 
 * Generate the emperical correction for unknown material in the
 * AliROOT geometric description.  This is done by taking the ratio
 * of the nominal vertex results to the results from satelitte events. 
 *
 * This can be done in two ways 
 *
 * - Either by using FMD results only 
 * - Or by using the full combined FMD+V0+SPD results 
 * 
 * The correction is written to the file "EmpiricalCorrection.root",
 * which will contain a TGraphErrors for each defined centrality
 * class, plus a TGraphErrors object that contains the average over
 * all centrality classes.
 *
 * Other results should be DIVIDED by the correction obtained from
 * this script.
 * 
 * @param fmdfmd If true, use FMD results only 
 */
void 
GenerateDispVtxCorr(Bool_t fmdfmd=true,
		    const char* fmdDispVtx="forward_dndetaDispVtx.root",
		    const char* fmdNomVtx="forward_dndetaNominalVtxZDC.root",
		    const char* fullDispVtx="dndetaPreliminaryQM12.root")
{
  
  TString outName    = "EmpiricalCorrection.root";
  TFile* fdispvtxFMD = OpenFile(fmdDispVtx);
  TFile* fNominalFMD = OpenFile(fmdNomVtx);
  TFile* fdispvtxAll = OpenFile(fullDispVtx);
  if (!fdispvtxFMD || !fNominalFMD || !fdispvtxAll) return;
 
  TGraphErrors* corrcent[4];
  TMultiGraph*  mg = new TMultiGraph();
  mg->SetName("corrections");
  
  Int_t limits[] = {     0,        5,      10,         20,    30  };
  Int_t colors[] = {kRed+2, kGreen+2, kBlue+1, kMagenta+1, kBlack };

  TGraphErrors* allsym     = 0;
  TGraphErrors* fmddispvtx = 0;
  TGraphErrors* fmdnominal = 0;

  TList* displist = static_cast<TList*>(fdispvtxFMD->Get("ForwardResults"));
  TList* nomlist  = static_cast<TList*>(fNominalFMD->Get("ForwardResults"));
    
  Int_t nMin = 1000;
  // Get each defined centrality and form the ratio 
  for(Int_t i=0; i<4; i++) {
    Int_t cl = limits[i];
    Int_t ch = limits[i+1];
    corrcent[i] = new TGraphErrors;
    corrcent[i]->SetName(Form("correction_cent_%03d_%03d",cl, ch));
    corrcent[i]->SetTitle(Form("%2d%% - %2d%%", cl, ch));
    corrcent[i]->GetHistogram()->SetXTitle("#eta");
    corrcent[i]->GetHistogram()->SetYTitle("#frac{dN_{ch}/d#eta|_{nominal}}{"
					   "dN_{ch}/d#eta|_{satelitte}}");  
    corrcent[i]->SetMarkerColor(colors[i]);
    corrcent[i]->SetLineColor(colors[i]);
    corrcent[i]->SetMarkerStyle(8);
    corrcent[i]->SetFillStyle(0);
    corrcent[i]->SetFillColor(0);

    mg->Add(corrcent[i]);
    TGraphErrors* allsym = 
      static_cast<TGraphErrors*>(fdispvtxAll->Get(Form("graphErrors_cent_%d_%d",
						       cl, ch)));
    
    TString folderName   = Form("cent%03d_%03d",cl, ch);
    TList*  dispcentlist = GetList(displist, folderName);
    TList*  nomcentlist  = GetList(nomlist, folderName);
    TString histName     = "dndetaForward_rebin05";
    TH1D*   hDisp        = GetHist(dispcentlist, histName);
    TH1D*   hNominal     = GetHist(nomcentlist, histName);
    
    // Make our temporary graph 
    if   (fmdfmd) fmddispvtx = new TGraphErrors(hDisp);
    else          fmddispvtx = static_cast<TGraphErrors*>(allsym->Clone());
    fmdnominal               = new TGraphErrors(hNominal);
    
    Int_t nPoints = 0;
    
    for(Int_t n=0; n<fmdnominal->GetN(); n++) {
      Double_t eta        = fmdnominal->GetX()[n];
      Double_t nommult    = fmdnominal->GetY()[n];
      Double_t nommulterr = fmdnominal->GetErrorY(n);
      Double_t etaerr     = fmdnominal->GetErrorX(n);

      // Ignore empty bins 
      if(nommult < 0.0001) continue;

      // Find the corresponding bin from the dispaclaced vertex analysis
      for(Int_t m=0; m < fmddispvtx->GetN(); m++) {
	Double_t eta1        = fmddispvtx->GetX()[m];
	Double_t dispmult    = fmddispvtx->GetY()[m];
	Double_t dispmulterr = fmddispvtx->GetErrorY(m);

	if(TMath::Abs(eta-eta1) >= 0.001) continue;

	Double_t corr  = nommult/dispmult;
	Double_t rd    = dispmulterr/dispmult;
	Double_t rn    = nommulterr/nommult;
	Double_t error = (1/corr)*TMath::Sqrt(TMath::Power(rd,2) + 
					      TMath::Power(rn,2));
	corrcent[i]->SetPoint(nPoints,eta,corr);
	corrcent[i]->SetPointError(nPoints,etaerr,error);
	nPoints++;
      
      }
      //std::cout<<eta<<"  "<<nommult<<std::endl;
    }
    nMin = TMath::Min(nPoints, nMin);
  }

  // Calulate the average 
  TGraphErrors* average = new TGraphErrors;
  average->SetName(Form("correction_average"));
  average->SetTitle(Form("%2d%% - %2d%%", limits[0], limits[4]));
  average->SetMarkerColor(colors[4]);
  average->SetLineColor(colors[4]);
  average->SetMarkerStyle(8);
  average->SetFillStyle(0);
  average->SetFillColor(0);
  mg->Add(average);

  for(Int_t k = 0; k < nMin; k++) {
    Double_t mean   = 0;
    Double_t error2 = 0;
    Double_t sumw2  = 0;
    
    // Loop over centralities 
    for(Int_t l=0; l<4; l++) {
      Double_t eta  = corrcent[l]->GetX()[k];
      Double_t corr = corrcent[l]->GetY()[k];
      Double_t err  = corrcent[l]->GetErrorY(k);
      Double_t err2 = err * err;
      Double_t w    = 1 / err2;
      sumw2         += w;
      mean          += w * corr;
#if 0
      mean   += corr;
      error2 += TMath::Power(err,2);
#endif
    }
    mean           /= sumw2;
    Double_t error = TMath::Sqrt(1. / sumw2);
#if 0
    mean           /= 4;
    Double_t error =  TMath::Sqrt(error2) / 4;
#endif
    
    average->SetPoint(k,eta,mean);
    average->SetPointError(k,0.125,error);    
  }

  Double_t min = +1000;
  Double_t max = -1000;
  TMultiGraph* ratios = new TMultiGraph;
  for(Int_t l=0; l<4; l++) {
    Int_t cl        = limits[l];
    Int_t ch        = limits[l+1];
    TGraphErrors* r = 
      static_cast<TGraphErrors*>(corrcent[l]->Clone(Form("ratio%02d%02d",
							cl,ch)));
    ratios->Add(r);
    for (Int_t k = 0; k < r->GetN(); k++) { 
      Double_t x = r->GetX()[k];
      Double_t y = r->GetY()[k];
      Double_t a = average->Eval(x);
      Double_t s = (y-a)/a;
      r->SetPoint(k, x, s);
      min  = TMath::Min(s, min);
      max  = TMath::Max(s, max);
    }
  }
  Printf("Min=%f max=%f", min, max);
  // Draw the results 
  TCanvas* canvas = new TCanvas("overview","overview",800,600);
  TPad* overview =  new TPad("overview", "Overview", 0, .3, 1, 1, 0, 0);
  overview->SetBottomMargin(0);
  overview->SetTopMargin(0.02);
  overview->SetRightMargin(0.02);
  overview->Draw();
  overview->SetGridx();
  overview->cd();

  mg->Draw("APEL");
  mg->GetXaxis()->SetTitle("#eta");
  mg->GetYaxis()->SetTitle("#frac{dN_{ch}/d#eta|_{nominal}}{"
			   "dN_{ch}/d#eta|_{satelitte}}");  
  overview->Clear();
  mg->DrawClone("APEC");
  TLegend* leg = overview->BuildLegend(0.35,0.6,0.6,0.975, "Centrality");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  canvas->cd();
  TPad* details =  new TPad("details", "Details", 0, 0, 1, .3, 0, 0);
  details->SetBottomMargin(0.15);
  details->SetTopMargin(0);
  details->SetRightMargin(0.02);
  details->Draw();
  details->SetGridx();
  details->cd();


  ratios->Draw("APEL");
  ratios->GetXaxis()->SetTitle("#eta");
  ratios->GetYaxis()->SetTitle("#frac{c-#LTc#GT}{#LTc#GT}");  
  ratios->GetYaxis()->SetTitleSize(0.09);
  ratios->GetYaxis()->SetTitleOffset(0.5);
  ratios->GetYaxis()->SetLabelSize(0.08);
  ratios->GetXaxis()->SetTitleSize(0.09);
  ratios->GetXaxis()->SetTitleOffset(0.5);
  ratios->GetXaxis()->SetLabelSize(0.08);
  ratios->GetYaxis()->SetNdivisions(10);

  Double_t x1 = ratios->GetXaxis()->GetXmin();
  Double_t x2 = ratios->GetXaxis()->GetXmax();
  Double_t y1 = min; // ratios->GetHistogram()->GetMinimum();
  Double_t y2 = max; // ratios->GetHistogram()->GetMaximum();
  TGraphAsymmErrors* band = new TGraphAsymmErrors(2);
  band->SetName("band");
  band->SetTitle(Form("Error band min=%4.2f%% max=%4.2f%%", -100*min, 100*max));
  band->SetPoint(0, x1, 0);
  band->SetPoint(1, x2, 0);
  band->SetPointError(0, 0, 0, -y1, y2);
  band->SetPointError(1, 0, 0, -y1, y2);
  band->SetFillColor(kYellow+2);
  band->SetFillStyle(3001);
  band->SetLineStyle(2);
  band->SetLineColor(kBlack);
  band->Draw("A E3 L");
  band->GetXaxis()->SetTitle("#eta");
  band->GetYaxis()->SetTitle("#frac{c-#LTc#GT}{#LTc#GT}");  
  band->GetYaxis()->SetTitleSize(0.09);
  band->GetYaxis()->SetTitleOffset(0.5);
  band->GetYaxis()->SetLabelSize(0.08);
  band->GetXaxis()->SetTitleSize(0.09);
  band->GetXaxis()->SetTitleOffset(0.5);
  band->GetXaxis()->SetLabelSize(0.08);
  band->GetYaxis()->SetNdivisions(10);

  details->Clear();
  // band->DrawClone("A E3 L");
  ratios->GetListOfGraphs()->AddFirst(band, "E3 L");
  ratios->DrawClone("PEC same");

  leg = details->BuildLegend(0.35,0.6,0.6,0.975, "Centrality");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  // Store the results 
  TFile* fout = TFile::Open(outName,"RECREATE");

  for(Int_t p=0; p<4; p++) corrcent[p]->Write();
  average->Write("average");
  mg->Write("all");
  ratios->Write("all");
  fout->ls();
  fout->Close();
  Info("GenerateDispVtxCorr", "Corrections written to %s", outName.Data());
}
//
// EOF
//
