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
 * Get an object from parent list 
 * 
 * @param l    Parent list 
 * @param name Name of object to find 
 * 
 * @return 0 in case of problems, object pointer otherwise 
 */
TObject* GetObject(const TDirectory* l, const char* name)
{
  if (!l) {
    Warning("GetObject", "No parent directory given");
    return 0;
  }
  TObject* o = l->Get(name);
  if (!o) { 
    Warning("GetObject", "Object named %s not found in directory %s", 
	    name, l->GetName());
    return 0;
  }
  return o;
}
Bool_t CheckType(const TObject* o, const TClass* cl)
{
  if (!o) return false;
  if (!o->IsA()->InheritsFrom(cl)) { 
    Warning("CheckType", "Object %s is not a %s, but a %s", 
	    o->GetName(), cl->GetName(), o->ClassName());
    return false;
  }
  return true;
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
  if (!CheckType(o, TList::Class())) return 0;
  return static_cast<TList*>(o);
}
/** 
 * Get a list from another list 
 * 
 * @param l    Parent list 
 * @param name Name of list to find 
 * 
 * @return 0 in case of problems, object pointer otherwise 
 */
TList* GetList(const TDirectory* l, const char* name)
{
  TObject* o = GetObject(l, name);
  if (!CheckType(o, TList::Class())) return 0;
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
  if (!CheckType(o, TH1::Class())) return 0;
  return static_cast<TH1*>(o);
} 	  
/** 
 * Get a histogram from a list 
 * 
 * @param l    Parent list 
 * @param name Name of histogram to find 
 * 
 * @return 0 in case of problems, object pointer otherwise 
 */
TH1* GetHist(const TDirectory* l, const char* name)
{
  TObject* o = GetObject(l, name);
  if (!CheckType(o, TH1::Class())) return 0;
  return static_cast<TH1*>(o);
} 	  
/** 
 * Get a histogram from a list 
 * 
 * @param l    Parent list 
 * @param name Name of histogram to find 
 * 
 * @return 0 in case of problems, object pointer otherwise 
 */
TGraphErrors* GetGraph(const TList* l, const char* name)
{
  TObject* o = GetObject(l, name);
  if (!CheckType(o, TGraphErrors::Class())) return 0;
  return static_cast<TGraphErrors*>(o);
} 	  
/** 
 * Get a histogram from a list 
 * 
 * @param l    Parent list 
 * @param name Name of histogram to find 
 * 
 * @return 0 in case of problems, object pointer otherwise 
 */
TGraphErrors* GetGraph(const TDirectory* l, const char* name)
{
  TObject* o = GetObject(l, name);
  if (!CheckType(o, TGraphErrors::Class())) return 0;
  return static_cast<TGraphErrors*>(o);
} 	  

TGraphErrors* Ratio(const TGraphErrors* num, const TGraphErrors* denom)
{
  TGraphErrors* r = static_cast<TGraphErrors*>(num->Clone("tmp"));
  for (Int_t k = 0; k < r->GetN(); k++) { 
    Double_t xn = r->GetX()[k];
    Double_t yn = r->GetY()[k];
    Double_t en = r->GetErrorY(k);
    // Double_t xd = denom->GetX()[k];
    Double_t yd = denom->GetY()[k];
    Double_t ed = denom->GetErrorY(k);
    Double_t s  = (yn-yd)/yd;
    Double_t e  = 0;
    r->SetPoint(k, xn, s);
    r->SetPointError(k, 0, e);
  }
  return r;
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
GenerateEmpirical(const char* fmdDispVtx="forward_dndetaDispVtx.root",
		  const char* fmdNomVtx="forward_dndetaNominalVtxZDC.root",
		  const char* fullDispVtx="dndetaPreliminaryQM12.root")
{
  
  TString outName    = "EmpiricalCorrection.root";
  TFile* fdispvtxFMD = OpenFile(fmdDispVtx);
  TFile* fNominalFMD = OpenFile(fmdNomVtx);
  TFile* fdispvtxAll = OpenFile(fullDispVtx);
  if (!fdispvtxFMD || !fNominalFMD || !fdispvtxAll) return;

  TFile* out = TFile::Open(outName, "RECREATE");
  if (!out) { 
    Error("GenerateEmpirical", "Failed to open output file %s", 
	  outName.Data());
    return;
  }

  Int_t limits[] = {     0,        5,      10,         20,    30  };
  Int_t colors[] = {kRed+2, kGreen+2, kBlue+1, kMagenta+1, kBlack };

  TGraphErrors* afmdfmd  = 0;
  TGraphErrors* afmdfull = 0;

  // --- Do two things fmd-to-fmd and fmd-to-full
  for (Int_t iM = 0; iM < 2; iM++) {
    Bool_t       fmdfmd   = (iM == 0);
    TMultiGraph* mg       = new TMultiGraph();
    TMultiGraph* mgfmdnom = (iM == 0 ? new TMultiGraph() : 0);
    TMultiGraph* mgref    = new TMultiGraph();
    mgref->SetTitle(Form("Satellite, FMD%s", fmdfmd ? "" : "+V0+Tracklets"));
    if (mgfmdnom) mgfmdnom->SetTitle("Nominal FMD");

    Info("", "Doing %s/FMD", fmdfmd ? "FMD" : "Full");
    TGraphErrors* corrcent[4];    
    TGraphErrors* allsym     = 0;
    TGraphErrors* fmddispvtx = 0;
    TGraphErrors* fmdnominal = 0;

    TList* displist = GetList(fdispvtxFMD,"ForwardResults");
    TList* nomlist  = GetList(fNominalFMD,"ForwardResults");
    
    Int_t nMin = 1000;
    // --- Get each defined centrality and form the ratio ------------
    for(Int_t iC=0; iC<4; iC++) {
      Int_t cl = limits[iC];
      Int_t ch = limits[iC+1];
      TString base; base.Form("_cent_%03d_%03d", cl, ch);
      corrcent[iC] = new TGraphErrors;
      corrcent[iC]->SetName(Form("correction%s",base.Data()));
      corrcent[iC]->SetTitle(Form("%2d%% - %2d%%", cl, ch));
      corrcent[iC]->GetHistogram()->SetXTitle("#eta");
      corrcent[iC]->GetHistogram()->SetYTitle("#frac{dN_{ch}/d#eta|_{nominal}}{"
					     "dN_{ch}/d#eta|_{satelitte}}");  
      corrcent[iC]->SetMarkerColor(colors[iC]);
      corrcent[iC]->SetLineColor(colors[iC]);
      corrcent[iC]->SetMarkerStyle(fmdfmd ? 20 : 24);
      corrcent[iC]->SetMarkerSize(fmdfmd ? 1 : 1.2);
      corrcent[iC]->SetFillStyle(0);
      corrcent[iC]->SetFillColor(0);
      
      mg->Add(corrcent[iC]);
      TGraphErrors* allsym = GetGraph(fdispvtxAll,Form("graphErrors_cent_%d_%d",
						       cl, ch));
    
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
      if (mgfmdnom) {
	TGraph* g = static_cast<TGraph*>(fmdnominal->Clone(Form("nominal%s", 
								base.Data())));
	g->SetMarkerColor(colors[iC]);
	g->SetLineColor(colors[iC]);
	g->SetMarkerStyle(21);
	g->SetTitle("Nominal FMD");
	mgfmdnom->Add(g);
      }
      TGraph* ref = static_cast<TGraph*>(fmddispvtx->Clone(Form("satellite%s",
								base.Data())));
      ref->SetTitle(Form("%2d - %2d", cl, ch));
      ref->SetMarkerColor(colors[iC]);
      ref->SetLineColor(colors[iC]);
      ref->SetMarkerStyle(fmdfmd ? 22 : 20);
      ref->SetMarkerSize(1);
      mgref->Add(ref);
    
      Int_t nPoints = 0;
    
      for(Int_t iN=0; iN<fmdnominal->GetN(); iN++) {
	Double_t eta        = fmdnominal->GetX()[iN];
	Double_t nommult    = fmdnominal->GetY()[iN];
	Double_t nommulterr = fmdnominal->GetErrorY(iN);
	Double_t etaerr     = fmdnominal->GetErrorX(iN);

	// Ignore empty bins 
	if(nommult < 0.0001) continue;

	// Find the corresponding bin from the dispaclaced vertex analysis
	for(Int_t iS=0; iS < fmddispvtx->GetN(); iS++) {
	  Double_t eta1        = fmddispvtx->GetX()[iS];
	  Double_t dispmult    = fmddispvtx->GetY()[iS];
	  Double_t dispmulterr = fmddispvtx->GetErrorY(iS);

	  if(TMath::Abs(eta-eta1) >= 0.001) continue;
	  
	  Double_t corr  = nommult/dispmult;
	  Double_t rd    = dispmulterr/dispmult;
	  Double_t rn    = nommulterr/nommult;
	  Double_t error = (1/corr)*TMath::Sqrt(TMath::Power(rd,2) + 
					      TMath::Power(rn,2));
	  corrcent[iC]->SetPoint(nPoints,eta,corr);
	  corrcent[iC]->SetPointError(nPoints,etaerr,error);
	  nPoints++;
	  
	}
	//std::cout<<eta<<"  "<<nommult<<std::endl;
      }
      nMin = TMath::Min(nPoints, nMin);
    }

    // --- Calulate the average --------------------------------------
    TGraphErrors* average = new TGraphErrors;
    average->SetName(Form("average"));
    average->SetTitle(Form("%2d%% - %2d%%", limits[0], limits[4]));
    average->SetMarkerColor(colors[4]);
    average->SetLineColor(colors[4]);
    average->SetMarkerStyle(fmdfmd ? 20 : 24);
    average->SetMarkerSize(fmdfmd ? 1 : 1.2);
    average->SetFillStyle(0);
    average->SetFillColor(0);
    mg->Add(average);
    
    for(Int_t iA = 0; iA < nMin; iA++) {
      Double_t mean   = 0;
      Double_t error2 = 0;
      Double_t sumw2  = 0;
      
      // Loop over centralities 
      for(Int_t iC=0; iC<4; iC++) {
	Double_t eta  = corrcent[iC]->GetX()[iA];
	Double_t corr = corrcent[iC]->GetY()[iA];
	Double_t err  = corrcent[iC]->GetErrorY(iA);
	Double_t err2 = err * err;
	Double_t w    = 1 / err2;
	sumw2         += w;
	mean          += w * corr;
      }
      mean           /= sumw2;
      Double_t error = TMath::Sqrt(1. / sumw2);
          average->SetPoint(iA,eta,mean);
      average->SetPointError(iA,0.125,error);    
    }
    if (fmdfmd) afmdfmd  = average;
    else        afmdfull = average;

    // --- Calculate ratios ------------------------------------------
    TMultiGraph* ratios = new TMultiGraph;
    ratios->SetName("ratios");
    for(Int_t iC=0; iC<4; iC++) {
      Int_t cl        = limits[iC];
      Int_t ch        = limits[iC+1];
      TGraphErrors* r = Ratio(corrcent[iC], average);
      r->SetName(Form("ratio%s",base.Data()));
      ratios->Add(r);
    }

    // --- Store the result ------------------------------------------
    TDirectory* d = out->mkdir(fmdfmd ? "fmdfmd" : "fmdfull", 
			       Form("Empirical corrections %s", 
				    fmdfmd ? "FMD/FMD" : "FULL/FMD"));
    d->cd();
    for(Int_t p=0; p<4; p++) corrcent[p]->Write();
    average->Write("average");
    mg->Write("all");
    ratios->Write("ratios");
    mgref->Write("satellite");
    out->cd();
    if (mgfmdnom) mgfmdnom->Write("nominal");
    out->ls();
  } // End method loop

  TGraphErrors* r = Ratio(afmdfull, afmdfmd);
  r->Write("ratio");

  out->Close();
}
//
// EOF
//

