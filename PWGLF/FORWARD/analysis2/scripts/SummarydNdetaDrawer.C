#include "SummaryDrawer.C"
#include <TPaveText.h>

class SummarydNdetaDrawer : public SummaryDrawer
{
public:
  SummarydNdetaDrawer() 
    : SummaryDrawer()
  {}
  //____________________________________________________________________
  void Run(const char* fname="forward_dndeta.root",
	   UShort_t flags=0xf)
  {
    // --- Open the file -----------------------------------------------
    TString filename(fname);
    TFile* file = TFile::Open(filename.Data(), "READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return;
    }
    // Options 
    Bool_t forward = flags & 0x1;
    Bool_t central = flags & 0x2;
    Bool_t sums    = flags & 0x4;
    Bool_t results = flags & 0x8;
    Bool_t onlyMB  = flags & 0x10;
    fPause         = flags & 0x40;
    
    
    // --- Make our canvas ---------------------------------------------
    TString pdfName(filename);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, flags & 0x20);

    // --- Do each sub-algorithm ---------------------------------------
    THStack* rF = 0;
    if (forward && sums)    DrawSums(file, "Forward", onlyMB);
    if (forward && results) rF = DrawRes(file, "Forward", onlyMB);
    
    THStack* rC = 0;
    if (central && sums)    DrawSums(file, "Central", onlyMB);
    if (central && results) rC = DrawRes(file, "Central", onlyMB);

    if (!(rC && rF) || !results) {
      CloseCanvas();
      return;
    }

    fBody->cd();
    Double_t y1 = fLandscape ? 0  : .3;
    Double_t x2 = fLandscape ? .7 : 1;
    Double_t x1 = fLandscape ? x2 : 0;
    Double_t y2 = fLandscape ? 1  : y1;
    TPad* p1 = new TPad("p1", "p1", 0,  y1, x2, 1,  0, 0);
    TPad* p2 = new TPad("p2", "p2", x1, 0,  1,  y2, 0, 0);

    fBody->cd();
    p1->Draw();
    p1->cd();

    TIter next(rF->GetHists());
    TH1*  h  = 0;
    while ((h = static_cast<TH1*>(next()))) rC->Add(h);
    
    rC->Draw("nostack");

    fBody->cd();
    p2->Draw();
    p2->cd();

    TLegend* l = new TLegend(0.01, 0.1, 0.99, 0.99, "Centralities");
    l->SetNColumns(fLandscape ? 1 : 2);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    TAxis* centAxis = 
      static_cast<TAxis*>(GetObject(GetCollection(file, "ForwardSums"), 
				    "centAxis"));
    CleanStack(rC, l, centAxis);
    l->Draw();

    PrintCanvas("Both");
  
    CloseCanvas();
  }

protected:
  //____________________________________________________________________
  void DrawSums(TDirectory* top, const TString& base, bool onlyMB)
  {
    TCollection* c = GetCollection(top, Form("%sSums", base.Data()));
    if (!c) return;
    
    TAxis* centAxis = static_cast<TAxis*>(GetObject(c, "centAxis"));
    
    fBody->Divide(1, 2);
    
    fBody->cd(1);
    Double_t save  = fParName->GetTextSize();
    Double_t xSave = fParVal->GetX();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);
    fParVal->SetX(.4);
    Double_t y = .8;
    for (Int_t i = 1; i <= centAxis->GetNbins(); i++) { 
      DrawParameter(y, (i == 1 ? "Centrality classes" : ""),
		    Form("%3d%% - %3d%%", 
			 Int_t(centAxis->GetBinLowEdge(i)), 
			 Int_t(centAxis->GetBinUpEdge(i))));
    }
    Int_t sys, sNN, scheme, trigger;
    GetParameter(c, "sNN",     sNN); 
    GetParameter(c, "sys",     sys); 
    GetParameter(c, "scheme",  scheme); 
    GetParameter(c, "trigger", trigger); 
    DrawParameter(y, "Collision system", (sys == 1 ? "pp" : 
					  (sys == 2 ? "PbPb" : 
					   (sys == 3 ? "pPb" : 
					    "unknown"))));
    DrawParameter(y, "#sqrt{s_{NN}}", Form("%4dGeV", sNN));
    DrawParameter(y, "Normalization scheme", Form("0x%x", scheme));
    DrawParameter(y, "Triggers",      Form("0x%x", trigger));
  
    TH1* cent = GetH1(c, "cent");
    cent->SetFillColor(kRed+1);
    cent->SetFillStyle(3001);
    cent->SetXTitle("Centrality [%]");
    cent->SetYTitle("Events");

    DrawInPad(fBody, 2, cent);
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);
    fParVal->SetX(xSave);

    PrintCanvas(Form("%s sums", base.Data()));

    DrawCentSum(c, base, centAxis->GetXmin(), -centAxis->GetXmax());
    if (onlyMB) return;

    for (Int_t i = 1; i <= centAxis->GetNbins(); i++) 
      DrawCentSum(c, base, centAxis->GetBinLowEdge(i), 
		  centAxis->GetBinUpEdge(i));
  }
  //____________________________________________________________________
  void DrawCentSum(const TCollection* sums, const TString& base, 
		   Int_t cLow, Int_t cHigh)
  {
    TString folder; 
    if (cLow < 0 || cHigh < 0 || cLow >= cHigh) {
      folder = "all";
      cHigh  *= -1;
    }
    else folder.Form("cent%03d_%03d", cLow, cHigh);
    
    TCollection* c = GetCollection(sums, folder);
    if (!c) return;
    
    fBody->Divide(2, 2);

    TH1* type = GetH1(c, Form("%sEvents",base.Data()));
    type->SetFillStyle(3001);
    type->SetFillColor(kGreen+1);
    TH2* bin0 = GetH2(c, Form("%s0", base.Data()));
    
    DrawInPad(fBody, 1, GetH1(c, "triggers"), "HIST TEXT");
    DrawInPad(fBody, 2, type,                 "HIST TEXT");
    DrawInPad(fBody, 3, GetH2(c, base.Data()),"colz");
    DrawInPad(fBody, 4, bin0,                 "colz");
    
    if (bin0->GetEntries() <= 0) {
      fBody->cd(4);
      TLatex* l = new TLatex(0.5, 0.5, "No 0-bin events");
      l->SetNDC();
      l->SetTextAlign(22);
      l->Draw();
    }
    PrintCanvas(Form("%s sums: %3d%% - %3d%%", base.Data(), cLow, cHigh));
  }
  //____________________________________________________________________
  THStack* DrawRes(TDirectory* top, const TString& base, Bool_t onlyMB)
  {
    TCollection* c = GetCollection(top, Form("%sResults", base.Data()));
    if (!c) return 0;

    TAxis* centAxis = static_cast<TAxis*>(GetObject(c, "centAxis"));
    
    fBody->cd();

    Double_t save = fParName->GetTextSize();
    Double_t xSave = fParVal->GetX();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);
    fParVal->SetX(.4);
    Double_t y = .9;
    for (Int_t i = 1; i <= centAxis->GetNbins(); i++) { 
      DrawParameter(y, (i == 1 ? "Centrality classes" : ""),
		  Form("%3d%% - %3d%%", 
		       Int_t(centAxis->GetBinLowEdge(i)), 
		       Int_t(centAxis->GetBinUpEdge(i))));
    }
    
    DrawParameter(y, "Collision system", GetObject(c, "sys")->GetTitle());
    DrawParameter(y, "#sqrt{s_{NN}}",GetObject(c,"sNN")->GetTitle());
    DrawParameter(y, "trigger",GetObject(c,"trigger")->GetTitle());
    DrawParameter(y, "scheme", GetObject(c,"scheme")->GetTitle());

    Double_t epsT, epsT0;
    GetParameter(c, "triggerEff",  epsT);
    GetParameter(c, "triggerEff0", epsT0);
    DrawParameter(y, "#epsilon_{T}", Form("%f", epsT));
    DrawParameter(y, "#epsilon_{T,zero bin}", Form("%f", epsT0));
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);
    fParVal->SetX(xSave);

    PrintCanvas(Form("%s results", base.Data()));
  
    fBody->Divide(1, 3, 0, 0);

    TLegend* l = new TLegend(0.1, 0.1, 0.9, 0.9, "Centralities");
    l->SetNColumns(2);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    THStack* dndeta_  = GetStack(c, "dndeta");
    THStack* dndeta   = CleanStack(dndeta_, l, centAxis);
    THStack* dndeta5_ = GetStack(c, "dndeta_rebin05");
    THStack* dndeta5  = CleanStack(dndeta5_, 0, 0);
    
    DrawInPad(fBody, 1, l, "");
    DrawInPad(fBody, 2, dndeta,  "nostack");
    DrawInPad(fBody, 3, dndeta5, "nostack");
    fBody->GetPad(2)->SetGridx();
    fBody->GetPad(3)->SetGridx();

    PrintCanvas(Form("%s results - stacks", base.Data()));
  
    DrawCentRes(c, base, centAxis->GetXmin(), -centAxis->GetXmax());
    if (onlyMB) return dndeta;

    for (Int_t i = 1; i <= centAxis->GetNbins(); i++) 
      DrawCentRes(c, base, centAxis->GetBinLowEdge(i), 
		  centAxis->GetBinUpEdge(i));
    
    return dndeta;
  }

  //____________________________________________________________________
  void DrawCentRes(const TCollection* sums, const TString& base, 
		   Int_t cLow, Int_t cHigh)
  {
    TString folder; 
    if (cLow < 0 || cHigh < 0 || cLow >= cHigh) {
      folder =  "all";
      cHigh  *= -1;
    }
    else folder.Form("cent%03d_%03d", cLow, cHigh);
  
    TCollection* c = GetCollection(sums, folder);
    if (!c) return;

    fBody->Divide(2, 3, 0.05, 0);

    Int_t        trP = 1;
    TVirtualPad* p   = fBody->GetPad(trP);
    p->SetBottomMargin(0.15);
    p->SetLeftMargin(0.15);
    if (trP > 2) p->SetTopMargin(0.05);
  
    TH1* norm = GetH1(c, Form("norm%s",base.Data()));
    norm->SetFillColor(kGreen+1);
    norm->SetFillStyle(3001);
    
    DrawInPad(fBody, trP, GetH1(c, "triggers"), "HIST TEXT");
    DrawInPad(fBody, 2,   norm);
    DrawInPad(fBody, 4,   GetH1(c, Form("dndeta%s",base.Data())));
    DrawInPad(fBody, 6,   GetH1(c, Form("dndeta%s_rebin05",base.Data())));
    DrawInPad(fBody, 5,   GetH2(c, Form("d2Ndetadphi%s", base.Data())),"colz");
  
    fBody->GetPad(2)->SetGridx(); fBody->GetPad(2)->SetLeftMargin(0.15);
    fBody->GetPad(4)->SetGridx(); fBody->GetPad(4)->SetLeftMargin(0.15);
    fBody->GetPad(6)->SetGridx(); fBody->GetPad(6)->SetLeftMargin(0.15);

    TObject*   normCalc = GetObject(c, "normCalc");
    TString    calc     = normCalc->GetTitle();
    TObjArray* lines    = calc.Tokenize("\n");
    TPaveText* disp     = new TPaveText(.1,.1,.9,.9, "NDC");
    TIter      next(lines);
    TObject*   line     = 0;
    while ((line = next())) disp->AddText(line->GetName());
    disp->SetBorderSize(0);
    disp->SetBorderSize(0);
    disp->SetFillStyle(0);
    DrawInPad(fBody, 3, disp);

    PrintCanvas(Form("%s result: %3d%% - %3d%%", base.Data(), cLow, cHigh));
  }  
  //____________________________________________________________________
  THStack* CleanStack(const THStack* stack, TLegend* l, const TAxis* axis)
  {
    THStack* ret   = new THStack(stack->GetName(), stack->GetTitle());
    TList*   hists = stack->GetHists();
    TIter    next(hists);
    TH1*     h     = 0;
    Int_t    j     = 0;
    Bool_t   ok    = false;
    while ((h = static_cast<TH1*>(next()))) {
      TString name(h->GetName());
      if (name.Contains("_mirror")) continue;
      if (l && !ok) { 
	j++;
	name.Form("%3d%% - %3d%%", 
		  Int_t(axis->GetBinLowEdge(j)), 
		  Int_t(axis->GetBinUpEdge(j)));
	ok = axis->GetBinUpEdge(j) > 100;
	TLegendEntry* e = l->AddEntry("dummy", name, "f");
	e->SetFillStyle(1001);
	e->SetFillColor(h->GetMarkerColor());
      }
      ret->Add(h);
    }
    return ret;
  }
};
//
// EOF
//
