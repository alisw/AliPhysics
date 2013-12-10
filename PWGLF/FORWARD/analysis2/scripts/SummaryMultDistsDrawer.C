#include "SummaryDrawer.C"
#include <TPaveText.h>
#include <TMultiGraph.h>

class SummaryMultDistsDrawer : public SummaryDrawer
{
public:
  enum { 
    kNormal    = 0xF
  };
  SummaryMultDistsDrawer() 
    : SummaryDrawer()
  {}
  //____________________________________________________________________
  void Run(const char* fname="forward_multdists.root", UShort_t flags=kNormal)
  {
    // --- Open the file -----------------------------------------------
    TString filename(fname);
    TFile* file = TFile::Open(filename.Data(), "READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return;
    }
    fPause         = flags & kPause;
    
    // --- Make our canvas ---------------------------------------------
    TString pdfName(filename);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, flags & kLandscape);

    // --- Force MB for pp ---------------------------------------------
    TCollection* c   = GetCollection(file, "ForwardMultSums");

    // --- Make a Title page -------------------------------------------
    DrawTitlePage(c);

    // --- Overview plots ----------------------------------------------
    fBody->Divide(1,3);
    DrawInPad(fBody, 1, GetH1(c, "triggers"),   "hist text30");
    DrawInPad(fBody, 2, GetH1(c, "status"),     "hist text30");
    DrawInPad(fBody, 3, GetH1(c, "diagnostics"),"colz text");
    PrintCanvas("Overview");

    DrawSumCollection(c, "symmetric");
    DrawSumCollection(c, "negative");
    DrawSumCollection(c, "positive");
    DrawSumCollection(c, "other");

    c   = GetCollection(file, "ForwardMultResults");
    if (!c) {
      CloseCanvas();
      return;
    }
    
    DrawResCollection(c, "symmetric");
    DrawResCollection(c, "negative");
    DrawResCollection(c, "positive");
    DrawResCollection(c, "other");

    CloseCanvas();
  }

protected:
  TCollection* GetEtaBin(TObject* o, Double_t& etaMin, Double_t& etaMax)
  {
    const char* re = "[pm][0-9]*d[0-9]*_[pm][0-9]*d[0-9]*";
    TRegexp     check(re);

    if (!o->IsA()->InheritsFrom(TCollection::Class())) {
      // Warning("GetEtaBin", "Don't know how to deal with %s - a %s",
      //         o->GetName(), o->ClassName());
      return 0;
    }
    TString oN(o->GetName());
    if (oN.Index(check) == kNPOS) { 
      // Warning("GetEtaBin", "Collection %s does not match %s",
      //         oN.Data(), re);
      return 0;
    }
    Int_t    ul     = oN.Index("_");
    TString  sMin   = oN(0, ul);
    TString  sMax   = oN(ul+1, oN.Length()-ul-1);
    sMin.ReplaceAll("p", "+");
    sMin.ReplaceAll("m", "-");
    sMin.ReplaceAll("d", ".");
    sMax.ReplaceAll("p", "+");
    sMax.ReplaceAll("m", "-");
    sMax.ReplaceAll("d", ".");
    etaMin = sMin.Atof();
    etaMax = sMax.Atof();

    return static_cast<TCollection*>(o);
  }
  //____________________________________________________________________
  void DrawSumCollection(TCollection* top, const TString& name)
  {
    TCollection* c = GetCollection(top, name, false);
    if (!c) return;

    PrintCanvas(Form("Sums - %s", name.Data()));
    
    TIter       next(c);
    TObject*    o = 0;
    while ((o = next())) { 
      Double_t etaMin = 999;
      Double_t etaMax = 999;
      TCollection* bin = GetEtaBin(o, etaMin, etaMax);
      if (!bin) continue;

      fBody->Divide(2,2);
      DrawInPad(fBody, 1, GetH1(bin, "rawDist"), "",          kLogy);
      DrawInPad(fBody, 1, GetH1(bin, "truthAccepted",false),"same", 
		kLogy|kSilent);
      DrawInPad(fBody, 1, GetH1(bin, "truth",false),   "same", 
		kLogy|kLegend|kSilent);
      DrawInPad(fBody, 2, GetH1(bin, "coverage"));
      DrawInPad(fBody, 3, GetH2(bin, "corr"),     "colz");
      DrawInPad(fBody, 4, GetH2(bin, "response",false), "colz",kLogz|kSilent);
      
      PrintCanvas(Form("%+5.1f < #eta < %+5.1f", etaMin, etaMax));
    }
  }
  //____________________________________________________________________
  void DrawResCollection(TCollection* top, const TString& name)
  {
    TCollection* c = GetCollection(top, name, false);
    if (!c) return;

    THStack* s = GetStack(c, "all");
    s->SetTitle("");
    DrawInPad(fBody, 0, s, "nostack", kLogy);
    TLegend* l = new TLegend(.5, .75, .98, .98, "P(#it{N}_{ch})");
    l->SetBorderSize(0);
    // l->SetBorderMode(0);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    TIter next(s->GetHists());
    TH1*  h = 0;
    Bool_t hasTrue = false;
    while ((h = static_cast<TH1*>(next()))) { 
      TString n(h->GetTitle());
      if (n.BeginsWith("True")) { hasTrue = true; continue; }
      n.ReplaceAll("Raw P(#it{N}_{ch}) in ", "");
      TLegendEntry* e = l->AddEntry("dummy", n, "p");
      e->SetMarkerStyle(h->GetMarkerStyle());
    }
    if (hasTrue) {
      TLegendEntry* e = l->AddEntry("dummy", "Raw", "p");
      e->SetMarkerStyle(20);
      e->SetMarkerColor(kRed+1);
      e = l->AddEntry("dummy", "MC truth", "p");
      e->SetMarkerStyle(24);
      e->SetMarkerColor(kBlue+1);
      e = l->AddEntry("dummy", "MC truth selected", "p");
      e->SetMarkerStyle(24);
      e->SetMarkerColor(kOrange+1);
    }
    fBody->cd();
    l->Draw();
    
    PrintCanvas(Form("%s results", name.Data()));

    // return;

    TIter       nextO(c);
    TObject*    o = 0;
    while ((o = nextO())) { 
      Double_t etaMin = 999;
      Double_t etaMax = 999;
      TCollection* bin = GetEtaBin(o, etaMin, etaMax);
      if (!bin) continue;

      fBody->Divide(2,3);
      DrawInPad(fBody, 1, GetH1(bin, "rawDist"),      "",     kLogy);
      DrawInPad(fBody, 1, GetH1(bin, "truthAccepted", false),
		"same", kSilent);
      DrawInPad(fBody, 1, GetH1(bin, "truth", false),"same", kSilent|kLegend);
      DrawInPad(fBody, 2, GetH1(bin, "coverage"));
      DrawInPad(fBody, 3, GetH2(bin, "corr"),            "colz");
      DrawInPad(fBody, 4, GetH2(bin, "response", false), "colz", kLogz|kSilent);
      DrawInPad(fBody, 5, GetH1(bin, "triggerVertex", false), "", kSilent);
      
      PrintCanvas(Form("%+5.1f < #eta < %+5.1f", etaMin, etaMax));
    }
  }
  //____________________________________________________________________
  void DrawTitlePage(const TCollection* c)
  {
    fBody->cd();
    
    Double_t y = .7;
    TLatex* ltx = new TLatex(.5, y, "AOD #rightarrow P(#it{N}_{ch} )");
    ltx->SetTextSize(0.07);
    ltx->SetTextFont(62);
    ltx->SetTextAlign(22);
    ltx->SetNDC();
    ltx->Draw();

    Bool_t mc = GetObject(c, "mcVertex", false) != 0;
    if (mc) { 
      y -= 0.05;
      TLatex* sub = new TLatex(.5, y, "(Simulation input)");
      sub->SetTextSize(0.04);
      sub->SetTextFont(42);
      sub->SetTextAlign(22);
      sub->SetNDC();
      sub->Draw();
    }

    Double_t save = fParName->GetTextSize();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);
    y = .6;
    
    UShort_t sys;
    UShort_t sNN;
    ULong_t  trig;
    Double_t minIpZ;
    Double_t maxIpZ;
    GetParameter(c, "sys",     sys);
    GetParameter(c, "sNN",     sNN);
    GetParameter(c, "trigger", trig);
    GetParameter(c, "minIpZ",  minIpZ);
    GetParameter(c, "maxIpZ",  maxIpZ);
    
    TString tT; TriggerString(trig, tT); DrawParameter(y, "Trigger", tT);
    TString tS; SysString(sys, tS);      DrawParameter(y, "System", tS);
    TString tE; SNNString(sNN, tE);      DrawParameter(y, "#sqrt{s_{NN}}", tE);
    DrawParameter(y, "IP_{z} range", Form("%+5.2fcm - %+5.2fcm", 
					  minIpZ, maxIpZ));
						       
    PrintCanvas("Title page");
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);

  }
};
//
// EOF
//
