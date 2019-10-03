#include "SummaryDrawer.C"
#include <TPaveText.h>
#include <TMultiGraph.h>

class SummarydNdetaDrawer : public SummaryDrawer
{
public:
  enum { 
    kForward   = 0x01,
    kCentral   = 0x02, 
    kSums      = 0x04, 
    kResults   = 0x08, 
    kMinBias   = 0x10, 
    kMC        = 0x80,
    kNormal    = 0x0F
  };
  SummarydNdetaDrawer() 
    : SummaryDrawer()
  {}
  const char* ColName(const char* prefix, bool results=false)
  {
    return Form("%sdNdeta%s", prefix, results ? "Results" : "Sums");
  }
  //____________________________________________________________________
  void Run(const char* fname="forward_dndeta.root", UShort_t flags=kNormal)
  {
    // --- Open the file -----------------------------------------------
    TString filename(fname);
    TFile* file = TFile::Open(filename.Data(), "READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return;
    }
    // Options 
    Bool_t forward = flags & kForward;
    Bool_t central = flags & kCentral;
    Bool_t sums    = flags & kSums;
    Bool_t results = flags & kResults;
    Bool_t onlyMB  = flags & kMinBias;
    Bool_t mc      = flags & kMC;
    fPause         = flags & kPause;
    
    // --- Force MB for pp ---------------------------------------------
    UShort_t     sys = 0;
    TCollection* c   = GetCollection(file, ColName("Forward", false));
    GetParameter(c, "sys", sys); 
    if (sys == 1) {
      // onlyMB = true;
      // Info("Run", "Found sys==1 -> Forcing MB");
    }

    // --- Test of MC --------------------------------------------------
    TCollection* mcC = GetCollection(file, ColName("MCTruth"), false);
    if (mcC) { 
      TCollection* mcAll = GetCollection(mcC, "all");
      if (mcAll && GetObject(mcAll, "sum")) {
	Info("Run", "Found MC truth output");
	mc = true;
      }
    }
    // --- Make our canvas ---------------------------------------------
    TString pdfName(filename);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, flags & kLandscape);

    // --- Make a Title page -------------------------------------------
    DrawTitlePage(file, mc, onlyMB);

    // --- Do each sub-algorithm ---------------------------------------
    THStack* rF = 0;
    if (forward && sums)    DrawSums(file, "Forward", onlyMB);
    if (forward && results) rF = DrawRes(file, "Forward", onlyMB);
    
    THStack* rC = 0;
    if (central && sums)    DrawSums(file, "Central", onlyMB);
    if (central && results) rC = DrawRes(file, "Central", onlyMB);

    THStack* rM = 0;
    if (mc && sums)    DrawSums(file, "MCTruth", onlyMB);
    if (mc && results) rM = DrawRes(file, "MCTruth", onlyMB);

    if (rC && rF && results) DrawBoth(rC, rF, rM, file, onlyMB);
  
    CloseCanvas();
  }

protected:
  //____________________________________________________________________
  TAxis* GetCentAxis(const TCollection* parent, Bool_t verbose=false)
  {
    TObject* cO = GetObject(parent, "centAxis", verbose);
    TAxis*   cA = 0;
    if (!cO) return 0;

    if (cO->IsA()->InheritsFrom(TAxis::Class())) 
      cA = static_cast<TAxis*>(cO);
    else if (cO->IsA()->InheritsFrom(TH1::Class())) {
      TH1*  cH = static_cast<TH1*>(cO);
      cA = cH->GetXaxis();
    }
    // if (cA) cA->Dump();
    if (!cA || !cA->GetXbins() || !cA->GetXbins()->GetArray() ||
	cA->GetXmin() > cA->GetXmax()) return 0;
    return cA;
  }
  //____________________________________________________________________
  TCollection* GetCentCollection(const TCollection* sums, 
				 const TString&     base, 
				 Float_t            cLow, 
				 Float_t            cHigh,
				 TString&           title)
  {
    TString folder; 
    title = TString::Format("%s %s: ", base.Data(), title.Data());
    if (cLow < 0 || cHigh < 0 || cLow >= cHigh) {
      folder = "all";
      title.Append("All selected events");
    }
    else {
      folder.Form("cent%03dd%02d_%03dd%02d",
		  Int_t(cLow), Int_t(cLow*100) % 100,
		  Int_t(cHigh), Int_t(cHigh*100) % 100);
      title.Append(Form("%6.2f%% - %6.2f%%", cLow, cHigh));
    }
    
    return GetCollection(sums, folder);
    
  }
  //____________________________________________________________________
  void DrawTitlePage(TFile* file, Bool_t mc, Bool_t onlyMB)
  {
    TCollection* c   = GetCollection(file, ColName("Forward", true));

    fBody->cd();
    
    Double_t y = .9;
    TLatex* ltx = new TLatex(.5, y, 
			     "#frac{1}{#it{N}}#kern[.1]"
			     "{#frac{d#it{N_{ch}}}{d#it{#eta}}}");
    ltx->SetTextSize(0.07);
    ltx->SetTextFont(42);
    ltx->SetTextAlign(22);
    ltx->SetNDC();
    ltx->Draw();
    y -= .075;

    TObject* tN = GetObject(file, "trainName");
    if (!tN) tN = GetObject(c, "trainName");

    if (tN) {
      TLatex* ltN = new TLatex(.5, y, Form("(%s)", tN->GetTitle()));
      ltN->SetTextSize(0.05);
      ltN->SetTextFont(82);
      ltN->SetTextAlign(22);
      ltN->SetNDC();
      ltN->Draw();
      y -= .055;
    }
    if (mc) {
      ltx = new TLatex(.5, y, "Simulation input");
      ltx->SetNDC();
      ltx->SetTextAlign(23);
      ltx->SetTextFont(42);
      ltx->SetTextSize(.03);
      ltx->Draw();
      y -= .035;
    }
    if (onlyMB) {
      ltx = new TLatex(.5, y, "No centrality");
      ltx->SetNDC();
      ltx->SetTextAlign(23);
      ltx->SetTextFont(42);
      ltx->SetTextSize(.03);
      ltx->Draw();
      y -= .035;
    }
    
    DrawResTitle(c, y, onlyMB);

    PrintCanvas("1/N dN/d#eta");
  }
  //____________________________________________________________________
  void DrawBoth(THStack* rC, THStack* rF, THStack* rM, 
		TFile* file, Bool_t onlyMB)
  {
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

    TH1*  h  = 0;
    if (rM) { 
      TIter nextM(rM->GetHists());
      while ((h = static_cast<TH1*>(nextM()))) rC->Add(h);
    }

    TIter next(rF->GetHists());
    while ((h = static_cast<TH1*>(next()))) rC->Add(h);


    rC->Draw("nostack");

    TCollection* fS = GetCollection(file, ColName("Forward", false));
    UShort_t sys, sNN;
    ULong_t trigger;
    GetParameter(fS, "sNN",     sNN); 
    GetParameter(fS, "sys",     sys); 
    GetParameter(fS, "trigger", trigger); 
    TAxis*   centAxis = GetCentAxis(fS);
    UShort_t cLow     = centAxis && !onlyMB ? centAxis->GetXmin() : 0;
    UShort_t cHigh    = centAxis && !onlyMB ? centAxis->GetXmax() : 100;

    // CompileScript("OtherData.C", "", "RefData", false);

    // If we have V0AND trigger, get NSD other data
    TMultiGraph* other = 0;
#if 0
    if (!centAxis) {
      Int_t   oT = (trigger == 0x2000) ? 0x4 : trigger;
      TString oC = Form("RefData::GetData(%hu,%hu,%hu,%hu,%hu,0xF)", 
			sys, sNN, oT, cLow, cHigh);
      other = reinterpret_cast<TMultiGraph*>(gROOT->ProcessLine(oC));
    }
    else { 
      other = new TMultiGraph("other", "");
      Int_t nCent = centAxis->GetNbins();
      for (Int_t i = 1; i <= nCent; i++) { 
	TString oC = Form("RefData::GetData(%hu,%hu,%hu,%hu,%hu,0xF)", 
			  sys, sNN, 0, UShort_t(centAxis->GetBinLowEdge(i)),
			  UShort_t(centAxis->GetBinUpEdge(i)));
	TMultiGraph* oM = 
	  reinterpret_cast<TMultiGraph*>(gROOT->ProcessLine(oC));
	if (oM) other->Add(oM);
      }
    }
#endif     
    if (other) {
      // p1->Clear();
      // other->Draw("ap");
      // Double_t oMax = other->GetHistogram()->GetMaximum();
      // Double_t rMax = rC->GetMaximum("nostack");
      // other->SetMaximum(1.2*TMath::Max(oMax, rMax));
      // rC->Draw("same nostack");
      TObject* g = 0;
      TIter    nextG(other->GetListOfGraphs());
      while ((g = nextG())) {
	// Printf("Drawing %s/%s", g->GetName(), g->GetTitle());
        g->DrawClone("same p");
      }
    }

    fBody->cd();
    p2->Draw();
    p2->cd();


    TLegend* l = new TLegend(0.01, 0.1, 0.99, 0.99, 
			     onlyMB || !centAxis ? "" : "Centralities");
    l->SetNColumns(fLandscape ? 1 : 2);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    CleanStack(rC, l, onlyMB ? 0 : centAxis);
    TString seen;
    if (other) { 
      TIter nextG(other->GetListOfGraphs());
      TObject* g = 0;
      while ((g = nextG())) {
	if (seen.Index(g->GetTitle()) != kNPOS) continue;
	seen.Append(Form("|%s", g->GetTitle()));
	TLegendEntry* e = l->AddEntry("dummy", g->GetTitle(), "p");
	TGraph* gg = static_cast<TGraph*>(g);
	e->SetMarkerStyle(gg->GetMarkerStyle());
	e->SetMarkerSize(gg->GetMarkerSize());
	e->SetMarkerColor(kBlack);
      }
    }
    l->Draw();

    PrintCanvas("Both");
  }
  //____________________________________________________________________
  void DrawSums(TDirectory* top, const TString& base, bool onlyMB)
  {
    TCollection* c = GetCollection(top, ColName(base));
    if (!c) return;
    
    TAxis* centAxis = (onlyMB ? 0 : GetCentAxis(c));
    if (centAxis && centAxis->GetNbins() < 1) centAxis = 0;

    Int_t    txtPad = 0;
    Double_t xSave  = fParVal->GetX();
    Double_t size   = 0.05;
    fParVal->SetX(.45);
    Double_t y = .8;
    
    if (!onlyMB && centAxis) {
      size = 0.03;
      fBody->Divide(1, 2);
      txtPad = 1;
    
      fBody->cd(1);
      for (Int_t i = 1; i <= centAxis->GetNbins(); i++) { 
	DrawParameter(y, (i == 1 ? "Centrality classes" : ""),
		      Form("%3d%% - %3d%%", 
			   Int_t(centAxis->GetBinLowEdge(i)), 
			   Int_t(centAxis->GetBinUpEdge(i))), size);
      }
    
      TH1* cent = GetH1(c, "cent");
      cent->GetXaxis()->SetRangeUser(0,100);
      cent->SetFillColor(kRed+1);
      cent->SetFillStyle(3002);
      cent->SetXTitle("Centrality [%]");
      cent->SetYTitle("Events");
      cent->SetMaximum(1.3*cent->GetMaximum());
      TH1* centAcc = GetH1(c, "centAcc");
      centAcc->SetFillStyle(3002);

      TLatex* overUnder = new TLatex(0.15, .88,
				     Form("#splitline{<0: %d}{>100: %d}",
					  int(cent->GetBinContent(1)),
					  int(cent->GetBinContent(102))));
      overUnder->SetTextColor(kRed+1);
      overUnder->SetNDC();
      overUnder->SetTextAlign(13);
      overUnder->SetTextFont(42);
      TLatex* overUnderAcc = new TLatex(0.3, .88,
					Form("#splitline{<0: %d}{>100: %d}",
					     int(centAcc->GetBinContent(1)),
					     int(centAcc->GetBinContent(102))));
      overUnderAcc->SetTextColor(kGreen+1);
      overUnderAcc->SetNDC();
      overUnderAcc->SetTextAlign(13);
      overUnderAcc->SetTextFont(42);
      
      DrawInPad(fBody, 2, cent);
      DrawInPad(fBody, 2, centAcc, "same", kLegend);
      DrawInPad(fBody, 2, overUnder, "same");
      DrawInPad(fBody, 2, overUnderAcc, "same");
      
    }
    fBody->cd(txtPad);
    
    UShort_t sys, sNN, scheme;
    ULong_t trigger;
    GetParameter(c, "sNN",     sNN); 
    GetParameter(c, "sys",     sys); 
    GetParameter(c, "scheme",  scheme); 
    GetParameter(c, "trigger", trigger); 

    TString schemeString;
    if (scheme == 0)   schemeString = "1/N_{accepted}";
    if (scheme & 0x1)  schemeString.Append("1/#epsilon_{V}1/#epsilon_{T}");
    if (scheme & 0x2)  schemeString.Append("Shape ");
    if (scheme & 0x4)  schemeString.Append("A+C-E ");
    if (scheme & 0x8)  schemeString.Append("#epsilon_{T,MC} ");
    if (scheme & 0x10) schemeString.Append("0-bin");

    TString trigString;   TriggerString(trigger, trigString);
    TString sysString;    SysString(sys, sysString);
    TString sNNString;    SNNString(sNN, sNNString);
    
    DrawParameter(y, "Collision system",     sysString,   size);
    DrawParameter(y, "#sqrt{s_{NN}}",        sNNString,   size);
    DrawParameter(y, "Normalization scheme", schemeString,size);
    DrawParameter(y, "Triggers",             trigString,  size);
    
    fParVal->SetX(xSave);

    PrintCanvas(Form("%s sums", base.Data()));

    Int_t cLow = centAxis  ?  centAxis->GetXmin() : 0;
    Int_t cHigh = centAxis ? -centAxis->GetXmax() : -100;
    DrawCentSum(c, base, cLow, cHigh);
    if (onlyMB || !centAxis) return;

    for (Int_t i = 1; i <= centAxis->GetNbins(); i++) 
      DrawCentSum(c, base, centAxis->GetBinLowEdge(i), 
		  centAxis->GetBinUpEdge(i));
  }
  //____________________________________________________________________
  void DrawCentSum(const TCollection* sums, const TString& base, 
		   Float_t cLow, Float_t cHigh)
  {
    // Info("DrawCentSum", "Drawing centrality sum [%d,%d] in %s (%s)",
    //      cLow, cHigh, sums->GetName(), base.Data());
    TString title("sums");
    TCollection* c = GetCentCollection(sums, base, cLow, cHigh, title);
    if (!c) return;
    
    TH2* bin  = GetH2(c, "sum");
    TH2* bin0 = GetH2(c, "sum0");
    TH1* type = GetH1(c, "events");
    TH1* trig = GetH1(c, "triggers");
    TH1* stat = GetH1(c, "status");
    if (!bin0 || !bin || !trig || !type) return;

    type->SetFillStyle(3001);
    type->SetFillColor(kGreen+1);

    fBody->Divide(2, 2);

    DrawInPad(fBody, 1, trig, "HIST TEXT");
    DrawInPad(fBody, 2, type, "HIST TEXT");
    DrawInPad(fBody, 3, bin,  "colz");
    
    if (bin0->GetEntries() <= 0) {
      DrawInPad(fBody, 4, stat, "HIST TEXT");
      // fBody->cd(4);
      // TLatex* l = new TLatex(0.5, 0.5, "No 0-bin events");
      // l->SetNDC();
      // l->SetTextAlign(22);
      // l->Draw();
    }
    else
      DrawInPad(fBody, 4, bin0, "colz");
      
    PrintCanvas(title);
  }
  //____________________________________________________________________
  void DrawResTitle(TCollection* c, Double_t& y, Bool_t onlyMB)
  {
    Double_t xSave = fParVal->GetX();
    Double_t size  = 0.05;
    fParVal->SetX(.5);
    // Double_t y = .9;
    TAxis*   centAxis = GetCentAxis(c);
    if (!onlyMB && centAxis) {
      size = 0.03;
      for (Int_t i = 1; i <= centAxis->GetNbins(); i++) { 
	DrawParameter(y, (i == 1 ? "Centrality classes" : ""),
		      Form("%6.2f%% - %6.2f%%", 
			   centAxis->GetBinLowEdge(i), 
			   centAxis->GetBinUpEdge(i)), size);
      }
    }
    TObject* oSNN = GetObject(c, "sNN");
    TString  tSNN; SNNString(oSNN->GetUniqueID(), tSNN);
    TObject* oTrg = GetObject(c,"trigger");
    DrawParameter(y, "Collision system", GetObject(c, "sys")->GetTitle(), size);
    DrawParameter(y, "#sqrt{s_{NN}}",tSNN, size);
    DrawParameter(y, "Trigger",(oTrg ? oTrg->GetTitle() : "?"), size);
    TObject* oscheme = GetObject(c,"scheme");
    TString  scheme  = oscheme ? oscheme->GetTitle() : "";
    if (scheme.IsNull()) scheme = "1/N_{accepted}";
    DrawParameter(y, "Normalization scheme", scheme, size);
    
    Double_t epsT = 0, epsT0 = 0;
    GetParameter(c, "triggerEff",  epsT);
    GetParameter(c, "triggerEff0", epsT0);
    DrawParameter(y, "#epsilon_{T}", Form("%5.3f", epsT), size);
    DrawParameter(y, "#epsilon_{T,zero bin}", Form("%5.3f", epsT0), size);
    Double_t deltaIP =0;
    GetParameter(c, "deltaIP", deltaIP);
    DrawParameter(y, "IP #delta_{xy}", Form("%5.3fmm", deltaIP), size);
    
    TObject*    options = GetObject(c, "options");
    TString     opts(options->GetTitle());
    TObjArray*  tokens = opts.Tokenize(",");
    TObjString* opt = 0;;
    TIter       oNext(tokens);
    Bool_t      first  = true;
    while ((opt = static_cast<TObjString*>(oNext()))) { 
      DrawParameter(y, (first ? "options" : ""), 
		    opt->String().Strip(TString::kBoth), size);
      first = false;
    }

    fParVal->SetX(xSave);      
  }

  //____________________________________________________________________
  THStack* DrawRes(TDirectory* top, const TString& base, Bool_t onlyMB)
  {
    // Info("DrawRes", "Drawing results for %s", base.Data());
    TCollection* c = GetCollection(top, ColName(base, true));
    if (!c) return 0;
    TCollection* s = GetCollection(top, ColName(base, false));
    // if (!s) return 0;

    fBody->cd();
    Double_t y = .9;
    DrawResTitle(c, y, onlyMB);
    PrintCanvas(Form("%s results", base.Data()));

    Int_t     sys = GetObject(c, "sys")->GetUniqueID();
    TObject*  emp = GetObject(c, "empirical");
    // TH1*      emp = GetH1(c, "empirical");
    TF1*      dc  = static_cast<TF1*>(GetObject(c,"deltaCorr"));
    TF1*      vw  = static_cast<TF1*>(GetObject(s,"ipZw"));
    TProfile* sc  = static_cast<TProfile*>(GetObject(s,"sumVsC"));
    if (vw) vw->SetRange(-4,6);
    Int_t     nPad = 0;
    if (emp)  nPad++;
    if (dc)   nPad++;
    if (vw)   nPad++;
    if (sc)   nPad++;
    if (nPad > 0) {
      fBody->Divide(nPad,1);
      Int_t iPad = 1;
      if (emp) DrawInPad(fBody, iPad++, emp, "", 0, "Empirical");
      if (dc)  DrawInPad(fBody, iPad++, dc, "", 0, "\\hbox{IP} \\delta_{xy}");
      if (vw) {
	DrawInPad(fBody, iPad++, vw, "", 0, "\\hbox{IP} \\delta_{z}");
	Double_t y = .95;
	DrawParameter(y, "#mu_{Z}", Form("%5.3f", vw->GetParameter(0)));
	DrawParameter(y, "#sigma_{Z}", Form("%5.3f", vw->GetParameter(1)));
	DrawParameter(y, "#mu_{Z,ref}", Form("%5.3f", vw->GetParameter(2)));
	DrawParameter(y, "#sigma_{Z,ref}", Form("%5.3f", vw->GetParameter(3)));
      }
      if (sc) DrawInPad(fBody, iPad, sc, "", 0, "#LT#Sigma signal#GT");
      PrintCanvas(Form("%s results - corrections", base.Data()));
    }
    
    TAxis*   centAxis = (onlyMB ? 0 : GetCentAxis(c));
    if (centAxis && centAxis->GetNbins() < 1) centAxis = 0;

    THStack* dndeta_  = GetStack(c, "dndeta");
    if (!dndeta_ || !dndeta_->GetHists() ||
	dndeta_->GetHists()->GetEntries() < 0) return 0;
    
    TLegend* l = new TLegend(0.1, 0.1, 0.9, 0.9, 
			     onlyMB || !centAxis? "" : "Centralities");
    l->SetNColumns(fLandscape ? 1 : 2);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    THStack* dndeta   = CleanStack(dndeta_, l, centAxis);

    THStack* dndetaEmp  = GetStack(c, "dndetaEmp");
    if (!dndetaEmp || !dndetaEmp->GetHists() ||
	dndetaEmp->GetHists()->GetEntries() < 0) dndetaEmp = 0;

    THStack* leftRight  = GetStack(c, "leftRight");
    if (!leftRight || !leftRight->GetHists() ||
	leftRight->GetHists()->GetEntries() < 0) leftRight = 0;
    if (leftRight) { leftRight->SetMinimum(0.8); leftRight->SetMaximum(1.2); }
    
    if (!onlyMB) {
      Double_t y1 = fLandscape ? 0  : .3;
      Double_t x2 = fLandscape ? .7 : 1;
      Double_t x1 = fLandscape ? x2 : 0;
      Double_t y2 = fLandscape ? 1  : y1;
      TPad* p1 = new TPad("p1", "p1", 0,  y1, x2, 1,  0, 0);
      TPad* p2 = new TPad("p2", "p2", x1, 0,  1,  y2, 0, 0);
      fBody->cd();
      p1->Draw();
      p1->cd();
      fBody->cd();
      p2->Draw();
      p2->cd();
      p1->Divide(1,3,0,0);

      // fBody->Divide(1, 3, 0, 0);
      
      DrawInPad(p2, 0, l, "");
      DrawInPad(p1, 1, dndeta,  "nostack", 0,
		"\\mathrm{d}N_{\\mathrm{ch}}/\\mathrm{d}\\eta|"
		"_{\\mathrm{incl}}");
      DrawInPad(p1, 2, dndetaEmp, "nostack", (sys == 2 ? kLogy : 0),
		"\\mathrm{d}N_{\\mathrm{ch}}/\\mathrm{d}\\eta|"
		"_{\\mathrm{prim}}");
      DrawInPad(p1, 3, leftRight, "nostack", 0, "Left/Right");
      p1->GetPad(1)->SetGridx();
      p1->GetPad(2)->SetGridx();
      p1->GetPad(3)->SetGridx();
      p1->GetPad(1)->SetGridy();
      p1->GetPad(2)->SetGridy();
      p1->GetPad(3)->SetGridy();
      
      PrintCanvas(Form("%s results - stacks", base.Data()));
    }

    Int_t cLow = centAxis  ?  centAxis->GetXmin() : 0;
    Int_t cHigh = centAxis ? -centAxis->GetXmax() : -100;
    DrawCentRes(c, base, cLow, cHigh);
    if (onlyMB || !centAxis) {
      // Info("", "Returning dndeta for MB");
      dndeta = MakeMBStack(c, base);
      return dndeta;
    }

    for (Int_t i = 1; i <= centAxis->GetNbins(); i++) 
      DrawCentRes(c, base, centAxis->GetBinLowEdge(i), 
		  centAxis->GetBinUpEdge(i));
    
    return dndeta;
  }
  //____________________________________________________________________
  THStack* MakeMBStack(const TCollection* sums, const TString& base)
  {
    TString title("results");
    TCollection* c = GetCentCollection(sums, base, 0, -1, title);
    if (!c) return 0;

    TH1* dndeta      = GetH1(c, Form("dndeta%s",base.Data()));
    if (!dndeta) return 0;

    THStack* ret = new THStack("dndetaMB", title);
    ret->Add(dndeta);

    if (base.EqualTo("MCTruth")) {
      dndeta      = GetH1(c, "dndetaTruth");
      if (dndeta) ret->Add(dndeta);
    }
    return ret;
  }
      
  //____________________________________________________________________
  void DrawCentRes(const TCollection* sums, const TString& base, 
		   Float_t cLow, Float_t cHigh)
  {
    // Info("DrawCentRes", "Drawing centrality results [%d,%d] in %s (%s)",
    //      cLow, cHigh, sums->GetName(), base.Data());
    TString title("results");
    TCollection* c = GetCentCollection(sums, base, cLow, cHigh, title);
    if (!c) return;


    TH1* trig        = GetH1(c, "triggers");
    TH1* norm        = GetH1(c, Form("norm%s",base.Data()));
    TH1* dndeta      = GetH1(c, Form("dndeta%s",base.Data()));
    TH1* dndetaEmp   = GetH1(c, Form("dndeta%sEmp",base.Data()));
    TH2* d2ndetadphi = GetH2(c, Form("d2Ndetadphi%s", base.Data()));
    if (!trig || !norm || !dndeta || !d2ndetadphi) return;
    if (norm->GetEntries() <= 0) return;

    norm->SetFillColor(kGreen+1);
    norm->SetFillStyle(3001);


    fBody->Divide(2, (dndetaEmp ? 4 : 3), 0.05, 0);

    Int_t        trP = 1;
    TVirtualPad* p   = fBody->GetPad(trP);
    p->SetBottomMargin(0.15);
    p->SetLeftMargin(0.15);
    if (trP > 2) p->SetTopMargin(0.05);
    
    DrawInPad(fBody, trP, trig,   "HIST TEXT");
    DrawInPad(fBody, 2,   d2ndetadphi, "colz", 0,
	      "d^{2}#it{N}/d#it{#eta}d#it{phi}|_{incl}");
    DrawInPad(fBody, 4,   norm,   "", 0, "Normalization");
    DrawInPad(fBody, 6,   dndeta, "", 0,
	      "d#it{N}_{ch}/d#it{#eta}|_{incl}");
  
    fBody->GetPad(2)->SetGridx(); 
    fBody->GetPad(4)->SetGridx(); 
    fBody->GetPad(6)->SetGridx(); 
    fBody->GetPad(2)->SetLeftMargin(0.15); 
    fBody->GetPad(4)->SetLeftMargin(0.15); 
    fBody->GetPad(6)->SetLeftMargin(0.15);
    fBody->GetPad(4)->SetRightMargin(0.15);
    fBody->GetPad(6)->SetRightMargin(0.15);
    if (dndetaEmp) {
      DrawInPad(fBody, 8,   dndetaEmp, "", 0,
		"d#it{N}_{ch}/d#it{#eta}|_{prim}");
      fBody->GetPad(8)->SetGridx();
      fBody->GetPad(8)->SetLeftMargin(0.15);
      fBody->GetPad(8)->SetRightMargin(0.15);
    }

    TObject*   normCalc = GetObject(c, "normCalc");
    TString    calc     = normCalc ? normCalc->GetTitle() : "?";

    // Beautify the text
    calc.ReplaceAll("beta", "#beta");
    calc.ReplaceAll("eps", "#varepsilon");
    const char* sufs[] = { "all", "acc", "trg", "vtx", "B", "A", "C", "E", 
			   "V", "T", 0 };
    const char** suf = sufs;
    while (*suf) { 
      calc.ReplaceAll(Form("_%s", *suf), Form("_{%s}", *suf));
      suf++;
    }

    p = fBody->cd(3);
    p->SetPad(p->GetXlowNDC(), 0, 
	      p->GetXlowNDC()+p->GetWNDC(), p->GetYlowNDC()+p->GetHNDC());
    fBody->GetPad(5)->Delete();
    if (dndetaEmp)  fBody->GetPad(7)->Delete();
    TObjArray* lines    = calc.Tokenize("\n");
    // TPaveText* disp     = new TPaveText(.1,.1,.9,.9, "NDC");
    TIter       next(lines);
    TObjString* sline     = 0;
    Double_t y = .95;
    Double_t xSave = fParName->GetX();
    Int_t    aSave = fParName->GetTextAlign();
    Double_t tSave = fParVal->GetTextSize();
    fParName->SetTextAlign(33);
    fParName->SetX(fParVal->GetX()-.05);
    while ((sline = static_cast<TObjString*>(next()))) {
      // disp->AddText(line->GetName());
      TString& line = sline->String();
      Ssiz_t   eq   = line.Last('=');
      if (eq == kNPOS) { 
	DrawParameter(y, line, "", .6*tSave);
	continue;
      }
      TString name = line(0, eq);
      TString val  = line(eq+1,line.Length()-eq-1);
      DrawParameter(y, name.Strip(TString::kBoth), 
		    val.Strip(TString::kBoth),
		    .6*tSave);
      
    }
    fParName->SetTextAlign(aSave);
    fParName->SetX(xSave);
    // disp->SetBorderSize(0);
    // disp->SetBorderSize(0);
    // disp->SetFillStyle(0);
    // DrawInPad(fBody, 3, disp);
    // fBody->cd();

    PrintCanvas(title);

    DrawCentResDetails(c, title);
  }  
  //____________________________________________________________________
  void DrawCentResDetails(const TCollection* sums, const TString& base)
  {
    TString title = TString::Format("%s - details: ", base.Data());

    TCollection* c = GetCollection(sums, "partial");
    if (!c) {
      Warning("", "Collection partical not found in %s", sums->GetName());
      sums->ls();
      return;
    }

    fBody->Divide(3, 1, 0.05, 0);
    
    const char* typs[] = { "", "0", "All" };
    const char* tits[] = { "Non-zero events", "Zero events", "Weighted sum" };
    for (Int_t i = 1; i <= 3; i++) {
      const char* suf = typs[i-1];
      TVirtualPad* p   = fBody->cd(i);
      p->SetTopMargin(0.10);

      TLatex* ltx = new TLatex(0.5, .99, tits[i-1]);
      ltx->SetNDC();
      ltx->SetTextAlign(23);
      ltx->SetTextSize(0.05);
      ltx->Draw();

      TH1* sum  = GetH2(c, Form("sum%s",     suf));
      TH1* norm = GetH1(c, Form("norm%s", suf));
      TH1* phi  = GetH1(c, Form("phi%s",  suf));
    
      norm->SetFillColor(kGreen+1);
      norm->SetFillStyle(3002);
      phi->SetFillColor(kBlue+1);
      phi->SetFillStyle(3001);
      
      p->Divide(1, 3, 0, 0);
      DrawInPad(p, 1, sum, sum->Integral()>0 ? "col" : "", 0, 
		"d^{2}#it{N}_{ch}/d#it{#varphi}d#it{#eta}");
      DrawInPad(p, 2, GetH1(c, Form("average%s", suf)), "", 0, 
		"d#it{N}_{ch}/d#it{#eta}");
      DrawInPad(p, 3, norm, "", 0, "#eta-coverage/#varphi-acceptance");
      DrawInPad(p, 3, phi, "same", kLegend);      
    }
    PrintCanvas(title);
  }

  //____________________________________________________________________
  THStack* CleanStack(const THStack* stack, TLegend* l, const TAxis* axis)
  {
    if (!stack) return 0;
    THStack* ret   = new THStack(stack->GetName(), stack->GetTitle());
    TList*   hists = stack->GetHists();
    TIter    next(hists);
    TH1*     h     = 0;
    Int_t    j     = 0;
    Bool_t   ok    = false;
    while ((h = static_cast<TH1*>(next()))) {
      TString name(h->GetTitle());
      TString nme(h->GetName());
      if (nme.Contains("_mirror", TString::kIgnoreCase)) {
	// Printf("Ignore %s/%s in stack", nme.Data(), name.Data());
	continue;
      }
      if (l && !ok) { 
	j++;
	if (axis) {
	  Int_t bin = j; // axis ? TMath::Min(j, axis->GetNbins()) : 1;
	  if (j >  axis->GetNbins())
	    name = "0% - 100%";
	  else 
	    name.Form("%3d%% - %3d%%", 
		      Int_t(axis->GetBinLowEdge(bin)), 
		      Int_t(axis->GetBinUpEdge(bin)));
	  ok = axis->GetBinUpEdge(bin) > 100;
	}
	else {
	  name.ReplaceAll("ALICE", "");
	  name.ReplaceAll("dNdeta", " - work in progress");
	}
	// Printf("Adding entry %d: %s/%s", j,  nme.Data(), name.Data());
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
