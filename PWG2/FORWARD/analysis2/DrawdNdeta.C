#include <TH1.h>
#include <THStack.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TError.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TImage.h>

/**
 * Class to draw dN/deta results 
 * 
 */
struct dNdetaDrawer 
{
  /**
   * POD of data for range zooming 
   */
  struct RangeParam 
  {
    TAxis*       fMasterAxis; // Master axis 
    TAxis*       fSlave1Axis; // First slave axis 
    TVirtualPad* fSlave1Pad;  // First slave pad 
    TAxis*       fSlave2Axis; // Second slave axis 
    TVirtualPad* fSlave2Pad;  // Second slave pad 
  };
  //__________________________________________________________________
  /** 
   * Constructor 
   * 
   */
  dNdetaDrawer()
    : fShowOthers(false),	// Bool_t 
      fShowAlice(false),        // Bool_t
      fShowRatios(false),	// Bool_t 
      fShowLeftRight(false),	// Bool_t
      fRebin(5),		// UShort_t
      fCutEdges(false),		// Bool_t
      fTitle(""),		// TString
      fHHDFile(""),		// TString
      fTrigString(0),		// TNamed*
      fSNNString(0),		// TNamed*
      fSysString(0),		// TNamed*
      fVtxAxis(0),		// TAxis*
      fForward(0),		// TH1*
      fForwardMC(0),		// TH1*
      fForwardHHD(0),		// TH1*
      fTruth(0),		// TH1*
      fCentral(0),		// TH1*
      fForwardSym(0),		// TH1*
      fForwardMCSym(0),		// TH1*
      fForwardHHDSym(0),	// TH1*
      fTriggers(0),		// TH1*
      fRangeParam(0)

  {
    fRangeParam = new RangeParam;
    fRangeParam->fMasterAxis = 0;
    fRangeParam->fSlave1Axis = 0;
    fRangeParam->fSlave1Pad  = 0;
    fRangeParam->fSlave2Axis = 0;
    fRangeParam->fSlave2Pad  = 0;
  }
  //==================================================================
  /** 
   * @{ 
   * @name Set parameters 
   */
  /** 
   * Show other (UA5, CMS, ...) data 
   * 
   * @param x Whether to show or not 
   */
  void SetShowOthers(Bool_t x)    { fShowOthers = x; }
  //__________________________________________________________________
  /** 
   * Show ALICE published data 
   * 
   * @param x Wheter to show or not 
   */
  void SetShowAlice(Bool_t x)     { fShowAlice = x; }
  //__________________________________________________________________
  /** 
   * Whether to show ratios or not.  If there's nothing to compare to,
   * the ratio panel will be implicitly disabled
   * 
   * @param x Whether to show or not 
   */
  void SetShowRatios(Bool_t x)    { fShowRatios = x; }
  //__________________________________________________________________
  /** 
   * 
   * Whether to show the left/right asymmetry 
   *
   * @param x To show or not 
   */
  void SetShowLeftRight(Bool_t x) { fShowLeftRight = x; }
  //__________________________________________________________________
  /** 
   * Set the rebinning factor 
   * 
   * @param x Rebinning factor (must be a divisor in the number of bins) 
   */
  void SetRebin(UShort_t x)       { fRebin = x; }
  //__________________________________________________________________
  /** 
   * Wheter to cut away the edges 
   * 
   * @param x Whether or not to cut away edges 
   */
  void SetCutEdges(Bool_t x)      { fCutEdges = x; }
  //__________________________________________________________________
  /** 
   * Set the title of the plot
   * 
   * @param x Title
   */
  void SetTitle(TString x)        { fTitle = x; }
  //__________________________________________________________________
  /** 
   * Set the file name of the file containing the HHD results
   * 
   * @param fn File name 
   */
  void SetHHDFile(const char* fn) { fHHDFile = fn; }
  /* @} */
  //==================================================================  
  /** 
   * @{ 
   * @name Override settings from input 
   */
  /** 
   * Override setting from file 
   * 
   * @param sNN Center of mass energy per nucleon pair (GeV)
   */
  void SetSNN(UShort_t sNN) 
  {
    fSNNString = new TNamed("sNN", Form("%04dGeV", sNN));
    fSNNString->SetUniqueID(sNN);
  }
  //__________________________________________________________________
  /** 
   * Set the collision system 
   * - 1: pp 
   * - 2: PbPb
   * 
   * @param sys collision system
   */
  void SetSys(UShort_t sys)
  {
    fSysString = new TNamed("sys", (sys == 1 ? "pp" : 
				    sys == 2 ? "PbPb" : "unknown"));
    fSysString->SetUniqueID(sys);
  }
  //__________________________________________________________________
  /** 
   * Set the vertex range in centimeters 
   * 
   * @param vzMin Min @f$ v_z@f$
   * @param vzMax Max @f$ v_z@f$
   */
  void SetVertexRange(Double_t vzMin, Double_t vzMax) 
  {
    fVtxAxis = new TAxis(10, vzMin, vzMax);
    fVtxAxis->SetName("vtxAxis");
    fVtxAxis->SetTitle(Form("v_{z}#in[%+5.1f,%+5.1f]cm", vzMin, vzMax));
  }
  //__________________________________________________________________
  void SetTrigger(UShort_t trig)
  {
    fTrigString = new TNamed("trigString", (trig & 0x1 ? "INEL" : 
					    trig & 0x2 ? "INEL>0" : 
					    trig & 0x4 ? "NSD" : 
					    "unknown"));
    fTrigString->SetUniqueID(trig);
  }


  //==================================================================
  /** 
   * @{ 
   * @name Main steering functions 
   */
  /** 
   * Run the code to produce the final result. 
   * 
   * @param filename  File containing the data 
   */
  void Run(const char* filename="forward_dndeta.C") 
  {
    if (!Open(filename)) return;

    Double_t max = 0;

    // Create our stack of results
    THStack* results = StackResults(max);

    // Create our stack of other results 
    TMultiGraph* other = 0;
    if (fShowOthers || fShowAlice) other = StackOther(max);
    
    Double_t smax = 0;
    THStack* ratios = 0;
    if (fShowRatios) ratios = StackRatios(other, smax);

    Double_t amax = 0;
    THStack* leftright = 0;
    if (fShowLeftRight) leftright = StackLeftRight(amax);

    Plot(results, other, max, ratios, smax, leftright, amax);
  }
    
  //__________________________________________________________________
  /** 
   * Open input file, and find data 
   * 
   * @param filename File name
   * 
   * @return true on success 
   */
  Bool_t Open(const char* filename)
  {
    TFile* file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("Open", "Cannot open %s", filename);
      return false;
    }
    
    TList* results = static_cast<TList*>(file->Get("ForwardResults"));
    if (!results) { 
      Error("Open", "Couldn't find list ForwardResults");
      return false;
    }

    fForward   = GetResult(results, "dndetaForward");
    fForwardMC = GetResult(results, "dndetaForwardMC");
    fTruth     = GetResult(results, "dndetaTruth");
    if (!fTruth) results->ls();

    TList* clusters = static_cast<TList*>(file->Get("CentralResults"));
    if (!clusters) 
      Warning("Open", "Couldn't find list CentralResults");
    else {
      fCentral   = GetResult(clusters, "dndetaCentral");
      if (fCentral) fCentral->SetMarkerColor(kGreen+1);
    }
    if (!fTrigString) 
      fTrigString = static_cast<TNamed*>(results->FindObject("trigString"));
    if (!fSNNString) 
      fSNNString  = static_cast<TNamed*>(results->FindObject("sNN"));
    if (!fSysString) 
      fSysString  = static_cast<TNamed*>(results->FindObject("sys"));
    if (!fVtxAxis)
      fVtxAxis    = static_cast<TAxis*>(results->FindObject("vtxAxis"));
    
    if (!fTrigString) fTrigString = new TNamed("trigString", "unknown");
    if (!fSNNString)  fSNNString  = new TNamed("sNN", "unknown");
    if (!fSysString)  fSysString  = new TNamed("sys", "unknown");
    if (!fVtxAxis) { 
      fVtxAxis    = new TAxis(1,0,0);
      fVtxAxis->SetName("vtxAxis");
      fVtxAxis->SetTitle("v_{z} range unspecified");
    }

    Info("Open", 
	 "Initialized for\n"
	 "   Trigger:    %s  (%d)\n"
	 "   sqrt(sNN):  %s  (%dGeV)\n"
	 "   System:     %s  (%d)\n"
	 "   Vz range:   %s  (%f,%f)",
	 fTrigString->GetTitle(), fTrigString->GetUniqueID(), 
	 fSNNString->GetTitle(),  fSNNString->GetUniqueID(), 
	 fSysString->GetTitle(),  fSysString->GetUniqueID(), 
	 fVtxAxis->GetTitle(), fVtxAxis->GetXmin(), fVtxAxis->GetXmax());

    TList* sums = static_cast<TList*>(file->Get("ForwardSums"));
    if (sums) 
      fTriggers = GetResult(sums, "triggers");

    if (!fForward) { 
      Error("Open", "Couldn't find the result of the forward analysis");
      return false;
    }
    file->Close();

    
    fForwardHHD = GetHHD();

    return true;
  }
  //__________________________________________________________________
  /** 
   * Make a histogram stack of results 
   * 
   * @param max On return, the maximum value in the stack 
   * 
   * @return Newly allocated stack
   */
  THStack* StackResults(Double_t& max)
  {
    THStack* stack = new THStack("results", "Stack of Results");
    max = TMath::Max(max, AddHistogram(stack, fTruth,      "e5 p"));
    max = TMath::Max(max, AddHistogram(stack, fForwardHHD, "", fForwardHHDSym));
    max = TMath::Max(max, AddHistogram(stack, fForwardMC,  "", fForwardMCSym));
    max = TMath::Max(max, AddHistogram(stack, fCentral,    ""));
    max = TMath::Max(max, AddHistogram(stack, fForward,    "", fForwardSym));
    return stack;
  }
  //__________________________________________________________________
  /** 
   * Make a histogram stack of results 
   * 
   * @param max On return, the maximum value in the stack 
   * 
   * @return Newly allocated stack
   */
  TMultiGraph* StackOther(Double_t& max) const
  {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/OtherData.C");
    Int_t    error = 0;
    Bool_t   onlya = (fShowOthers ? false : true);
    Int_t    trg   = (fTrigString ? fTrigString->GetUniqueID() : 0);
    UShort_t snn   = (fSNNString  ? fSNNString->GetUniqueID() : 0);
    Long_t   ret   = gROOT->ProcessLine(Form("GetData(%d,%d,%d);",
					     snn,trg,onlya));
    if (error) { 
      Error("StackOther", "Failed to execute GetData(%d,%d,%d)", 
	    snn, trg, onlya);
      return 0;
    }
    if (!ret) { 
      Warning("StackOther", "No other data found for sNN=%d, trigger=%d", 
	      snn, trg);
      return 0;
    }
    TMultiGraph* other = reinterpret_cast<TMultiGraph*>(ret);

    TGraphAsymmErrors* o      = 0;
    TIter              next(other->GetListOfGraphs());
    while ((o = static_cast<TGraphAsymmErrors*>(next()))) 
      max = TMath::Max(max, TMath::MaxElement(o->GetN(), o->GetY()));

    return other;
  }
  //__________________________________________________________________
  /** 
   * Make a histogram stack of ratios of results to other data
   * 
   * @param max On return, the maximum value in the stack 
   * 
   * @return Newly allocated stack
   */
  THStack* StackRatios(TMultiGraph* others, Double_t& max) 
  {
    THStack* ratios = new THStack("ratios", "Ratios");

    if (others) {
      TGraphAsymmErrors* ua5_1  = 0;
      TGraphAsymmErrors* ua5_2  = 0;
      TGraphAsymmErrors* alice  = 0;
      TGraphAsymmErrors* cms    = 0;
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(others->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG()))) {
	ratios->Add(Ratio(fForward,          o, max));
	ratios->Add(Ratio(fForwardSym,       o, max));
	ratios->Add(Ratio(fForwardHHD,       o, max));
	ratios->Add(Ratio(fForwardHHDSym,    o, max));
	ratios->Add(Ratio(fCentral,          o, max));
	TString oName(o->GetName());
	oName.ToLower();
	if (oName.Contains("ua5"))  { if (ua5_1) ua5_2 = o; else ua5_1 = o; }
	if (oName.Contains("alice")) alice = o;
	if (oName.Contains("cms"))   cms = o;
      }
      if (ua5_1 && alice) ratios->Add(Ratio(alice, ua5_1, max));
      if (ua5_2 && alice) ratios->Add(Ratio(alice, ua5_2, max));
      if (cms   && alice) ratios->Add(Ratio(alice, cms,   max));
    }

    // Check if we have a primaries from MC 
    if (fTruth) {
      ratios->Add(Ratio(fForward,    fTruth, max));
      ratios->Add(Ratio(fForwardSym, fTruth, max));
      ratios->Add(Ratio(fCentral,    fTruth, max));
    }

    // If we have data from HHD's analysis, then do the ratio of 
    // our result to that data. 
    if (fForwardHHD) { 
      ratios->Add(Ratio(fForward,    fForwardHHD,    max));
      ratios->Add(Ratio(fForwardSym, fForwardHHDSym, max));
    }

    // Do comparison to MC 
    if (fForwardMC) { 
      ratios->Add(Ratio(fForward,    fForwardMC,    max));
      ratios->Add(Ratio(fForwardSym, fForwardMCSym, max));
    }

    // Check if we have ratios 
    if (!ratios->GetHists() || 
	(ratios->GetHists()->GetEntries() <= 0)) { 
      delete ratios; 
      ratios = 0; 
    }
    return ratios;
  }
  //__________________________________________________________________
  /** 
   * Make a histogram stack of the left-right asymmetry 
   * 
   * @param max On return, the maximum value in the stack 
   * 
   * @return Newly allocated stack
   */
  THStack* StackLeftRight(Double_t& max)
  {
    THStack* ret = new THStack("leftright", "Left-right asymmetry");
    ret->Add(Asymmetry(fForward,    max));
    ret->Add(Asymmetry(fForwardHHD, max));
    ret->Add(Asymmetry(fForwardMC,  max));

    if (!ret->GetHists() || 
	(ret->GetHists()->GetEntries() <= 0)) { 
      delete ret; 
      ret = 0; 
    }
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Plot the results
   * 
   * @param results    Results
   * @param others     Other data
   * @param max        Max value 
   * @param ratios     Stack of ratios (optional)
   * @param rmax       Maximum diviation from 1 of ratios 
   * @param leftright  Stack of left-right asymmetry (optional)	 
   * @param amax       Maximum diviation from 1 of asymmetries 
   */
  void Plot(THStack*     results,    
	    TMultiGraph* others, 
	    Double_t     max, 
	    THStack*     ratios,     
	    Double_t     rmax,
	    THStack*     leftright, 
	    Double_t     amax)
  {
    gStyle->SetOptTitle(0);
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetLabelFont(132, "xyz");
    
    Int_t    h = 800;
    Int_t    w = 800; // h / TMath::Sqrt(2);
    if (!ratios) w *= 1.4;
    if (!leftright) w *= 1.4;
    TCanvas* c = new TCanvas("Results", "Results", w, h);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);

    Double_t y1 = 0;
    Double_t y2 = 0;
    Double_t y3 = 0;
    if (ratios)    y1 = 0.3;
    if (leftright) { 
      if (y1 > 0.0001) {
	y2 = 0.2;
	y1 = 0.4;
      }
      else {
	y1 = 0.2;
	y2 = y1;
      }
    }
    PlotResults(results, others, max, y1);
    c->cd();

    PlotRatios(ratios, rmax, y2, y1);
    c->cd( );

    PlotLeftRight(leftright, amax, y3, y2);
    c->cd();

    
    Int_t   vMin = fVtxAxis->GetXmin();
    Int_t   vMax = fVtxAxis->GetXmax();    
    TString trg(fTrigString->GetTitle());
    Int_t   nev  = fTriggers->GetBinContent(fTriggers->GetNbinsX());
    trg          = trg.Strip(TString::kBoth);
    TString base(Form("dndeta_%s_%s_%s_%c%02d%c%02dcm_%09dev",
		      fSysString->GetTitle(), 
		      fSNNString->GetTitle(), 
		      trg.Data(),
		      vMin < 0 ? 'm' : 'p',  TMath::Abs(vMin),
		      vMax < 0 ? 'm' : 'p',  TMath::Abs(vMax),
		      nev));
    c->SaveAs(Form("%s.png",  base.Data()));
    c->SaveAs(Form("%s.root", base.Data()));
    c->SaveAs(Form("%s.C",    base.Data()));
  }
  //__________________________________________________________________
  /** 
   * Plot the results
   *    
   * @param results   Results
   * @param others    Other data
   * @param max       Maximum 
   * @param yd        Bottom position of pad 
   */
  void PlotResults(THStack* results, TMultiGraph* others, 
		   Double_t max, Double_t yd) 
  {
    // Make a sub-pad for the result itself
    TPad* p1 = new TPad("p1", "p1", 0, yd, 1.0, 1.0, 0, 0, 0);
    p1->SetTopMargin(0.05);
    p1->SetBorderSize(0);
    p1->SetBorderMode(0);
    p1->SetBottomMargin(yd > 0.001 ? 0.001 : 0.1);
    p1->SetRightMargin(0.05);
    p1->SetGridx();
    p1->SetTicks(1,1);
    p1->SetNumber(1);
    p1->Draw();
    p1->cd();
    
    results->SetMaximum(1.15*max);
    results->SetMinimum(yd > 0.00001 ? -0.1 : 0);

    FixAxis(results, 1/(1-yd)/1.7, "#frac{1}{N} #frac{dN_{ch}}{d#eta}");

    p1->Clear();
    results->DrawClone("nostack e1");

    fRangeParam->fSlave1Axis = results->GetXaxis();
    fRangeParam->fSlave1Pad  = p1;

    // Draw other data
    if (others) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(others->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG())))
        o->DrawClone("same p");
    }

    // Make a legend in the main result pad
    TLegend* l = p1->BuildLegend(.15,p1->GetBottomMargin()+.01,.90,.35);
    l->SetNColumns(2);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(132);

    // Put a title on top
    TLatex* tit = new TLatex(0.10, 0.95, fTitle.Data());
    tit->SetNDC();
    tit->SetTextFont(132);
    tit->SetTextSize(0.05);
    tit->Draw();

    // Put a nice label in the plot
    TString     eS;
    UShort_t    snn = fSNNString->GetUniqueID();
    const char* sys = fSysString->GetTitle();
    if (snn > 1000) eS = Form("%4.2fTeV", float(snn)/1000);
    else            eS = Form("%3dGeV", snn);
    TLatex* tt = new TLatex(.93, .93, Form("%s #sqrt{s}=%s, %s", 
					   sys, 
					   eS.Data(), 
					   fTrigString->GetTitle()));
    tt->SetNDC();
    tt->SetTextFont(132);
    tt->SetTextAlign(33);
    tt->Draw();

    // Put number of accepted events on the plot
    Int_t nev = fTriggers->GetBinContent(fTriggers->GetNbinsX());
    TLatex* et = new TLatex(.93, .83, Form("%d events", nev));
    et->SetNDC();
    et->SetTextFont(132);
    et->SetTextAlign(33);
    et->Draw();

    // Put number of accepted events on the plot
    if (fVtxAxis) { 
      TLatex* vt = new TLatex(.93, .88, fVtxAxis->GetTitle());
      vt->SetNDC();
      vt->SetTextFont(132);
      vt->SetTextAlign(33);
      vt->Draw();
    }
    // results->Draw("nostack e1 same");

    fRangeParam->fSlave1Axis = FindXAxis(p1, results->GetName());
    fRangeParam->fSlave1Pad  = p1;


    // Mark the plot as preliminary
    TLatex* pt = new TLatex(.12, .93, "Preliminary");
    pt->SetNDC();
    pt->SetTextFont(22);
    pt->SetTextSize(0.07);
    pt->SetTextColor(kRed+1);
    pt->SetTextAlign(13);
    pt->Draw();

    if (!gSystem->AccessPathName("ALICE.png")) { 
      TPad* logo = new TPad("logo", "logo", .12, .65, .25, .85, 0, 0, 0);
      logo->SetFillStyle(0);
      logo->Draw();
      logo->cd();
      TImage* i = TImage::Create();
      i->ReadImage("ALICE.png");
      i->Draw();
    }
    p1->cd();
  }
  //__________________________________________________________________
  /** 
   * Plot the ratios 
   * 
   * @param ratios  Ratios to plot (if any)
   * @param max     Maximum diviation from 1 
   * @param y1      Lower y coordinate of pad
   * @param y2      Upper y coordinate of pad
   */
  void PlotRatios(THStack* ratios, Double_t max, Double_t y1, Double_t y2) 
  {
    if (!ratios) return;
    bool isBottom = (y1 < 0.0001);
    Double_t yd = y2 - y1;
    // Make a sub-pad for the result itself
    TPad* p2 = new TPad("p2", "p2", 0, y1, 1.0, y2, 0, 0, 0);
    p2->SetTopMargin(0.001);
    p2->SetRightMargin(0.05);
    p2->SetBottomMargin(isBottom ? 1/yd * 0.07 : 0.0001);
    p2->SetGridx();
    p2->SetTicks(1,1);
    p2->SetNumber(2);
    p2->Draw();
    p2->cd();

    // Fix up axis
    FixAxis(ratios, 1/yd/1.7, "Ratios", 7);

    ratios->SetMaximum(1+TMath::Max(.22,1.05*max));
    ratios->SetMinimum(1-TMath::Max(.32,1.05*max));
    p2->Clear();
    ratios->DrawClone("nostack e1");

    
    // Make a legend
    TLegend* l2 = p2->BuildLegend(.15,p2->GetBottomMargin()+.01,.9,
				  isBottom ? .6 : .4);
    l2->SetNColumns(2);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);
    l2->SetBorderSize(0);
    l2->SetTextFont(132);

    // Make a nice band from 0.9 to 1.1
    TGraphErrors* band = new TGraphErrors(2);
    band->SetPoint(0, fForwardSym->GetXaxis()->GetXmin(), 1);
    band->SetPoint(1, fForward->GetXaxis()->GetXmax(), 1);
    band->SetPointError(0, 0, .1);
    band->SetPointError(1, 0, .1);
    band->SetFillColor(kYellow+2);
    band->SetFillStyle(3002);
    band->SetLineStyle(2);
    band->SetLineWidth(1);
    band->Draw("3 same");
    band->DrawClone("X L same");
    
    // Replot the ratios on top
    ratios->DrawClone("nostack e1 same");

    if (isBottom) {
      fRangeParam->fMasterAxis = FindXAxis(p2, ratios->GetName());
      p2->AddExec("range", Form("RangeExec((dNdetaDrawer::RangeParam*)%p)", 
				fRangeParam));
    }
    else { 
      fRangeParam->fSlave2Axis = FindXAxis(p2, ratios->GetName());
      fRangeParam->fSlave2Pad  = p2;
    }
  }
  //__________________________________________________________________
  /** 
   * Plot the asymmetries
   * 
   * @param ratios  Asymmetries to plot (if any)
   * @param max     Maximum diviation from 1 
   * @param y1      Lower y coordinate of pad
   * @param y2      Upper y coordinate of pad
   */
  void PlotLeftRight(THStack* leftright, Double_t max, 
		     Double_t y1, Double_t y2) 
  {
    if (!leftright) return;
    bool isBottom = (y1 < 0.0001);
    Double_t yd = y2 - y1;
    // Make a sub-pad for the result itself
    TPad* p3 = new TPad("p3", "p3", 0, y1, 1.0, y2, 0, 0, 0);
    p3->SetTopMargin(0.001);
    p3->SetRightMargin(0.05);
    p3->SetBottomMargin(isBottom ? 1/yd * 0.07 : 0.0001);
    p3->SetGridx();
    p3->SetTicks(1,1);
    p3->SetNumber(2);
    p3->Draw();
    p3->cd();

    TH1* dummy = 0;
    if (leftright->GetHists()->GetEntries() == 1) { 
      // Add dummy histogram
      dummy = new TH1F("dummy","", 10, -6, 6);
      dummy->SetLineColor(0);
      dummy->SetFillColor(0);
      dummy->SetMarkerColor(0);
      leftright->Add(dummy);
    }

    // Fix up axis
    FixAxis(leftright, 1/yd/1.7, "Right/Left", 4);

    leftright->SetMaximum(1+TMath::Max(.12,1.05*max));
    leftright->SetMinimum(1-TMath::Max(.15,1.05*max));
    p3->Clear();
    leftright->DrawClone("nostack e1");

    
    // Make a legend
    TLegend* l2 = p3->BuildLegend(.15,p3->GetBottomMargin()+.01,.9,.5);
    l2->SetNColumns(2);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);
    l2->SetBorderSize(0);
    l2->SetTextFont(132);
#ifndef __CINT__
    if (dummy) {
      TList* prims = l2->GetListOfPrimitives();
      TIter next(prims);
      TLegendEntry* o = 0;
      while ((o = static_cast<TLegendEntry*>(next()))) { 
	TString lbl(o->GetLabel());
	if (lbl != "dummy") continue; 
	prims->Remove(o);
	break;
      }
    }
#endif
    // Make a nice band from 0.9 to 1.1
    TGraphErrors* band = new TGraphErrors(2);
    band->SetPoint(0, fForwardSym->GetXaxis()->GetXmin(), 1);
    band->SetPoint(1, fForward->GetXaxis()->GetXmax(), 1);
    band->SetPointError(0, 0, .05);
    band->SetPointError(1, 0, .05);
    band->SetFillColor(kYellow+2);
    band->SetFillStyle(3002);
    band->SetLineStyle(2);
    band->SetLineWidth(1);
    band->Draw("3 same");
    band->DrawClone("X L same");

    leftright->DrawClone("nostack e1 same");
    if (isBottom) {
      fRangeParam->fMasterAxis = FindXAxis(p3, leftright->GetName());
      p3->AddExec("range", Form("RangeExec((dNdetaDrawer::RangeParam*)%p)", 
				fRangeParam));
    }
    else { 
      fRangeParam->fSlave2Axis = FindXAxis(p3, leftright->GetName());
      fRangeParam->fSlave2Pad  = p3;
    }
  }
  /** @} */
  //==================================================================
  /** 
   * @{ 
   * @name Data utility functions 
   */
  /** 
   * Get a result from the passed list
   * 
   * @param list List to search 
   * @param name Object name to search for 
   * 
   * @return 
   */
  TH1* GetResult(TList* list, const char* name) const 
  {
    if (!list) return 0;
    
    TH1* ret = static_cast<TH1*>(list->FindObject(name));
    if (!ret) 
      Warning("GetResult", "Histogram %s not found", name);
    
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Get the result from previous analysis code 
   * 
   * @param fn  File to open 
   * @param nsd Whether this is NSD
   * 
   * @return null or result of previous analysis code 
   */
  TH1* GetHHD() 
  {
    if (fHHDFile.IsNull()) return 0;
    const char* fn = fHHDFile.Data();
    Bool_t nsd = (fTrigString ? fTrigString->GetUniqueID() & 0x4 : false);
    TDirectory* savdir = gDirectory;
    if (gSystem->AccessPathName(fn)) { 
      Warning("GetHHD", "Output of HHD analysis (%s) not available", fn);
      return 0;
    }
    TFile* file = TFile::Open(fn, "READ");
    if (!file) { 
      Warning("GetHHD", "couldn't open HHD file %s", fn);
      return 0;
    }
    TString hist(Form("dNdeta_dNdeta%s", (nsd ? "NSD" : "")));
    TH1* h = static_cast<TH1*>(file->Get(hist.Data()));
    if (!h) { 
      Warning("GetHHD", "Couldn't find HHD histogram %s in %s", 
	      hist.Data(), fn);
      file->Close();
      savdir->cd();
      return 0;
    }
    TH1* r = static_cast<TH1*>(h->Clone("dndeta_hhd"));
    r->SetTitle("ALICE Forward (HHD)");
    r->SetFillStyle(0);
    r->SetFillColor(0);
    r->SetMarkerStyle(21);
    r->SetMarkerColor(kPink+1);
    r->SetDirectory(0);

    file->Close();
    savdir->cd();
    return r;
  }
  //__________________________________________________________________
  /** 
   * Add a histogram to the stack after possibly rebinning it  
   * 
   * @param stack   Stack to add to 
   * @param hist    histogram
   * @param option  Draw options 
   * 
   * @return Maximum of histogram 
   */
  Double_t AddHistogram(THStack* stack, TH1* hist, Option_t* option) const 
  {
    // Check if we have input 
    if (!hist) return 0;

    // Rebin if needed 
    Rebin(hist);

    stack->Add(hist, option);
    return hist->GetMaximum();
  }
  //__________________________________________________________________
  /** 
   * Add a histogram to the stack after possibly rebinning it  
   * 
   * @param stack   Stack to add to 
   * @param hist    histogram
   * @param option  Draw options 
   * @param sym     On return, the data symmetriced (added to stack)
   * 
   * @return Maximum of histogram 
   */
  Double_t AddHistogram(THStack* stack, TH1* hist, Option_t* option, 
			TH1*& sym) const 
  {
    // Check if we have input 
    if (!hist) return 0;

    // Rebin if needed 
    Rebin(hist);
    stack->Add(hist, option);

    // Now symmetrice the histogram 
    sym = Symmetrice(hist);
    stack->Add(sym, option);

    return hist->GetMaximum();
  }

  //__________________________________________________________________
  /** 
   * Rebin a histogram 
   * 
   * @param h     Histogram to rebin
   * @param rebin Rebinning factor 
   * 
   * @return 
   */
  virtual void Rebin(TH1* h) const
  { 
    if (fRebin <= 1) return;

    Int_t nBins = h->GetNbinsX();
    if(nBins % fRebin != 0) {
      Warning("Rebin", "Rebin factor %d is not a devisor of current number "
	      "of bins %d in the histogram %s", fRebin, nBins, h->GetName());
      return;
    }
    
    // Make a copy 
    TH1* tmp = static_cast<TH1*>(h->Clone("tmp"));
    tmp->Rebin(fRebin);
    tmp->SetDirectory(0);
    tmp->Reset();
    // The new number of bins 
    Int_t nBinsNew = nBins / fRebin;
    for(Int_t i = 1;i<= nBinsNew; i++) {
      Double_t content = 0;
      Double_t sumw    = 0;
      Double_t wsum    = 0;
      Int_t    nbins   = 0;
      for(Int_t j = 1; j<=fRebin;j++) {
	Int_t    bin = (i-1)*fRebin + j;
	Double_t c   =  h->GetBinContent(bin);

	if (c <= 0) continue;

	if (fCutEdges) {
	  if (h->GetBinContent(bin+1)<=0 || 
	      h->GetBinContent(bin-1)) {
	    Warning("Rebin", "removing bin %d=%f of %s (%d=%f,%d=%f)", 
		    bin, c, h->GetName(), 
		    bin+1, h->GetBinContent(bin+1), 
		    bin-1, h->GetBinContent(bin-1));
	    continue;
	  }	
	}
	Double_t e =  h->GetBinError(bin);
	Double_t w =  1 / (e*e); // 1/c/c
	content    += c;
	sumw       += w;
	wsum       += w * c;
	nbins++;
      }
      
      if(content > 0 && nbins > 1 ) {
	tmp->SetBinContent(i, wsum / sumw);
	tmp->SetBinError(i,1./TMath::Sqrt(sumw));
      }
    }

    // Finally, rebin the histogram, and set new content
    h->Rebin(fRebin);
    h->Reset();
    for(Int_t i = 1; i<= nBinsNew; i++) {
      h->SetBinContent(i,tmp->GetBinContent(i));
      h->SetBinError(i,  tmp->GetBinError(i));
    }
    
    delete tmp;
  }
  //__________________________________________________________________
  /** 
   * Make an extension of @a h to make it symmetric about 0 
   * 
   * @param h Histogram to symmertrice 
   * 
   * @return Symmetric extension of @a h 
   */
  TH1* Symmetrice(const TH1* h) const
  {
    Int_t nBins = h->GetNbinsX();
    TH1* s = static_cast<TH1*>(h->Clone(Form("%s_mirror", h->GetName())));
    s->SetTitle(Form("%s (mirrored)", h->GetTitle()));
    s->Reset();
    s->SetBins(nBins, -h->GetXaxis()->GetXmax(), -h->GetXaxis()->GetXmin());
    s->SetMarkerColor(h->GetMarkerColor());
    s->SetMarkerSize(h->GetMarkerSize());
    s->SetMarkerStyle(h->GetMarkerStyle()+4);
    s->SetFillColor(h->GetFillColor());
    s->SetFillStyle(h->GetFillStyle());
    s->SetDirectory(0);

    // Find the first and last bin with data 
    Int_t first = nBins+1;
    Int_t last  = 0;
    for (Int_t i = 1; i <= nBins; i++) { 
      if (h->GetBinContent(i) <= 0) continue; 
      first = TMath::Min(first, i);
      last  = TMath::Max(last,  i);
    }
    
    Double_t xfirst = h->GetBinCenter(first-1);
    Int_t    f1     = h->GetXaxis()->FindBin(-xfirst);
    Int_t    l2     = s->GetXaxis()->FindBin(xfirst);
    for (Int_t i = f1, j=l2; i <= last; i++,j--) { 
      s->SetBinContent(j, h->GetBinContent(i));
      s->SetBinError(j, h->GetBinError(i));
    }
    // Fill in overlap bin 
    s->SetBinContent(l2+1, h->GetBinContent(first));
    s->SetBinError(l2+1, h->GetBinError(first));
    return s;
  }
  //__________________________________________________________________
  /** 
   * Calculate the left-right asymmetry of input histogram 
   * 
   * @param h   Input histogram
   * @param max On return, the maximum distance from 1 of the histogram
   * 
   * @return Asymmetry 
   */
  TH1* Asymmetry(TH1* h, Double_t& max)
  {
    if (!h) return 0;

    TH1* ret = static_cast<TH1*>(h->Clone(Form("%s_leftright", h->GetName())));
    // Int_t    oBins = h->GetNbinsX();
    // Double_t high  = h->GetXaxis()->GetXmax();
    // Double_t low   = h->GetXaxis()->GetXmin();
    // Double_t dBin  = (high - low) / oBins;
    // Int_t    tBins = Int_t(2*high/dBin+.5);
    // ret->SetBins(tBins, -high, high);
    ret->Reset();
    ret->SetTitle(Form("%s (+/-)", h->GetTitle()));
    ret->SetYTitle("Right/Left");
    Int_t nBins = h->GetNbinsX();
    for (Int_t i = 1; i <= nBins; i++) { 
      Double_t x = h->GetBinCenter(i);
      if (x > 0) break;
      
      Double_t c1 = h->GetBinContent(i);
      Double_t e1 = h->GetBinError(i);
      if (c1 <= 0) continue; 
      
      Int_t    j  = h->FindBin(-x);
      if (j <= 0 || j > nBins) continue;

      Double_t c2 = h->GetBinContent(j);
      Double_t e2 = h->GetBinError(j);

      Double_t c12 = c1*c1;
      Double_t e   = TMath::Sqrt((e2*e2*c1*c1+e1*e1*c2*c2)/(c12*c12));
      
      Int_t    k   = ret->FindBin(x);
      ret->SetBinContent(k, c2/c1);
      ret->SetBinError(k, e);
    }
    max = TMath::Max(max, RatioMax(ret));

    return ret;
  }
  //__________________________________________________________________
  /** 
   * Transform a graph into a histogram 
   * 
   * @param g 
   * 
   * @return 
   */
  TH1* Graph2Hist(const TGraphAsymmErrors* g) const
  {
    Int_t    nBins = g->GetN();
    TArrayF  bins(nBins+1);
    Double_t dx = 0;
    for (Int_t i = 0; i < nBins; i++) { 
      Double_t x   = g->GetX()[i];
      Double_t exl = g->GetEXlow()[i];
      Double_t exh = g->GetEXhigh()[i];
      bins.fArray[i]   = x-exl;
      bins.fArray[i+1] = x+exh;
      Double_t dxi = exh+exl;
      if (i == 0) dx  = dxi;
      else if (dxi != dx) dx = 0;
    }
    TString name(g->GetName());
    TString title(g->GetTitle());
    TH1D* h = 0;
    if (dx != 0) {
      h = new TH1D(name.Data(), title.Data(), nBins, bins[0], bins[nBins]);
    }
    else {
      h = new TH1D(name.Data(), title.Data(), nBins, bins.fArray);
    }
    h->SetMarkerStyle(g->GetMarkerStyle());
    h->SetMarkerColor(g->GetMarkerColor());
    h->SetMarkerSize(g->GetMarkerSize());
    
    return h;
  }
  /* @} */
  //==================================================================
  /** 
   * @{ 
   * @name Ratio utility functions 
   */
  /** 
   * Get the maximum diviation from 1 in the passed ratio
   * 
   * @param h Ratio histogram
   * 
   * @return Max diviation from 1 
   */
  Double_t RatioMax(TH1* h) const
  {
    Int_t    nBins = h->GetNbinsX();
    Double_t ret   = 0;
    for (Int_t i = 1; i <= nBins; i++) { 
      Double_t c = h->GetBinContent(i);
      if (c == 0) continue;
      Double_t e = h->GetBinError(i);
      Double_t d = TMath::Abs(1-c-e);
      ret        = TMath::Max(d, ret);
    }
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Compute the ratio of @a h to @a g.  @a g is evaluated at the bin
   * centers of @a h 
   * 
   * @param h  Numerator 
   * @param g  Divisor 
   * 
   * @return h/g 
   */
  TH1* Ratio(const TH1* h, const TGraph* g, Double_t& max) const 
  {
    if (!h || !g) return 0;

    TH1* ret = static_cast<TH1*>(h->Clone("tmp"));
    ret->SetName(Form("%s_over_%s", h->GetName(), g->GetName()));
    ret->SetTitle(Form("%s / %s", h->GetTitle(), g->GetTitle()));
    ret->Reset();
    ret->SetMarkerStyle(g->GetMarkerStyle());
    ret->SetMarkerColor(h->GetMarkerColor());
    ret->SetMarkerSize(0.9*g->GetMarkerSize());
    Double_t xlow  = g->GetX()[0];
    Double_t xhigh = g->GetX()[g->GetN()-1];
    if (xlow > xhigh) { Double_t t = xhigh; xhigh = xlow; xlow = t; }

    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c <= 0) continue;

      Double_t x = h->GetBinCenter(i);
      if (x < xlow || x > xhigh) continue; 

      Double_t f = g->Eval(x);
      if (f <= 0) continue; 

      ret->SetBinContent(i, c / f);
      ret->SetBinError(i, h->GetBinError(i) / f);
    }
    if (ret->GetEntries() <= 0) { 
      delete ret; 
      ret = 0; 
    }
    else 
      max = TMath::Max(RatioMax(ret), max);
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Make the ratio of h1 to h2 
   * 
   * @param h1 First histogram (numerator) 
   * @param h2 Second histogram (denominator)
   * 
   * @return h1 / h2
   */
  TH1* Ratio(const TH1* h1, const TH1* h2, Double_t& max) const
  {
    if (!h1 || !h2) return 0;
    TH1* t1 = static_cast<TH1*>(h1->Clone(Form("%s_%s", 
					       h1->GetName(), 
					       h2->GetName())));
    t1->SetTitle(Form("%s / %s", h1->GetTitle(), h2->GetTitle()));
    t1->Divide(h2);
    t1->SetMarkerColor(h1->GetMarkerColor());
    t1->SetMarkerStyle(h2->GetMarkerStyle());
    max = TMath::Max(RatioMax(t1), max);
    return t1;
  }
  //__________________________________________________________________
  /** 
   * Calculate the ratio of two graphs - g1 / g2
   * 
   * @param g1 Numerator 
   * @param g2 Denominator
   * 
   * @return g1 / g2 in a histogram 
   */
  TH1* Ratio(const TGraphAsymmErrors* g1, 
	     const TGraphAsymmErrors* g2, Double_t& max) const
  {
    Int_t    nBins = g1->GetN();
    TArrayF  bins(nBins+1);
    Double_t dx   = 0;
    for (Int_t i = 0; i < nBins; i++) {
      Double_t x   = g1->GetX()[i];
      Double_t exl = g1->GetEXlow()[i];
      Double_t exh = g1->GetEXhigh()[i];
      bins.fArray[i]   = x-exl;
      bins.fArray[i+1] = x+exh;
      Double_t dxi = exh+exl;
      if (i == 0) dx  = dxi;
      else if (dxi != dx) dx = 0;
    }
    TString name(Form("%s_%s", g1->GetName(), g2->GetName()));
    TString title(Form("%s / %s", g1->GetTitle(), g2->GetTitle()));
    TH1* h = 0;
    if (dx != 0) {
      h = new TH1F(name.Data(), title.Data(), nBins, bins[0], bins[nBins]);
    }
    else {
      h = new TH1F(name.Data(), title.Data(), nBins, bins.fArray);
    }
    h->SetMarkerStyle(g2->GetMarkerStyle());
    h->SetMarkerColor(g1->GetMarkerColor());
    h->SetMarkerSize(0.9*g2->GetMarkerSize());

    Double_t low  = g2->GetX()[0];
    Double_t high = g2->GetX()[g2->GetN()-1];
    if (low > high) { Double_t t = low; low = high; high = t; }
    for (Int_t i = 0; i < nBins; i++) { 
      Double_t x  = g1->GetX()[i];
      if (x < low-dx || x > high+dx) continue;
      Double_t c1 = g1->GetY()[i];
      Double_t e1 = g1->GetErrorY(i);
      Double_t c2 = g2->Eval(x);
      
      h->SetBinContent(i+1, c1 / c2);
      h->SetBinError(i+1, e1 / c2);
    }
    max = TMath::Max(RatioMax(h), max);
    return h;
  }  
  /* @} */
  //==================================================================
  /** 
   * @{ 
   * @name Graphics utility functions 
   */
  /** 
   * Find an X axis in a pad 
   * 
   * @param p     Pad
   * @param name  Histogram to find axis for 
   * 
   * @return Found axis or null
   */
  TAxis* FindXAxis(TVirtualPad* p, const char* name)
  {
    TObject* o = p->GetListOfPrimitives()->FindObject(name);
    if (!o) { 
      Warning("FindXAxis", "%s not found in pad", name);
      return 0;
    }
    THStack* stack = dynamic_cast<THStack*>(o);
    if (!stack) { 
      Warning("FindXAxis", "%s is not a THStack", name);
      return 0;
    }
    if (!stack->GetHistogram()) { 
      Warning("FindXAxis", "%s has no histogram", name);
      return 0;
    }
    TAxis* ret = stack->GetHistogram()->GetXaxis();
    return ret;
  }

  //__________________________________________________________________
  /**
   * Fix the apperance of the axis in a stack
   *
   * @param stack  stack of histogram
   * @param s      Scaling factor
   * @param ytitle Y axis title
   * @param force  Whether to draw the stack first or not
   * @param ynDiv  Divisions on Y axis
   */
  void FixAxis(THStack* stack, Double_t s, const char* ytitle,
               Int_t ynDiv=210, Bool_t force=true)
  {
    if (force) stack->Draw("nostack e1");

    TH1* h = stack->GetHistogram();
    if (!h) return;

    h->SetXTitle("#eta");
    h->SetYTitle(ytitle);
    TAxis* xa = h->GetXaxis();
    TAxis* ya = h->GetYaxis();
    if (xa) {
      xa->SetTitle("#eta");
      // xa->SetTicks("+-");
      xa->SetTitleSize(s*xa->GetTitleSize());
      xa->SetLabelSize(s*xa->GetLabelSize());
      xa->SetTickLength(s*xa->GetTickLength());
    }
    if (ya) {
      ya->SetTitle(ytitle);
      ya->SetDecimals();
      // ya->SetTicks("+-");
      ya->SetNdivisions(ynDiv);
      ya->SetTitleSize(s*ya->GetTitleSize());
      ya->SetTitleOffset(ya->GetTitleOffset()/s);
      ya->SetLabelSize(s*ya->GetLabelSize());
    }
  }
  /* @} */



  //__________________________________________________________________
  Bool_t      fShowOthers;   // Show other data
  Bool_t      fShowAlice;    // Show ALICE published data
  Bool_t      fShowRatios;   // Show ratios 
  Bool_t      fShowLeftRight;// Show asymmetry 
  UShort_t    fRebin;        // Rebinning factor 
  Bool_t      fCutEdges;     // Whether to cut edges
  TString     fTitle;        // Title on plot
  TString     fHHDFile;      // File name of old results
  TNamed*     fTrigString;   // Trigger string (read, or set)
  TNamed*     fSNNString;    // Energy string (read, or set)
  TNamed*     fSysString;    // Collision system string (read or set)
  TAxis*      fVtxAxis;      // Vertex cuts (read or set)
  TH1*        fForward;      // Results
  TH1*        fForwardMC;    // MC results
  TH1*        fForwardHHD;   // Old results
  TH1*        fTruth;        // MC truth
  TH1*        fCentral;      // Central region data
  TH1*        fForwardSym;   // Symmetric extension
  TH1*        fForwardMCSym; // Symmetric extension
  TH1*        fForwardHHDSym;// Symmetric extension
  TH1*        fTriggers;     // Number of triggers
  RangeParam* fRangeParam;   // Parameter object for range zoom 
  
};

//=== Stuff for auto zooming =========================================
void UpdateRange(dNdetaDrawer::RangeParam* p)
{
  if (!p) { 
    Warning("UpdateRange", "No parameters %p", p);
    return;
  }
  if (!p->fMasterAxis) { 
    Warning("UpdateRange", "No master axis %p", p->fMasterAxis);
    return;
  }
  Int_t    first = p->fMasterAxis->GetFirst();
  Int_t    last  = p->fMasterAxis->GetLast();
  Double_t x1    = p->fMasterAxis->GetBinCenter(first);
  Double_t x2    = p->fMasterAxis->GetBinCenter(last);
  //Info("UpdateRange", "Range set to [%3d,%3d]->[%f,%f]", first, last, x1,x2);

  if (p->fSlave1Axis) { 
    Int_t i1 = p->fSlave1Axis->FindBin(x1);
    Int_t i2 = p->fSlave1Axis->FindBin(x2);
    p->fSlave1Axis->SetRange(i1, i2);
    p->fSlave1Pad->Modified();
    p->fSlave1Pad->Update();
  }
  if (p->fSlave2Axis) { 
    Int_t i1 = p->fSlave2Axis->FindBin(x1);
    Int_t i2 = p->fSlave2Axis->FindBin(x2);
    p->fSlave2Axis->SetRange(i1, i2);
    p->fSlave2Pad->Modified();
    p->fSlave2Pad->Update();
  }
  TCanvas*  c = gPad->GetCanvas();
  c->cd();
}
  
//____________________________________________________________________
void RangeExec(dNdetaDrawer::RangeParam* p)
{
  // Event types: 
  //  51:   Mouse move 
  //  53:   
  //  1:    Button down 
  //  21:   Mouse drag
  //  11:   Mouse release 
  // dNdetaDrawer::RangeParam* p = 
  //   reinterpret_cast<dNdetaDrawer::RangeParam*>(addr);
  Int_t event     = gPad->GetEvent();
  TObject *select = gPad->GetSelected();
  if (event == 53) { 
    UpdateRange(p);
    return;
  }
  if (event != 11 || !select || select != p->fMasterAxis) return;
  UpdateRange(p);
}

//=== Steering function ==============================================  
void
DrawdNdeta(const char* filename="forward_dndeta.root", 
	   Int_t       flags=0xf,
	   const char* title="",
	   UShort_t    rebin=5, 
	   Bool_t      cutEdges=false,
	   UShort_t    sNN=0, 
	   UShort_t    sys=0,
	   UShort_t    trg=1,
	   Float_t     vzMin=999, 
	   Float_t     vzMax=-999)
{
  dNdetaDrawer d;
  d.SetRebin(rebin);
  d.SetCutEdges(cutEdges);
  d.SetTitle(title);
  d.SetHHDFile("");
  d.SetShowOthers(flags & 0x1);
  d.SetShowAlice(flags & 0x2);
  d.SetShowRatios(flags & 0x4);
  d.SetShowLeftRight(flags & 0x8);
  // Do the below if your input data does not contain these settings 
  if (sNN > 0) d.SetSNN(sNN);     // Collision energy per nucleon pair (GeV)
  if (sys > 0) d.SetSys(sys);     // Collision system (1:pp, 2:PbPB)
  if (trg > 0) d.SetTrigger(trg); // Collision trigger (1:INEL, 2:INEL>0, 4:NSD)
  if (vzMin < 999 && vzMax > -999) 
    d.SetVertexRange(vzMin,vzMax); // Collision vertex range (cm)
  d.Run(filename);
}
//____________________________________________________________________
//
// EOF
//

