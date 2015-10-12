#ifndef __CINT__
# include "SummaryDrawer.C"
# include "AliFMDCorrAcceptance.h"
# include "AliFMDCorrSecondaryMap.h"
# include "AliFMDCorrELossFit.h"
# include "AliFMDMultCuts.h"
# include "AliForwardUtil.h"
# include "AliForwardCorrectionManager.h"
# include "AliLog.h"
# include <TGraphErrors.h>
# include <TMultiGraph.h>
# include <TString.h>
# include <TError.h>
#else
class SummaryDrawer;
class TObject;
class AliFMDMultCuts;
class THStack;
class TH1;
class TMultiGraph;
class TGraphErrors;
#include <TString.h>
#endif

struct MultCutDrawer : public SummaryDrawer
{
  UShort_t fMinQuality;
  TList    fCuts;
  TList    fStacks;
  TList    fMultiGraphs;
  Bool_t   fMC;
  //__________________________________________________________________
  /** 
   * Constructor
   */
  MultCutDrawer() 
    : fMinQuality(AliFMDCorrELossFit::kDefaultQuality),
      fCuts(),
      fStacks(),
      fMultiGraphs(),
      fMC(false)
  {
    // Rough equvilance: 
    // 
    //   Cut name  |           Parameter values   
    //  -----------+----------------------------------------
    //    mpv      | 0.85     0.7      0.4     0.15
    //    xi       | 1        2.5      4.5     6.8
    //    sig      | .5       1        2       2.9
    //    prob     | 1e-1     2.5e-2   5e-4    2.5e-6
    //  -----------+----------------------------------------
    //   Cut name  |           Mean values   
    //  -----------+----------------------------------------
    //   mpv       | 0.43    0.36      0.20    0.08
    //   xi        | 0.49    0.36      0.23    0.09
    //   sig       | 0.44    0.37      0.22    0.10
    //   prob      | 0.43    0.35      0.21    0.09
    // 
    fCuts.Add(new TNamed("mpv",  "0.85  0.7    0.4  0.15"));
    fCuts.Add(new TNamed("xi",   "1     2.5    4.5  6.8"));
    fCuts.Add(new TNamed("sig",  ".8    .9     1    1.5"));
    fCuts.Add(new TNamed("prob", "0.01  0.025  0.04 0.06"));
    // fCuts.Add(new TNamed("prob", "1e-2 1e-3 1e-5 1e-7"));
  }
  //__________________________________________________________________
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  MultCutDrawer(const MultCutDrawer& o) 
    : fMinQuality(o.fMinQuality),
      fCuts(),
      fStacks(),
      fMultiGraphs(),
      fMC(false)
  {}
  //__________________________________________________________________
  /** 
   * Assignment operator 
   * 
   * @param o Obejct to assign from
   * 
   * @return Reference to this 
   */
  MultCutDrawer& operator=(const MultCutDrawer& o) 
  {
    if (&o == this) return *this;
    fMinQuality = o.fMinQuality;
    fCuts.AddAll(&o.fCuts);
    fStacks.AddAll(&o.fStacks);
    return *this;
  }
  //__________________________________________________________________
  /** 
   * Destructor
   */
  ~MultCutDrawer()
  {
    CloseCanvas();
  }
  //__________________________________________________________________
  /** 
   * Run the class
   * 
   * @param runNo   Run number (or 999 for don't care)
   * @param sys     System (or 0 for don't care)
   * @param sNN     Collision energy in GeV (or 0 for don't care)
   * @param field   L3 Field in kG (or 999 for don't care)
   * @param mc      True of MC
   * @param local   Possible local database 
   */  
  void Run(ULong_t       runNo=999, 
	   UShort_t      sys=0, 
	   UShort_t      sNN=0, 
	   Short_t       field=999, 
	   Bool_t        mc=false, 
	   const Char_t* local=0)
  {
    Bool_t sat = false;
    if (!Init(runNo, sys, sNN, field, mc, sat, local)) return;

    Double_t savX = fParVal->GetX();
    Double_t savY = fParVal->GetY();
    fParVal->SetX(.4);
    fParVal->SetY(.4);
    // fPause = true;

    TIter    iCut(&fCuts);
    TObject* pCut = 0;
    while ((pCut = iCut())) { 
      TString      method(pCut->GetName());
      TString      sP(pCut->GetTitle());
      TObjArray*   aP = sP.Tokenize(" ");
      TIter        iP(aP);
      TObjString*  pP = 0;
      TString      tM;
      TMultiGraph* sum = AllSummary(method, sP);
      fBody->SetBottomMargin(0.20);
      fBody->SetLeftMargin(0.06);
      fBody->Divide(1, aP->GetEntries(), 0, 0);
      Int_t iPad = 1;
      while ((pP = static_cast<TObjString*>(iP()))) {
	THStack*     all   = AllStack(iPad-1);
	Double_t     p     = pP->String().Atof();
	Double_t     vP[]  = { p, p, p, p, p };
	THStack*     stack = CutStack(method, vP, all, sum);
	if (tM.IsNull()) tM = stack->GetTitle();
	// Kill title on all but first sub-panel
	stack->SetTitle("");
	DrawInPad(fBody, iPad, stack, "nostack p");
	stack->GetYaxis()->SetTitleSize(0.12);
	stack->GetYaxis()->SetTitleOffset(0.2);
	stack->GetYaxis()->SetLabelSize(0.07);
	if (iPad == 1) stack->GetYaxis()->SetTitle(tM);
	stack->GetXaxis()->SetTitle("#eta");
	stack->GetXaxis()->SetTitleSize(0.12);
	stack->GetXaxis()->SetTitleOffset(0.6);
	stack->GetXaxis()->SetTitleColor(kBlack);
	stack->GetXaxis()->SetLabelSize(0.07);

	if (iPad == 1) {
	  Color_t    col   = kBlack;
	  Double_t   hLtx  = 0.07;
	  Double_t   yLtx  = 7*(hLtx+.005)+0.01;
	  TLatex*    nLtx  = new TLatex(-0.75, yLtx, "Ring");
	  TLatex*    pLtx  = new TLatex(-0.7,  yLtx, "Param");
	  TLatex*    vLtx  = new TLatex(+0, yLtx, "Mean#pmVar_{min}^{max}");
	  nLtx->SetTextAlign(31);pLtx->SetTextAlign(11);
	  nLtx->SetTextSize(hLtx);
	  pLtx->SetTextSize(hLtx);
	  vLtx->SetTextSize(hLtx);
	  nLtx->SetTextColor(col);
	  pLtx->SetTextColor(col);
	  vLtx->SetTextColor(col);
	  nLtx->Draw();
	  pLtx->Draw();
	  vLtx->Draw();
	}
	// if (iPad == 1) { 
	//   fBody->cd(1);
	//   DrawRingLegend(0.4, 0.4, 0.7, 0.9);
	// }
	iPad++;
      }
      PrintCanvas(Form("%s   X={%s}", tM.Data(), sP.Data()));
    }

    Int_t nAll = fStacks.GetEntries();
    fBody->SetBottomMargin(0.20);
    fBody->SetLeftMargin(0.06);
    fBody->Divide(1, nAll, 0, 0);
    for (Int_t iAll = 0; iAll < nAll; iAll++) {
      THStack* all = AllStack(iAll);
      DrawInPad(fBody, iAll+1, all, "nostack hist p");
      all->GetYaxis()->SetTitleSize(0.12);
      all->GetYaxis()->SetTitleOffset(0.2);
      all->GetYaxis()->SetLabelSize(0.07);
      if (iAll == 0) all->GetYaxis()->SetTitle("c");
      all->GetXaxis()->SetTitle("#eta");
      all->GetXaxis()->SetTitleSize(0.12);
      all->GetXaxis()->SetTitleOffset(0.6);
      all->GetXaxis()->SetTitleColor(kBlack);
      all->GetXaxis()->SetLabelSize(0.07);
      
      TVirtualPad* p = fBody->GetPad(iAll+1);
      p->cd();
      Double_t yT = 1-p->GetTopMargin();
      if      (iAll == 0) DrawRingLegend(p, kNorth|kCenter); 
      else if (iAll == 1) DrawMethodLegend(0.35, 0.4, 0.55,yT);
      
      Double_t y1 = ((iAll + 2 >= nAll) ? yT - .3 : p->GetBottomMargin());
      Double_t y2 = ((iAll + 2 >= nAll) ? yT      : 0.3);
      DrawValueLegend(all, 0.2, y1, 0.9, y2);
       
    }
    PrintCanvas("Comparisons");

    fParVal->SetX(savX);
    fParVal->SetY(savY);

    Int_t nSum = fMultiGraphs.GetEntries();
    fBody->Divide(1,nSum);
    for (Int_t i = 1; i <= nSum; i++) {
      DrawInPad(fBody, i, fMultiGraphs.At(i-1), "apl");
    }
    PrintCanvas("Trends");
    
    
    TFile* out = TFile::Open("cutMethods.root", "RECREATE");
    fStacks.Write("stacks", TObject::kSingleKey);
    fMultiGraphs.Write("trends", TObject::kSingleKey);
    TIter nextMG(&fMultiGraphs);
    TMultiGraph* mg = 0;
    while ((mg = static_cast<TMultiGraph*>(nextMG()))) {
      TDirectory* dir = out->mkdir(mg->GetName());
      dir->cd();
      mg->GetListOfGraphs()->Write();
      out->cd();
    }
    out->Write();
    // out->Close();
    
    CloseCanvas();
  }
  //__________________________________________________________________
  /** 
   * Draw legend of methods 
   * 
   * @param x1    Left x
   * @param y1    Lower y
   * @param x2    Right x
   * @param y2    Upper y
   */
  void DrawMethodLegend(Double_t x1, Double_t y1, 
			Double_t x2, Double_t y2) 
  {
    TLegend* l = new TLegend(x1,y1,x2,y2);
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetFillStyle(0);

    TIter    iCut(&fCuts);
    TObject* pCut = 0;
    while ((pCut = iCut())) { 
      TString                 method(pCut->GetName());
      AliFMDMultCuts::EMethod m = AliFMDMultCuts::String2Method(method);
      TString                 title(AliFMDMultCuts::Method2String(m,true));
      Style_t                 style = CutStyle(m);
      TLegendEntry*           e = l->AddEntry("dummy", title, "p");
      e->SetMarkerStyle(style);
      e->SetMarkerColor(kBlack);
    }
    l->Draw();
  }
  //__________________________________________________________________
  /** 
   * Draw a value legend 
   * 
   * @param stack Stack to take values from 
   * @param x1    Left x
   * @param y1    Lower y
   * @param x2    Right x
   * @param y2    Upper y
   */
  void DrawValueLegend(THStack* stack, 
		       Double_t x1, Double_t y1, 
		       Double_t x2, Double_t y2) 
  {
    TLegend* l = new TLegend(x1,y1,x2,y2);
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetFillStyle(0);

    TString seen;
    TIter   iHist(stack->GetHists());
    TH1*    pHist = 0;
    Int_t   nHist = 0;
    while ((pHist = static_cast<TH1*>(iHist()))) {
      TString name(pHist->GetName());
      if (seen.Contains(name)) continue; 
      seen.Append(Form(" %s", name.Data()));
      nHist++;
      
      AliFMDMultCuts::EMethod m = AliFMDMultCuts::String2Method(name);
      Style_t                 s = CutStyle(m);
      TLegendEntry*           e = l->AddEntry("dummy", pHist->GetTitle(), "p");
      e->SetMarkerStyle(s);
      e->SetMarkerColor(kBlack);
    }
    if (nHist < 5) l->SetNColumns(nHist);
    else           l->SetNColumns(nHist/2);

    l->Draw();
  }
  //__________________________________________________________________
  /** 
   * Get stack at @a i.  If the stack doesn't exist, make it
   * 
   * @param i Location (0-based)
   * 
   * @return Stack 
   */
  THStack* AllStack(Int_t i) 
  {
    TObject* o = fStacks.At(i);
    if (o) return static_cast<THStack*>(o);
    THStack* s = new THStack(Form("all%02d", i), "");
    fStacks.AddAt(s, i);
    return s;
  }
  //__________________________________________________________________
  /** 
   * Get multigraph at @a i.  If the stack doesn't exist, make it
   * 
   * @param method Method used 
   * @param title  Title on plot
   * 
   * @return Stack 
   */
  TMultiGraph* AllSummary(const TString& method, const TString& title)
  {
    TObject* o = fMultiGraphs.FindObject(method);
    if (o) return static_cast<TMultiGraph*>(o);
    TMultiGraph* mg = new TMultiGraph(method, title);
    fMultiGraphs.Add(mg);
    return mg;
  }
  /** 
   * Find summary 
   * 
   * @param summaries List of summaries
   * @param n         Which name 
   * @param method    Which method 
   * @param col       The color to use 
   * @param style     The style to use 
   * 
   * @return Found graph or null
   */
  TGraphErrors* FindSummary(TMultiGraph* summaries,
			    const TString& n,
			    const TString& method,
			    Color_t        col,
			    Style_t        style)
  {
    
    TObject* o = (summaries->GetListOfGraphs() ?
		  summaries->GetListOfGraphs()->FindObject(n) : 0);
    if (o) return static_cast<TGraphErrors*>(o);
    TGraphErrors* summary = new TGraphErrors;
    summary->SetName(n);
    summary->SetTitle(method);
    summary->SetLineColor(col);
    summary->SetMarkerColor(col);
    summary->SetFillColor(col);
    summary->SetMarkerStyle(style);
    summary->SetFillStyle(0);
    summaries->Add(summary);
    
    return summary;
  }

  //__________________________________________________________________
  /** 
   * Get the marker styoe associated with a cut
   * 
   * @param m Cut identifier 
   * 
   * @return Marker style 
   */
  static Style_t CutStyle(UShort_t m) 
  {
    switch (m) { 
    case AliFMDMultCuts::kFixed:            return kFullStar;
    case AliFMDMultCuts::kMPVFraction:      return kOpenCircle;
    case AliFMDMultCuts::kFitRange:         return 33; // Diamond
    case AliFMDMultCuts::kLandauWidth:      return 34; // Cross
    case AliFMDMultCuts::kLandauSigmaWidth: return kOpenSquare;
    case AliFMDMultCuts::kProbability:      return kFullTriangleDown;
    }
    return kFullDotMedium;
  }
  //__________________________________________________________________
  /** 
   * Update statistics
   * 
   * @param y     Current value 
   * @param cnt   Current count
   * @param mean  Current mean
   * @param var   Current variance 
   */
  static void Statistics(Double_t  y,
			 Int_t&    cnt,
			 Double_t& mean, 
			 Double_t& var) 
  {
    cnt        += 1;
    mean       += (y - mean) / cnt;
    var        += (cnt > 1 ? (TMath::Power(y-mean,2)/(cnt-1)-var/cnt) : 0);
  }
  //__________________________________________________________________
  /** 
   * Calculate statistics for a histogram
   * 
   * @param h     Histogram
   * @param mean  On return, the mean of y
   * @param var   On return, the variance in y
   * @param min   On return, the least y
   * @param max   On return, the largest y
   * @param rCnt  In/out: Current count
   * @param rMean In/out: Current mean 
   * @param rVar  In/out: Current variance 
   */
  static void HistStatistics(const TH1* h, 
			     Double_t& mean, 
			     Double_t& var, 
			     Double_t& min, 
			     Double_t& max,
			     Int_t&    rCnt,
			     Double_t& rMean,
			     Double_t& rVar)
  {
    mean      = 0;
    var       = 0;
    min       = +100000;
    max       = -100000;
    Int_t cnt = 0;
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t y = h->GetBinContent(i);
      if (TMath::Abs(y) <= 1e-9) continue;
      min        =  TMath::Min(min, y);
      max        =  TMath::Max(max, y);
      Statistics(y, cnt, mean, var);
      Statistics(y, rCnt, rMean, rVar);
    }
    // Info("", "Stats for %s:  mean=%f +/- %f [%f,%f]",
    //      h->GetTitle(), mean, var, min, max);
  }
  //__________________________________________________________________
  /** 
   * Create a stack from cuts
   * 
   * @param method Method to use 
   * @param param  Parameters 
   * @param all    Stack to add for this set of parameters
   * @param summaries List of summaries
   * 
   * @return Newly created stack 
   */    
  THStack* CutStack(const TString& method, Double_t* param,
		    THStack* all, TMultiGraph* summaries)
  {
    AliFMDMultCuts::EMethod m = AliFMDMultCuts::String2Method(method);
    Info("CutStack", "Method %s -> %d", method.Data(), m);
    AliFMDMultCuts* cut = new AliFMDMultCuts(m, 
					     param[0],
					     param[1], 
					     param[2],
					     param[3],
					     param[4]);
    // cut->Print();
    
    TH2* hist = new TH2D("cut", cut->GetMethodString(true),
			 200, -4, 6, 5, 0.5, 5.5);
    hist->GetYaxis()->SetBinLabel(1, "FMD1i");
    hist->GetYaxis()->SetBinLabel(2, "FMD2i");
    hist->GetYaxis()->SetBinLabel(3, "FMD2o");
    hist->GetYaxis()->SetBinLabel(4, "FMD3o");
    hist->GetYaxis()->SetBinLabel(5, "FMD3i");
    
    // Info("DrawMultCuts", "Filling histogram");
    cut->FillHistogram(hist);
    // Info("DrawMultCuts", "Done filling");

    Style_t  style = CutStyle(m);
    THStack* stack = new THStack(hist, "x");
    TList*   hists = stack->GetHists();
    Double_t rMin  = +1000000;
    Double_t rMax  = -1000000;
    Double_t rAvg  = 0;
    Double_t rVar  = 0;
    Int_t    rCnt  = 0;
    printf(" %6s %7.4f | ", method.Data(), param[0]);
    TH1*     first = 0;
    for (Int_t i = 1; i <= 5; i++) { 
      TH1*     h = static_cast<TH1*>(hists->At(i-1));
      TString  n(hist->GetYaxis()->GetBinLabel(i)); 
      TString  nn(n); nn.Remove(0,3);
      UShort_t det = nn.Atoi();
      Char_t   rng = nn[1];
      Color_t  col = AliForwardUtil::RingColor(det, rng);
      if (!first) first = h;
      h->SetName(method);
      h->SetTitle(Form("%f", param[i-1]));
      h->SetYTitle(cut->GetMethodString(true));
      h->SetXTitle("#eta");
      h->SetMarkerColor(col);
      h->SetFillColor(col);
      h->SetLineColor(col);
      h->SetMarkerStyle(style);
      h->SetFillStyle(0);
      Double_t avg, var, min, max;
      
      HistStatistics(h, avg, var, min, max, rCnt, rAvg, rVar);
      rMin = TMath::Min(min, rMin);
      rMax = TMath::Max(max, rMax);
      all->Add(h);
      Double_t   hLtx  = 0.07;
      Double_t   yLtx  = i*(hLtx+.005)+0.01;
      TObjArray* lines = new TObjArray(3);
      TLatex*    nLtx  = new TLatex(-0.75, yLtx, n);
      TLatex*    pLtx  = new TLatex(-0.7,  yLtx, Form("X=%g", param[i-1]));
      TLatex*    vLtx  = new TLatex(+0, yLtx, 
				    Form("%5.3f#pm%6.4f_{%6.4f}^{%6.4f}",
					 avg, var, max-avg, avg-min));
      nLtx->SetTextAlign(31);pLtx->SetTextAlign(11);
      nLtx->SetTextSize(hLtx);pLtx->SetTextSize(hLtx),vLtx->SetTextSize(hLtx);
      nLtx->SetTextColor(col);pLtx->SetTextColor(col);vLtx->SetTextColor(col);
      lines->Add(nLtx);lines->Add(pLtx);lines->Add(vLtx);
      h->GetListOfFunctions()->Add(lines);
      printf("%5.3f+/-%6.4f ", avg, var);

      TGraphErrors* summary = FindSummary(summaries, n, method, col, style);
      Int_t nSum = summary->GetN();
      summary->SetPoint(nSum, param[i-1], avg);
      summary->SetPointError(nSum, 0, var);
    }
    TLatex* rLtx = new TLatex(6, fMC ? 0.65 : 0.55, 
			      Form("All: %5.3f#pm%6.4f_{%6.4f}^{%6.4f}",
				   rAvg, rVar, rMin, rMax));
    rLtx->SetTextSize(0.05);
    rLtx->SetTextAlign(31);
    first->GetListOfFunctions()->Add(rLtx);
    Printf("-> %5.3f+/-%6.4f", rAvg, rVar);
    stack->SetTitle(cut->GetMethodString(true)); // hist->GetTitle());
    stack->SetMinimum(0); // 0.98*min);
    stack->SetMaximum(fMC ? 0.7 : 0.6); // 1.02*max);
    all->SetMinimum(0);
    all->SetMaximum(fMC ? 0.7 : 0.6);

    summaries->SetTitle(cut->GetMethodString(true));
    TGraphErrors* summary = FindSummary(summaries, "all", method,kBlack,style);
    Int_t nSum = summary->GetN();
    summary->SetPoint(nSum, param[0], rAvg);
    summary->SetPointError(nSum, 0, rVar);
    
    delete hist;
    return stack;
  }
  //__________________________________________________________________
  /** 
   * Initialize 
   * 
   * @param runNo   Run number (or 999 for don't care)
   * @param sys     System (or 0 for don't care)
   * @param sNN     Collision energy in GeV (or 0 for don't care)
   * @param field   L3 Field in kG (or 999 for don't care)
   * @param mc      True of MC
   * @param sat     True for including satellite collisions 
   * @param local   Possible local database 
   * 
   * @return true on sucess 
   */
  Bool_t Init(ULong_t       runNo=999, 
	      UShort_t      sys=0, 
	      UShort_t      sNN=0, 
	      Short_t       field=999, 
	      Bool_t        mc=false, 
	      Bool_t        sat=false,
	      const Char_t* local=0)
  {
    fMC = mc;
    AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
    mgr.SetDebug(true);
    UShort_t flags = AliForwardCorrectionManager::kELossFits;
  
    if (local && local[0] != '\0') mgr.SetELossFitsPath(local);
  
    if (!mgr.Init(runNo, sys, sNN, field, mc, false, flags, true)) {
      Error("DrawMultCuts", "Failed to initialize for flags=0x%02x, "
	    "run=%lu, sys=%hu, sNN=%hu, field=%hd, mc=%s, sat=%s",
	    flags, runNo, sys, sNN, field, mc ? "true" : "false", "false");
      return false;
    }
    const AliFMDCorrELossFit* cFit = mgr.GetELossFit();
    AliFMDCorrELossFit*       fit  = const_cast<AliFMDCorrELossFit*>(cFit);
    fit->CacheBins(8);

    CreateCanvas("multCuts.pdf", true);

    fBody->cd();
    
    Double_t y = .85;
    TLatex* title = new TLatex(.5, y, "#Delta Cuts");
    title->SetTextAlign(23);
    title->SetTextFont(42);
    title->SetTextSize(0.1);
    title->Draw();
    
    y -= 0.11;
    DrawParameter(y, "Run #", Form("%lu", runNo));
    DrawParameter(y, "System", AliForwardUtil::CollisionSystemString(sys));
    DrawParameter(y, "#sqrt{s_{NN}}", 
		  AliForwardUtil::CenterOfMassEnergyString(sNN));
    DrawParameter(y, "L3 field", AliForwardUtil::MagneticFieldString(field));
    DrawParameter(y, "Simulation", Form("%s", mc ? "yes" : "no"));
    DrawParameter(y, "Satellite", Form("%s", sat ? "yes" : "no"));
    PrintCanvas("Delta cuts");

    return true;
  }
};
