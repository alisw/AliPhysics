/**
 * @file   QAPlotter.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 12:08:04 2011
 * 
 * @brief  Class to plot QA trends
 * 
 * @ingroup pwglf_forward_qa_scripts
 * 
 */
#ifndef __CINT__
# include "QABase.h"
# include <TTree.h>
# include <TGraph.h>
# include <TGraphErrors.h>
# include <TGraphAsymmErrors.h>
# include <TMultiGraph.h>
# include <TH1.h>
# include <TLegend.h>
# include <TCanvas.h>
# include <TLine.h>
# include <TArrow.h>
# include <TArrayI.h>
# include <TMath.h>
# include <TChain.h>
#else 
class QABase;
class QARing;
class TTree;
class TGraphAsymmErrors;
class TGraph;
class TGraphErrors;
class TMultiGraph;
class RingQuantity;
class TArrayI;
class TH1;
#endif


/**
 * Class to plot QA trends 
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
struct QAPlotter : public QABase
{
  /**
   * Ring class 
   */
  struct Ring : public QARing
  {
    /** 
     * Constuctor
     * 
     * @param d Detector 
     * @param r Ring 
     */
    Ring(UShort_t d, Char_t r, Bool_t useVar=false)
      : QARing(d, r),
        fGChi2(0),
        fGC(0),
        fGDelta(0),
        fGXi(0),
        fGSigma(0),
        fGLow(0),
        fGSingles(0),
        fGLoss(0),
        fGBeta(0),
        fGOccupancy(0),
	fUseVar(useVar)
    {
      fGChi2           = new TGraphAsymmErrors;
      fGC              = new TGraphAsymmErrors;
      fGDelta          = new TGraphAsymmErrors;
      fGXi             = new TGraphAsymmErrors;
      fGSigma          = new TGraphAsymmErrors;
      fGLow            = new TGraph;
      fGSingles        = new TGraph;
      fGLoss           = new TGraph;
      fGBeta           = new TGraph;
      fGOccupancy      = new TGraphAsymmErrors;

      SetAtt(fGChi2,		"chi2",     "#LT#chi^{2}/#nu#GT");
      SetAtt(fGC,		"c",        "Constant");
      SetAtt(fGDelta,		"delta",    "#LT#Delta_{p}#GT");
      SetAtt(fGXi,		"xi",       "#LT#xi#GT");
      SetAtt(fGSigma,		"sigma",    "#LT#sigma#GT");
      SetAtt(fGLow,		"low",      "# of low statistics bins");
      SetAtt(fGSingles,		"singles",  "Fraction of single hits");
      SetAtt(fGLoss,            "loss",     "Data lossed due to cuts");
      SetAtt(fGOccupancy,	"occupancy","#LTOccupancy#GT");
      SetAtt(fGBeta,		"beta",     "Correlation of methods");

    }
    /** 
     * Static member function to get the ring color
     * 
     * @param d Detector number
     * @param r Ring identifer 
     * 
     * @return Color
     */
    static Color_t RingColor(UShort_t d, Char_t r)
    { 
      return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	      + ((r == 'I' || r == 'i') ? 2 : -3));
    }
    /** 
     * Set graph attributes
     * 
     * @param g      Graph
     * @param name   Name of graph
     */
    void SetAtt(TGraph* g, const char* name, const char* /*title=""*/) 
    {
      Color_t c = RingColor(fD, fR);
      g->SetName(Form("FMD%d%c_%s", fD, fR, name));
      // g->SetTitle(Form("FMD%d%c %s", fD, fR, 
      // !title || title[0] == '\0' ? name : title));
      Int_t marker = 20+(fD-1) + (fR == 'I' ? 0 : 4);
      g->SetTitle(Form("FMD%d%c", fD, fR));
      g->SetLineColor(c);
      g->SetFillColor(c);
      g->SetMarkerColor(c);
      g->SetMarkerStyle(marker);
      g->SetLineWidth(2);
      switch (marker) { 
      case 20: g->SetMarkerSize(1.2); break;
      case 21: g->SetMarkerSize(1.2); break;
      case 22: g->SetMarkerSize(1.3); break;
      case 26: g->SetMarkerSize(1.1); break;
      }
    }
    /** 
     * Update a graph from a RingQuantity 
     * 
     * @param g      Graph to update
     * @param q      Quantity 
     * @param runNo  Run number 
     */
    void UpdateGraph(TGraphAsymmErrors* g, RingQuantity& q, 
		     UInt_t runNo, UInt_t /* n */) 
    {
      Double_t y  = q.mean;
      Double_t el = y-q.min;
      Double_t eh = q.max-y;
      if (fUseVar) { 
	Info("UpdateGraph", "Setting errors on %s to variance",
	     g->GetName());
	el = q.var; 
	eh = q.var; 
      }
      if (TMath::Abs(y) < 1e-6) return;
      Int_t    i  = g->GetN();
      g->SetPoint(i, runNo, y);
      g->SetPointError(i, 0, 0, el, eh);
    }
    /** 
     * Update all graphs
     * 
     * @param n     Entry number (not used)
     * @param runNo Run number 
     * 
     * @return true on success
     */
    Bool_t Update(UInt_t n, UInt_t runNo)
    {
      UpdateGraph(fGChi2,      *fChi2, 		runNo, n);
      UpdateGraph(fGC,         *fC, 		runNo, n);
      UpdateGraph(fGDelta,     *fDelta, 	runNo, n);
      UpdateGraph(fGXi,        *fXi, 		runNo, n);
      UpdateGraph(fGSigma,     *fSigma, 	runNo, n);
      UpdateGraph(fGOccupancy, *fOccupancy, 	runNo, n);

      fGLow->SetPoint(n, runNo, fFitStatus->nLow);
      
      if (fMerge->one > 1e-6)
	fGSingles->SetPoint(fGSingles->GetN(), runNo, fMerge->one);
      if (fCorrelation->beta > 1e-6)
	fGBeta->SetPoint(fGBeta->GetN(), runNo, fCorrelation->beta);
      if (-fDataLoss->full > 1e-6)
	fGLoss->SetPoint(fGLoss->GetN(), runNo, -fDataLoss->full);
      return true;
    }
    TGraphAsymmErrors* fGChi2;	   // Graph of ELoss reduced chi-square
    TGraphAsymmErrors* fGC;	   // Graph of ELoss constant 
    TGraphAsymmErrors* fGDelta;    // Graph of ELoss MPV
    TGraphAsymmErrors* fGXi;       // Graph of ELoss Landau width
    TGraphAsymmErrors* fGSigma;    // Graph of ELoss Gaus width
    TGraph*            fGLow;      // Graph of bins with low statistics
    TGraph*            fGSingles;  // Graph of fraction of singles
    TGraph*            fGLoss;     // Graph of 'lost' data 
    TGraph*            fGBeta;     // Graph of Poisson vs ELoss correlation
    TGraphAsymmErrors* fGOccupancy;// Graph of mean occupancy              
    Bool_t             fUseVar;    // Use variance 
  };
  /** 
   * Constructor 
   */
  QAPlotter(Int_t prodYear=0, Char_t prodLetter='\0', Bool_t useVar=false) 
    : QABase(false, prodYear, prodLetter),
      fNAccepted(0),
      fVz(0),
      fUseVar(useVar)
  {
    Info("QAPlotter", "Do we use variance? %s", fUseVar ? "yes" : "no");
    fFMD1i = new Ring(1, 'I', useVar); 
    fFMD2i = new Ring(2, 'I', useVar); 
    fFMD2o = new Ring(2, 'O', useVar); 
    fFMD3i = new Ring(3, 'I', useVar); 
    fFMD3o = new Ring(3, 'O', useVar); 
    fNAccepted = new TGraph;
    fNAccepted->SetName("nAccepted");
    fNAccepted->SetMarkerStyle(20);
    fNAccepted->SetLineWidth(2);

    fVz = new TGraphErrors;
    fVz->SetName("vz");
    fVz->SetMarkerStyle(20);
    fVz->SetLineWidth(2);
  }
  /** 
   * Add a file to be processed
   * 
   * @param filename Name of file 
   */
  void AddFile(const char* filename)
  {
    fFiles.Add(new TObjString(filename));
  }
  /** 
   * Make a tree 
   * 
   * @param read If true, read from file 
   * 
   * @return True on success 
   */
  Bool_t MakeTree(bool read)
  {
    if (fFiles.GetEntriesFast() <= 0) return QABase::MakeTree(read);

    TChain* chain = new TChain("T", "T");
    if (!chain->AddFileInfoList(&fFiles)) return false;
    
    fTree   = chain;
    return true;
  }

  /** 
   * Run the job
   * 
   */
  void Run()
  {
    Init(true);
    if (!fTree) {
      Error("Run", "No input tree");
      return;
    }
      
    fFirst = 0xFFFFFFFF;
    fLast  = 0;

    UInt_t nEntries = fTree->GetEntries();
    UInt_t j = 0;
    fRuns.Set(nEntries);
    Info("Run", "Got %d runs", nEntries);
    for (UInt_t i = 0; i < nEntries; i++) {
      fTree->GetEntry(i);

      UInt_t run = fGlobal->runNo;
      UInt_t nev = fGlobal->nAccepted;

      fFirst = TMath::Min(run, fFirst);
      fLast  = TMath::Max(run, fLast);
      fRuns[i] = run;

      Info("Run", "Got run %d with %d accepted events", run, nev);
      fNAccepted->SetPoint(i, run, nev);
      fVz->SetPoint(i, run, fGlobal->meanVz);
      fVz->SetPointError(i, 0, fGlobal->sigmaVz);

      if (nev <= 100) continue;
      static_cast<Ring*>(fFMD1i)->Update(j, run);
      static_cast<Ring*>(fFMD2i)->Update(j, run);
      static_cast<Ring*>(fFMD2o)->Update(j, run);
      static_cast<Ring*>(fFMD3i)->Update(j, run);
      static_cast<Ring*>(fFMD3o)->Update(j, run);
      j++;
    }
    
    Plot();
  }
  /** 
   * Plot results
   * 
   */
  void Plot()
  {
    fTeXName = Form("trend_%09d_%09d", fFirst, fLast);
    TString title;
    if (fYear != 0 && fLetter != '\0') 
      title.Form("QA trends for LHC%d%c runs %d --- %d", 
		 fYear, fLetter, fFirst, fLast);
    else
      title.Form("QA trends for runs %d --- %d", fFirst, fLast);
    MakeCanvas(title);

    CanvasTitle("# of accepted events");
    fNAccepted->Draw("apl");
    PutCanvasTitle("# of accepted events");
    AddRuns(fNAccepted->GetHistogram(), "# of accepted events");

    TLine* l = new TLine(fFirst, 100, fLast, 100);
    l->SetLineColor(kRed+2);
    l->SetLineStyle(2);
    l->Draw();
    PrintCanvas("nAccepted");

    CanvasTitle("#LTv_{z}#GT");
    fVz->Draw("apl");
    PutCanvasTitle("Mean z coordinate of interaction point");
    AddRuns(fVz->GetHistogram(), "#LTv_{z}#GT");
    PrintCanvas("vz");

    TMultiGraph* chi2		= new TMultiGraph;		
    TMultiGraph* c		= new TMultiGraph;		
    TMultiGraph* delta		= new TMultiGraph;		
    TMultiGraph* xi		= new TMultiGraph;		
    TMultiGraph* sigma		= new TMultiGraph;		
    TMultiGraph* low		= new TMultiGraph;		
    TMultiGraph* singles	= new TMultiGraph;	
    TMultiGraph* loss           = new TMultiGraph;
    TMultiGraph* occ	= new TMultiGraph;	
    TMultiGraph* beta		= new TMultiGraph;		
    chi2	->SetName("chi2");	       
    c		->SetName("c");		       
    delta	->SetName("delta");	       
    xi		->SetName("xi");	       
    sigma	->SetName("sigma");	       
    low		->SetName("low");	       
    singles	->SetName("singles");    
    loss        ->SetName("loss");
    beta	->SetName("beta");          
    occ	->SetName("occupancy");


    AddToMulti(fFMD1i,chi2, c, delta, xi, sigma, low, singles, loss, beta, occ);
    AddToMulti(fFMD2i,chi2, c, delta, xi, sigma, low, singles, loss, beta, occ);
    AddToMulti(fFMD2o,chi2, c, delta, xi, sigma, low, singles, loss, beta, occ);
    AddToMulti(fFMD3i,chi2, c, delta, xi, sigma, low, singles, loss, beta, occ);
    AddToMulti(fFMD3o,chi2, c, delta, xi, sigma, low, singles, loss, beta, occ);

    PlotMulti(chi2, 	"#LT#chi^{2}/#nu#GT from #Delta fits", true);
    PlotMulti(c, 	"#LTc#GT from #Delta fits");
    PlotMulti(delta, 	"#LT#Delta_{p}#GT from #Delta fits");
    PlotMulti(xi, 	"#LT#xi#GT from #Delta fits");
    PlotMulti(sigma, 	"#LT#sigma#GT from #Delta fits");
    PlotMulti(low, 	"Bins with too low statistics");
    PlotMulti(singles, 	"Fraction of single hits");
    PlotMulti(loss,     "% of hits 'lost' due to merging+cuts");
    PlotMulti(occ,      "#LTOccupancy#GT [%]", true);
    PlotMulti(beta, 	"Correlation of methods");

    fStore->cd();
    fNAccepted->Write();
    chi2->Write();
    c->Write();
    delta->Write();
    sigma->Write();
    low->Write();
    singles->Write();
    loss->Write();
    occ->Write();
    beta->Write();

    Close(false); // Do not delete PNGs
  }
  /** 
   * Add graphs from a ring to the multi graphs
   * 
   * @param qr        Ring
   * @param chi2      ELoss reduced chi-square	   
   * @param c 	      ELoss constant 		   
   * @param delta     ELoss MPV		   
   * @param xi 	      ELoss Landau width	   
   * @param sigma     ELoss Gaus width		   
   * @param low       bins with low statistics	   
   * @param singles   fraction of singles	   
   * @param loss      'lost' data 		   
   * @param beta      Poisson vs ELoss correlation
   * @param occupancy mean occupancy              
   */
  void AddToMulti(QARing* qr, 
		  TMultiGraph* chi2,
		  TMultiGraph* c,
		  TMultiGraph* delta,
		  TMultiGraph* xi,
		  TMultiGraph* sigma,
		  TMultiGraph* low,
		  TMultiGraph* singles,
		  TMultiGraph* loss,
		  TMultiGraph* beta,
		  TMultiGraph* occupancy)
  {
    Ring* r = static_cast<Ring*>(qr);
    chi2	->Add(r->fGChi2);
    c		->Add(r->fGC);
    delta	->Add(r->fGDelta);
    xi		->Add(r->fGXi);
    sigma	->Add(r->fGSigma);
    low		->Add(r->fGLow);
    singles	->Add(r->fGSingles);
    loss        ->Add(r->fGLoss);
    occupancy	->Add(r->fGOccupancy);
    beta	->Add(r->fGBeta);
  }
  /** 
   * Plot a multi-graph
   * 
   * @param mg     Multi graph
   * @param title  Title
   * @param logy   If true, make @f$\log@f$ scale 
   */
  void PlotMulti(TMultiGraph* mg, const char* title, Bool_t logy=false)
  {
    CanvasTitle(title);
    // fCanvas->SetBottomMargin(.15);
    fCanvas->SetLeftMargin(.08);
    fCanvas->SetTopMargin(.1);
    fCanvas->SetLogy(logy);
    // fCanvas->SetRightMargin(.2);
    mg->Draw("apl");
    
    TH1*     h   = mg->GetHistogram();
    Double_t max = h->GetMaximum();
    Double_t min = h->GetMinimum();
    if (h->GetMinimum() == 0) {
      min = min - .1*(max-min);
      h->SetMinimum(min);
    }
    Int_t    x1  = h->GetXaxis()->GetXmin();
    Int_t    x2  = h->GetXaxis()->GetXmax();
    
    TLegend* l = new TLegend(.1, .91, .97, .95);
    l->SetNColumns(5);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(42);
    TIter next(mg->GetListOfGraphs());
    mg->GetListOfGraphs();
    TGraph* g = 0;

    // Get the runs we have here 
    TArrayI runs(fRuns.GetSize());
    runs.Reset(INT_MAX);
    Int_t   j = 0;
    while ((g = static_cast<TGraph*>(next()))) { 
      l->AddEntry(g, g->GetTitle(), "lp");
      Double_t* xs = g->GetX();
      Int_t     n  = g->GetN();
      for (Int_t i = 0; i < n; i++) {
	if (FindRun(runs, Int_t(xs[i])) >= 0) continue;
	runs.SetAt(xs[i], j++);
      }
      Double_t ymean = g->GetMean(2);
      Double_t xh    = x2 - .03 * Double_t(x2-x1); 
      TLine* lm = new TLine(x1, ymean, xh, ymean);
      lm->SetLineColor(g->GetLineColor());
      lm->SetLineStyle(2);
      lm->SetLineWidth(1);
      lm->Draw();
      TLatex* al = new TLatex(xh, ymean, g->GetTitle());
      al->SetTextColor(g->GetLineColor());
      al->SetTextFont(42);
      al->SetTextSize(.02);
      al->SetTextAlign(12);
      al->Draw();
    }
    l->Draw();

    TList areas;
    areas.SetOwner();
    AddRuns(h, title, &runs, &areas);

    PrintCanvas(mg->GetName(), &areas);
  }
  /** 
   * Find a run 
   * 
   * @param runs List of runs 
   * @param run  Run to find
   * 
   * @return Index of run in run list
   */
  Int_t FindRun(const TArrayI& runs, Int_t run) 
  {
    std::sort(&(runs.fArray[0]), &(runs.fArray[runs.GetSize()]));
    Int_t idx = TMath::BinarySearch(runs.GetSize(), runs.fArray, run);
    if (idx >= runs.GetSize() || idx < 0 || runs[idx] != run) return -1;
    return idx;
  }
  /** 
   * Add run labels at appropriate places on the plot
   * 
   * @param h      Frame histogram
   * @param title  Title 
   * @param runs   List of runs, if any
   * @param areas  Other areas 
   */
  void AddRuns(TH1* h, const char* title, TArrayI* runs=0, 
	       TList* areas=0)
  {
    h->GetXaxis()->SetNoExponent();
    // h->GetXaxis()->SetTitleOffset(1);
    TString ytitle(title);
    if (fUseVar) ytitle.Append(" (errors: variance)");
    else         ytitle.Append(" (errors: min/max)");
    h->SetYTitle(ytitle.Data());
    h->SetXTitle("Run #");
    Info("AddRuns", "%s: %s vs %s", h->GetName(), 
	 h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());

    Int_t    r1  = h->GetXaxis()->GetXmin();
    Int_t    r2  = h->GetXaxis()->GetXmax();
    Double_t lx  = 0;
    Double_t tx  = .025; // (r2 - r1) / 18;
    Double_t wx  = 1 - fCanvas->GetLeftMargin() - fCanvas->GetRightMargin();
    Double_t dy  = .025;
    Double_t y   = fCanvas->GetBottomMargin()+dy;
    for (Int_t i = 0; i < fRuns.GetSize(); i++) {
      Int_t    r = fRuns[i];
      Double_t x = fCanvas->GetLeftMargin() + wx*Double_t(r-r1)/(r2-r1);

      // Skip runs out of range 
      if (r < r1 || r > r2) continue;

      // Skip runs not in the graphs 
      if (runs) {
	if (FindRun(*runs, r) < 0) { 
	  continue;
	}
      }
	  

      if (TMath::Abs(x - lx) < tx) y += dy;
      else                         y =  fCanvas->GetBottomMargin() + dy;


      // Info(h->GetName(), "%6d (x,y)=(%12f,%12f) |lx-x|=|%f-x|=%f", 
      //      r, x, y, lx, TMath::Abs(lx-x));
      lx = x;
      
      Double_t* xa = fNAccepted->GetX();
      Int_t     na = fNAccepted->GetN();
      Int_t idx = TMath::BinarySearch(na, xa, Double_t(r));
      Color_t color = kBlue+3;
      if (idx >= 0  && idx < na  && r == xa[idx] && 
	  fNAccepted->GetY()[idx] < 10000) 
	color = kRed+3;

      TLatex* ll = new TLatex(x, y, Form("%d", r));
      ll->SetNDC();
      ll->SetTextAlign(21);
      ll->SetTextSize(0.02);
      ll->SetTextColor(color);
      ll->SetTextFont(42);
      // ll->SetTextAngle(90);
      ll->Draw();
      TLine* tl = new TLine(x, y, x, 1-fCanvas->GetTopMargin());
      tl->SetBit(TLine::kLineNDC);
      tl->SetLineStyle(3);
      tl->SetLineColor(color);
      tl->Draw();

      if (!areas) continue;
      
      TObjString* area = new TObjString;
      TString&    spec = area->String();
      spec.Form("<span style=\"left: %d%%; bottom: %d%%;\" "
		"onClick='window.location=\"qa_%09d.html\"' "
		"onMouseOver='this.style.cursor=\"pointer\"' "
		"alt=\"%d\""
		">%d</span>",
		UInt_t(100*x)-2, UInt_t(100*y), r, r, r);
      
#if 0
      spec.Form("<area shape='rect' alt='%d' title='%d' href='qa_%09d.html' "
		"coords='%d,%d,%d,%d'>", r, r, r, 
		UInt_t(cw*(x-tx)), 0, UInt_t(cw*(x+tx)), ch);
#endif
      areas->Add(area);
    }
  }
  /** 
   * Output list of runs 
   * 
   * @param o Output stream
   */
  void WriteRuns(std::ostream& o)
  {
    o << "<div class='jobid'><!--JOBID--></div>\n"
      << "<div class='runs'>\n"
      << "  Runs:<br> ";
    for (Int_t i = 0; i < fRuns.GetSize(); i++) {
      o << "<a href='qa_" << Form("%09d", fRuns[i]) << ".html'>"
	<< fRuns[i] << "</a> " << std::flush;
    }
    o << "\n"
      << "</div>" << std::endl;
  }
  /** 
   * Write page footer 
   * 
   */
  void WriteFooter()
  {
    WriteRuns(*fHtml);
    QABase::WriteFooter();
  }
  /** 
   * Write out image footer 
   * 
   * @param o Output stream
   * @param pngName Name of the PNG file 
   */
  void WriteImageFooter(std::ostream& o, const char* pngName)
  {
    WriteRuns(o);
    QABase::WriteImageFooter(o, pngName);
  }
  TGraph*       fNAccepted; // Graph of number of accepted events
  TGraphErrors* fVz;        // Graph of mean vertex
  UInt_t        fFirst;     // First run
  UInt_t        fLast;      // Last run
  TArrayI       fRuns;      // Seen runs 
  TObjArray     fFiles;
  Bool_t        fUseVar;    // Use variance rather than min/max 
};
 
//
// EOF
//

  
  
