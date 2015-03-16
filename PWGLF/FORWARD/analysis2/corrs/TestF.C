#include "AliLandauGaus.h"
#include <TSystem.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TFitResult.h>
#include <TNtuple.h>
#include <TAxis.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPaveStats.h>

//====================================================================
/** 
 * Base class for tests
 * 
 */
struct TestF
{
  //__________________________________________________________________
  /** 
   * Base class for individual tests
   */
  struct Test { 
    /** Source distribution */
    TF1* src;
    /** 
     * Constructor 
     */
    Test() : src(0) {}
    /** 
     * Copy constructor 
     */
    Test(const Test& o) 
      : src(static_cast<TF1*>(o.src ? o.src->Clone() : 0)) {}
    /** 
     * Destructor 
     */
    virtual ~Test() { delete src; }
    /** 
     * Assignment operator  
     */
    Test& operator=(const Test& o) {
      if (&o == this) return *this;
      if (o.src) src = static_cast<TF1*>(o.src->Clone());
      return *this;
    }
    /** 
     * Run the test
     * 
     * @param c     Constant
     * @param delta Most probable value 
     * @param xi    Landau width
     * @param sigma Gaussian smear
     * @param n     Maximum number of particles
     * @param a     Particle weights
     * @param xMin  Least energy loss 
     * @param xMax  Largest energy loss
     */
    void Run(const Double_t  c, 
	     const Double_t  delta, 
	     const Double_t  xi, 
	     const Double_t  sigma, 
	     const Int_t     n, 
	     const Double_t* a,
	     const Double_t  xMin, 
	     const Double_t  xMax) {
      printf("%8.4f | %8.4f | %8.4f | ", delta, xi, sigma);
      fflush(stdout);
      src = AliLandauGaus::MakeFn(c,delta,xi,sigma,0,n,&(a[1]),xMin,xMax);
      src->SetName("src");
      src->SetLineColor(kBlack);
      src->SetNpx(10000);
      // If the constant is 1, then this is a fake distribtion and we
      // need to scale it. 
      if (c >= 1) { 
	Double_t max = src->GetMaximum();
	src->SetParameter(0, c/max);
      }

      // src->Print();
      DoRun();
    }
    /** 
     * @return Source parameter @f$ C@f$
     */
    Double_t C()      const { return src->GetParameter(AliLandauGaus::kC); }
    /** 
     * @return Source parameter @f$\Delta_p@f$ 
     */
    Double_t Delta()  const { return src->GetParameter(AliLandauGaus::kDelta); }
    /** 
     * @return Source parameter @f$\xi@f$ 
     */
    Double_t Xi()     const { return src->GetParameter(AliLandauGaus::kXi); }
    /** 
     * @return Source parameter @f$ \sigma@f$ 
     */
    Double_t Sigma()  const { return src->GetParameter(AliLandauGaus::kSigma); }
    /** 
     * @return Source parameter @f$ \sigma_n@f$
     */
    Double_t SigmaN() const { return src->GetParameter(AliLandauGaus::kSigmaN);}
    /** 
     * @return Source parameter @f$ N@f$
     */
    Double_t N()      const { return src->GetParameter(AliLandauGaus::kN); }
    /** 
     * @return Source parameter @f$ a_i@f$ 
     */
    Double_t A(Int_t i) const { 
      return (i <= 1 ? 1 : src->GetParameter(AliLandauGaus::kN+(i-1))); }
    /** 
     * @return Least energy loss
     */
    Double_t XMin() const { return src->GetXmin(); }
    /** 
     * @return Largest energy loss
     */
    Double_t XMax() const { return src->GetXmax(); }
    /** 
     * Generate single component function for @a i particles
     * 
     * @param c      Constant 
     * @param delta  @f$\Delta_p@f$ 
     * @param xi     @f$\xi@f$ 
     * @param sigma  @f$\sigma@f$ 
     * @param sigmaN @f$\sigma_n@f$
     * @param i      Number of particles 
     * 
     * @return Pointer to function object 
     */
    TF1* Fi(Double_t c, Double_t delta, Double_t xi, Double_t sigma, 
	    Double_t sigmaN, Int_t i) const { 
      return AliLandauGaus::MakeFi(c,delta,xi,sigma,sigmaN,i,XMin(),XMax());
    }
    /** 
     * Generate single component function for @a i particles using
     * source parameters
     * 
     * @param i      Number of particles 
     * 
     * @return Pointer to function object 
     */
    TF1* Fi(Int_t i) { 
      return Fi(C()*A(i),Delta(),Xi(),Sigma(),SigmaN(),i);
    }
    
    /** 
     * @{ 
     * @name Overload-able interface 
     */
    /** 
     * Run the test
     */
    virtual void DoRun() = 0;
    /** 
     * Write results to disk 
     * 
     * @param d Directory to write to 
     */
    virtual void WriteOut(TDirectory* d) = 0;
    /** 
     * Draw results in a pad 
     * 
     * @param p Pad to draw in 
     */
    virtual void DrawInPad(TVirtualPad* p) = 0;
    /* @} */
  };
  //__________________________________________________________________
  /** 
   * @{
   * @name Utilities 
   */
  /** 
   * Draw a vertical line 
   * 
   * @param i    Line style
   * @param col  Line color
   * @param x    X value 
   * @param low  Least Y value 
   * @param high Largest Y value 
   * 
   * @return The line 
   */
  static TLine*
  VerticalLine(Int_t i, Color_t col, Double_t x, Double_t low, Double_t high)
  {
    TLine* l = new TLine(x,low,x,high);
    l->SetLineStyle(i);
    l->SetLineColor(col);
    return l;
  }
  //__________________________________________________________________
  virtual ~TestF() {}
  /** 
   * Make a test object - overload
   * 
   * @return The test object
   */
  virtual Test* MakeTest() = 0;
  /** 
   * Do one scan 
   * 
   * @param c        Constant
   * @param delta    Most probable value @f$\Delta_p@f$
   * @param xi       Landau width @f$\xi@f$
   * @param sigma    Gaussian spread @f$\sigma@f$
   * @param n        Maximum number of particles
   * @param a        Particle weights 
   * @param xMin     Least @f$\Delta @f$ 
   * @param xMax     Largest @f$\Delta@f$
   * @param pad      Pad to draw in 
   * 
   * @return The test object
   */
  Test* DoOne(const Double_t  c, 
	      const Double_t  delta, 
	      const Double_t  xi, 
	      const Double_t  sigma, 
	      const Int_t     n, 
	      const Double_t* a, 
	      const Double_t  xMin, 
	      const Double_t  xMax,
	      TVirtualPad*    pad=0)
  {
    Test* t = MakeTest();
    t->Run(c, delta, xi, sigma, n, a, xMin, xMax);

    TString name(Form("sigma%04d_xi%04d", Int_t(100*sigma), Int_t(100*xi)));
    TString title(Form("Sigma=%f xi=%f", sigma, xi));
    if (!pad) { 
      TCanvas* can = new TCanvas(name, title, 1200, 900);
      pad = can;
    }
    else 
      pad->Clear();
    t->DrawInPad(pad);

    return t;
  }
  /** 
   * Do a unit scan 
   * 
   * @param sigma @f$\sigma@f$ of the Gaussian
   * @param xi    @f$\xi@f$ width of Landau
   * @param n     Number of particles 
   * @param pad   Possibly existing canvas to plot in 
   *
   * @return Scan object
   */
  Test* DoUnit(const Double_t sigma, 
	       const Double_t xi=1, 
	       const Int_t    n=10,
	       TVirtualPad*   pad=0) 
  { 
    const Double_t c      = 1;
    const Double_t delta  = 10;
    Double_t       idelta = delta;
    Double_t       ixi    = xi;
    Double_t       isigma = sigma;
    AliLandauGaus::IPars(n+1,idelta,ixi, isigma);
    Double_t       up     = idelta + ixi + isigma; // 8*n
    Double_t       a[n];
    a[0] = 1;
    for (Int_t i = 1; i < n; i++) 
      a[i] = TMath::Exp(0.5 - TMath::Sqrt(2)*(i+1));

    return DoOne(c,delta,xi,sigma,n,a,-3,up, pad);
  }
  /** 
   * Do a realistic scan 
   * 
   * @param pad   Possibly existing canvas to plot in 
   * 
   * @return Scan object
   */
  Test* DoRealistic(TVirtualPad* pad=0) 
  {
    Double_t c     = 0.4070;
    Double_t delta = 0.5651;
    Double_t xi    = 0.0643;
    Double_t sigma = 0.0611;
    Int_t    n     = 5;
    Double_t a[]   = {1, 0.1179, 0.0333, 0.0068, 0.0012 };

    return DoOne(c, delta, xi, sigma, n, a, 0.02, n, pad);
  }
  /** 
   * Prefix on output 
   * 
   * 
   * @return 
   */
  virtual const char* Prefix() const = 0;
  /** 
   * Make an output ntuple 
   * 
   * @return Possibly newly allocated NTuple 
   */
  virtual TNtuple* MakeNTuple() { return 0; }
  /** 
   * Setup for scanning 
   * 
   * @param mode   Scanning mode 
   * @param file   Output file 
   * @param canvas Canvas 
   * @param nt     Possible Ntuple 
   */
  virtual void SetupScan(UShort_t  mode, 
			 TFile*&   file, 
			 TCanvas*& canvas, 
			 TNtuple*& nt) 
  {
    TString  base = Form("%s%s", Prefix(), 
			 (mode == 0 ? "SigmaXi" : 
			  mode == 1 ? "Sigma" : "Xi"));
    file   = TFile::Open(Form("%s.root", base.Data()), "RECREATE");
    nt     = MakeNTuple();
    canvas = new TCanvas(base, Form("%s.pdf", base.Data()),
			       TMath::Sqrt(2)*900, 900);
    TString title("pdf Landscape Title:Scan");
    gSystem->RedirectOutput("/dev/null");
    canvas->Print(Form("%s[", canvas->GetTitle()), title);
    gSystem->RedirectOutput(0);
  }
  /** 
   * Print canvas to output PDF 
   * 
   * @param c      Canvas
   * @param xi     @f$ \xi@f$ 
   * @param sigma  @f$ \sigma@f$ 
   * 
   * @return The full title of the page 
   */
  void PrintCanvas(TCanvas* c, Double_t xi, Double_t sigma) {
    TString tit = "pdf Landscape Title:";
    Bool_t  cls = false;
    if (xi > 0 && sigma > 0) {
      tit.Append(Form("xi=%8.4f sigma=%8.4f", xi, sigma));
    }
    else { 
      tit.Append("Results");
      cls = true;
    }
    gSystem->RedirectOutput("/dev/null");
    c->Print(c->GetTitle(), tit);
    if (cls) {
      c->Clear();
      c->Print(Form("%s]", c->GetTitle()), tit);
    }
    gSystem->RedirectOutput(0);
    if (cls)
      ::Info("", "Plots saved in %s", c->GetTitle());
  }
  /** 
   * Executed before loop 
   * 
   * @param mode  Mode
   * @param nVal  Number of values  
   */
  virtual void PreLoop(UShort_t mode, Int_t nVal) = 0;
  /** 
   * Process a single step 
   * 
   * @param mode   Mode 
   * @param i      First index 
   * @param j      Second index (if any)
   * @param n      Number of values
   * @param test   Test result
   * @param nt     NTutple 
   */
  virtual void Step(UShort_t mode, 
		    Int_t    i, 
		    Int_t    j, 
		    Int_t    n, 
		    Test*    test,
		    TNtuple* nt) = 0;
  /** 
   * Called after looping 
   * 
   * @param mode Execution mode 
   * @param can  Canvas
   * @param nt   NTuple
   * @param out  Output file 
   */
  virtual void PostLoop(UShort_t mode, 
			TCanvas* can, 
			TNtuple* nt, 
			TFile*   out) = 0;
  /** 
   * Scan over a parameter (@f$\xi@f$ or @f$\sigma@f$) 
   * 
   * @param scanSigma   If true, scan @f$\sigma@f$, otherwise @f$\xi@f$ 
   * @param n           Number of values
   * @param values      Values of the parameters
   * @param maxN        Number of particles 
   */
  void ScanOne(Bool_t          scanSigma, 
	       Int_t           n, 
	       const Double_t* values, 
	       Int_t           maxN=10) 
  {
    UShort_t mod = scanSigma ? 1 : 2;
    Scan(mod, n, values, maxN);
  }
  /** 
   * Scan over a parameter (@f$\xi@f$ or @f$\sigma@f$) 
   * 
   * @param scanSigma   If true, scan @f$\sigma@f$, otherwise @f$\xi@f$ 
   * @param values      Values of the parameters
   * @param maxN        Number of particles 
   */
  void ScanOne(Bool_t         scanSigma, 
	       const TArrayD& values, 
	       Int_t          maxN)
  {
    ScanOne(scanSigma, values.GetSize(), values.GetArray(), maxN);
  }
  /** 
   * Scan over both parameter (@f$\xi@f$ and @f$\sigma@f$) 
   * 
   * @param n       Number of parameters 
   * @param values  Values of the parameters
   * @param maxN    Number of particles 
   */
  void ScanTwo(Int_t           n, 
	       const Double_t* values, 
	       Int_t           maxN=10)
  {
    Scan(0, n, values, maxN);
  }
  /** 
   * Scan over both parameter (@f$\xi@f$ and @f$\sigma@f$) 
   * 
   * @param values  Values of the parameters
   * @param maxN    Number of particles 
   */
  void ScanTwo(const TArrayD& values, Int_t maxN)
  {
    ScanTwo(values.GetSize(), values.GetArray(), maxN);
  }
  /** 
   * Do the scan 
   * 
   * @param mode    Mode of operation (0: both, 1: @f$\sigma@f$, or @f$\xi@f$)
   * @param n       Number of values 
   * @param values  Values 
   * @param maxN    Maximum number of particles 
   */
  void Scan(UShort_t        mode, 
	    Int_t           n, 
	    const Double_t* values, 
	    Int_t           maxN)
  {
    TFile*   out = 0;
    TNtuple* nt  = 0;
    TCanvas* can = 0;
    SetupScan(mode, out, can, nt);
    PreLoop(mode, n);

    Test* t  = DoRealistic(can);
    t->WriteOut(out->mkdir("realistic"));
    can->Write();
    PrintCanvas(can, t->Xi(), t->Sigma());
    
    if (mode == 0) LoopTwo(n, values, maxN, can, nt, out);
    else           LoopOne(mode, n, values, maxN, can, nt, out);
    
    PostLoop(mode, can, nt, out);
    out->Write();
    out->Close();

  }
  /** 
   * Loop over one of @f$\xi@f$ or @f$\sigma@f$ 
   * 
   * @param mode    Mode (1: @f$\sigma@f$, otherwise @f$\xi@f$)
   * @param n       Number of values
   * @param values  Values
   * @param maxN    Maximum number of particles 
   * @param can     Canvas to draw in 
   * @param nt      NTuple to fill 
   * @param out     Output directory 
   */
  void LoopOne(UShort_t        mode, 
	       Int_t           n, 
	       const Double_t* values, 
	       Int_t           maxN,
	       TCanvas*        can, 
	       TNtuple*        nt, 
	       TFile*          out) {
    Bool_t      scanSigma = mode == 1;
    const char* var       = (scanSigma ? "sigma" : "xi");
	
    for (Int_t i = 0; i < n; i++) { 
      Double_t    v     = values[i];
      Double_t    sigma = (scanSigma  ? v : 1);
      Double_t    xi    = (!scanSigma ? v : 1);

      StepIt(mode, xi, sigma, maxN, 0, i, n, can, nt, out, var, v);
    }
  }
  /** 
   * Loop over both @f$\xi@f$ and @f$\sigma@f$ 
   * 
   * @param n       Number of values
   * @param values  Values
   * @param maxN    Maximum number of particles 
   * @param can     Canvas to draw in 
   * @param nt      NTuple to fill 
   * @param out     Output directory 
   */
  void LoopTwo(Int_t           n, 
	       const Double_t* values, 
	       Int_t           maxN,
	       TCanvas*        can, 
	       TNtuple*        nt, 
	       TFile*          out) 
  { 
    UShort_t mode = 0;
    for (Int_t i = 0; i < n; i++) { 
      Double_t      xi  = values[i];
      TDirectory*   d   = out->mkdir(Form("xi%04d", Int_t(100*xi)));
      for (Int_t j = 0; j < n; j++) { 
	Double_t    sigma = values[j];
	StepIt(mode, xi, sigma, maxN, i, j, n, can, nt, d, "sigma", sigma);
      }
    }
  }
  /** 
   * 
   * 
   * @param mode   Mode
   * @param xi     @f$\xi@f$ 
   * @param sigma  @f$\sigma@f$ 
   * @param maxN   At most this many particles
   * @param i      First index
   * @param j      Last index 
   * @param n      Number of values
   * @param can    Canvas 
   * @param nt     NTuple 
   * @param p      Parent directory
   * @param pre    Prefix on directory
   * @param v      Current value 
   */
  void StepIt(UShort_t    mode, 
	      Double_t    xi, 
	      Double_t    sigma, 
	      Int_t       maxN, 
	      Int_t       i, 
	      Int_t       j, 
	      Int_t       n, 
	      TCanvas*    can, 
	      TNtuple*    nt, 
	      TDirectory* p, 
	      const char* pre, 
	      Double_t    v)
  {
    TDirectory* d     = p->mkdir(Form("%s%04d", pre, Int_t(100*v)));
    Test*       t     = DoUnit(sigma, xi, maxN, can);
    t->WriteOut(d);
    can->Write();
    PrintCanvas(can, xi, sigma);
      
    Step(mode, i, j, n, t, nt);
  }
};

//====================================================================
#ifdef TEST_SHIFT
#include "WrappedGraph.C"
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TMultiGraph.h>
#include <Math/RootFinder.h>

struct TestShift : public TestF
{
  /** 
   * Our fit function
   *
   * @f[
   *   \Delta_{i,p} - \Delta_{i,r} = 
   *     f(i;c,p,\xi,\sigma) = 
   *     \frac{c u \sigma}{(1+1/i)^{p u^{3/2}}
   * @f]
   *
   * where @f$ u=\sigma/\xi@f$ 
   * 
   * @param xp Independent variables 
   * @param pp Parameters
   * 
   * @return @f$ f(i;c,p,\xi,\sigma)@f$ 
   */
  static Double_t fitFunc(Double_t* xp, Double_t* pp) {
    Double_t x     = *xp;
    Double_t c     = pp[0];
    Double_t p     = pp[1];
    Double_t xi    = pp[2];
    Double_t sigma = pp[3];
    Double_t u     = sigma / xi;
    Double_t a     = c * u * sigma;
    Double_t b     = p * u * TMath::Sqrt(u);
      
    return a / TMath::Power(1+1./x, b);
  }
  /** 
   * Scale a graph 
   * 
   * @param g     Graph to scale
   * @param scale Scale 
   */
  static void ScaleGraph(TGraph* g, Double_t scale)
  {
    // Info("", "Scaling graph %p with %f", g, scale);
    if (scale == 1) return;
    for (Int_t i = 0; i < g->GetN(); i++) 
      g->SetPoint(i,g->GetX()[i], g->GetY()[i] * scale);
  }
  //__________________________________________________________________
  struct Test : public TestF::Test 
  {
    TGraph*  zeros;       // Zeros of the derivatives
    TGraph*  diff;        // Difference from calculated Delta to roots
    TList    comps;       // Components 
    TList    derivs;      // Derivatives 
    TF1*     fit;         // Fit to difference 
    Test() : TestF::Test(), zeros(0), diff(0), fit(0) 
    {
      comps.SetName("components");
      derivs.SetName("derivatives");

    }
    /** 
     * Run the actual test
     */
    void DoRun() 
    { 
      AliLandauGaus::EnableSigmaShift(false); // Disable sigma-shift
      const Int_t n = N();
      zeros = new TGraph(n);
      zeros->SetName("zeros");
      zeros->SetTitle("Zeros of derivatives");
      zeros->SetMarkerStyle(25);
      zeros->SetMarkerSize(1.2);
      zeros->SetMarkerColor(kRed+1);

      diff = new TGraph(n);
      diff->SetName("diff");
      diff->SetTitle("#delta#Delta_{p}");
      diff->SetMarkerStyle(20);
      diff->SetMarkerColor(kBlue+1);

      // Make the components 
      Int_t j = 0;
      for (Int_t i = 1; i <= n; i++, j++) 
	if (!MakeComp(i)) break;
      zeros->Set(j);
      diff->Set(j);

      // Now fit 
      fit = new TF1("fit", &fitFunc, 1, j, 4);
      fit->SetParNames("c", "p", "#xi", "#sigma");
      fit->SetParameter(0, .5);
      fit->SetParameter(1, .5);
      fit->FixParameter(2, Xi());
      fit->FixParameter(3, Sigma());
      
      // fit->SetParNames("c", "#xi", "#sigma");
      // fit->FixParameter(1, xi);


      Int_t ifit = 4;
      do { 
	gSystem->RedirectOutput("/dev/null");
	diff->Fit(fit, Form("MR%s", ifit == 0 ? "Q" : "Q"));
	gSystem->RedirectOutput(0);
	ifit--;
      } while (ifit >= 0);
      Printf("%8.4f | %8.4f %9.5f | %8.4f %9.5f",
	     fit->GetChisquare() / fit->GetNDF(),
	     fit->GetParameter(0), fit->GetParError(0),
	     fit->GetParameter(1), fit->GetParError(1));

    }
    /** 
     * Make a single component 
     * 
     * @param i        Number of particles
     */
    Bool_t MakeComp(const Int_t    i) {
      TF1* fi = Fi(i);
      fi->SetNpx(src->GetNpx());
      // fi->Print();

      TGraph* dF = new TGraph(fi, "d");
      dF->SetName(Form("d%s", fi->GetName()));
      dF->SetTitle(Form("d%s/d#Delta", fi->GetTitle()));
      dF->SetMarkerStyle(1); // A dot 
      // dF->SetMarkerStyle(20+i-1);
      dF->SetMarkerColor(fi->GetLineColor());
      dF->SetLineColor(fi->GetLineColor());

      Double_t max = TMath::MaxElement(dF->GetN(), dF->GetY());
      if (max < 1e-6) return false;

      ScaleGraph(dF, 1./max);
      Double_t delta  = Delta();
      Double_t xi     = Xi();
      Double_t deltaI = i * (delta + xi * TMath::Log(i));
      Double_t xiI    = i * xi;
      Double_t xR     = FindRoot(dF, deltaI-2*xiI, deltaI+2*xiI);
      
      // Printf("Component %2d: Exp: %8.4f  Got: %8.4f  Diff: %9.5f", 
      //        i, deltaI, xR, (xR-deltaI));
      // printf(".");
      comps.Add(fi);
      derivs.Add(dF);
      zeros->SetPoint(i-1, xR, dF->Eval(xR));
      diff->SetPoint(i-1, i, xR - deltaI);

      return true;
    }
    /** 
     * Find the root of a graph
     * 
     * @param g     Graph
     * @param xMin  Least value of scan range 
     * @param xMax  Largest value of scan range 
     * 
     * @return Zero of the graph
     */
    Double_t FindRoot(TGraph* g, Double_t xMin, Double_t xMax) { 
      WrappedGraph* wG = new WrappedGraph(g);
      ROOT::Math::RootFinder rfb(ROOT::Math::RootFinder::kBRENT);
      rfb.SetFunction(*wG, xMin, xMax); 
      rfb.Solve();
      Double_t xR = rfb.Root();
      return xR;
    }
    /** 
     * Draw the results in a pad 
     * 
     * @param body Pad to draw in 
     */
    void DrawInPad(TVirtualPad* body) { 
      body->SetTopMargin(0.01);
      body->SetRightMargin(0.01);
      body->Divide(2,1);
      
      TVirtualPad* p = body->cd(1);
      p->SetTopMargin(0.01);
      p->SetRightMargin(0.01);
      p->Divide(1,2,0,0);
      
      TVirtualPad* q = p->cd(1);
      q->SetRightMargin(0.01);
      q->SetGridx();
      q->SetLogy();
      src->Draw();
      Double_t fmax = src->GetHistogram()->GetMaximum();
      Double_t fmin = src->GetHistogram()->GetMinimum();
      fmin          = TMath::Max(fmin, fmax/1e8);
      src->GetHistogram()->SetMinimum(fmin);
      TLatex*  ltx  = new TLatex(q->GetLeftMargin()+.02, 
				 1-q->GetTopMargin()-.02,
				 Form("#Delta_{p}=%f, #xi=%f, #sigma=%f", 
				      Delta(), Xi(), Sigma()));
      // Printf("%s: Max=%f scale=%f", src->GetName(), fmax, 1/fmax);
      ltx->SetTextFont(42);
      ltx->SetNDC();
      ltx->SetTextAlign(13);
      ltx->Draw();

      q = p->cd(2);
      q->SetRightMargin(0.01);
      q->SetGridx();
      
      Double_t x = 1-q->GetRightMargin()-.02;
      Double_t y = 1-q->GetTopMargin()-.1;
      for (Int_t i = 1; i <= comps.GetEntries(); i++) { 
	p->cd(1);
	TF1* fi = static_cast<TF1*>(comps.At(i-1));
	fi->Draw("same");

	Double_t deltaR = zeros->GetX()[i-1];
	Double_t deltaI = i * (Delta() + Xi() * TMath::Log(i));

	Color_t col = fi->GetLineColor();
	VerticalLine(1,col,deltaI,fmin,fmax)->Draw();
	VerticalLine(3,col,deltaR,fmin,fmax)->Draw();
	
	ltx = new TLatex(x,y,Form("#Delta_{p,%2d}=%8.4f %8.4f",
				  i,deltaI,deltaR));
	ltx->SetTextAlign(33);
	ltx->SetTextFont(42);
	ltx->SetNDC();
	ltx->Draw();
	y -= ltx->GetTextSize() + .01;

	p->cd(2);

	TGraph* df = static_cast<TGraph*>(derivs.At(i-1));
	df->Draw(i == 1 ? "alp" : "lp");
	df->GetHistogram()->GetXaxis()->SetRangeUser(XMin(), XMax());	
	VerticalLine(1,col,deltaI,-.6,1.1)->Draw();
	VerticalLine(3,col,deltaR,-.6,1.1)->Draw();
      }
      p->cd(2);
      zeros->Draw("p");
      
      gStyle->SetOptFit(111111);
      gStyle->SetStatY(.6);
      p = body->cd(2);
      p->SetRightMargin(0.02);
      p->SetTopMargin(0.02);
      diff->Draw("alp");
      diff->GetHistogram()->GetListOfFunctions()->Add(fit);
      diff->GetHistogram()->SetXTitle("N_{particle}");
      diff->GetHistogram()->SetYTitle("#Delta_{r}-#Delta_{p}");
      p->Clear();
      diff->Draw("alp");

      body->Modified();
      body->Update();
      body->cd();
    }
    /** 
     * Write results to disk 
     * 
     * @param dir Directory
     */
    void WriteOut(TDirectory* dir) { 
      dir->cd();
      src->Write();
      zeros->Write();
      diff->Write();
      comps.Write(comps.GetName(),  TObject::kSingleKey);
      derivs.Write(derivs.GetName(),  TObject::kSingleKey);
      fit->Write();
    }
  };
  //__________________________________________________________________
  TGraphErrors*   cs;
  TGraphErrors*   ps;
  TGraph2DErrors* c2s;
  TGraph2DErrors* p2s;
  TMultiGraph*    mcs;
  TMultiGraph*    mps;
  
  /** 
   * Make a test object - overload
   * 
   * @return The test object
   */
  virtual TestF::Test* MakeTest() { return new Test(); }
  /** 
   * Prefix on output 
   * 
   * @return The string "shift"
   */
  virtual const char* Prefix() const { return "shift"; }
  /** 
   * Make an output ntuple 
   * 
   * @return Possibly newly allocated NTuple 
   */
  virtual TNtuple* MakeNTuple () 
  { 
    return new TNtuple("shift","Shift scan",
		       "Delta:xi:sigma:chi2nu:p:ep:c:ec");
  }
  /** 
   * Set-up for looping 
   * 
   * @param mode 
   * @param nVal 
   */
  virtual void PreLoop(UShort_t mode, Int_t n) 
  { 
    if (mode == 0) {
      c2s = new TGraph2DErrors(n*n);c2s->SetName("cs");
      p2s = new TGraph2DErrors(n*n);p2s->SetName("ps");
      mcs = new TMultiGraph("mcs", "C as function of #sigma/#xi");
      mps = new TMultiGraph("mps", "P as function of #sigma/#xi");
      c2s->SetDirectory(0);
      p2s->SetDirectory(0);
      c2s->SetMarkerStyle(20);
      p2s->SetMarkerStyle(20);
      c2s->SetTitle("c(#xi,#sigma)");
      p2s->SetTitle("p(#xi,#sigma)");
    }
    else {
      const char* var = (mode == 1 ? "#sigma" : "#xi");
      cs = new TGraphErrors(n); cs->SetName("cs");
      ps = new TGraphErrors(n); ps->SetName("ps");
      cs->SetMarkerStyle(20);
      ps->SetMarkerStyle(20);
      cs->SetTitle(Form("c(%s)", var));
      ps->SetTitle(Form("p(%s)", var));
    }
    Printf("%-8s | %-8s | %-8s | %-8s | %-8s %-9s | %-8s %-9s",
	   "Delta_p", "xi", "sigma", "chi^2/nu", "c", "error", "p", "error");
  }
  /** 
   * Fill multi graph
   * 
   * @param mg - multi-graph
   * @param i    First index
   * @param j    Second index
   * @param n    Number of steps in each directon 
   * @param u    @f$\sigma/\xi@f$ 
   * @param y    Value 
   * @param e    Error
   */
  void FillMG(TMultiGraph* mg,
	      Int_t        i, 
	      Int_t        j, 
	      Int_t        n, 
	      Double_t     u, 
	      Double_t     y, 
	      Double_t     e) 
  { 
    TList* l = mg->GetListOfGraphs();
    TGraphErrors* ge = (l ? static_cast<TGraphErrors*>(l->At(i)) : 0);
    if (!ge) { 
      ge = new TGraphErrors(n);
      ge->SetName(Form("%s%02d", mg->GetName(), i));
      ge->SetTitle(Form(Form("%d", i)));
      ge->SetMarkerColor(AliLandauGaus::GetIColor(i));
      ge->SetLineColor(AliLandauGaus::GetIColor(i));
      ge->SetMarkerStyle(20);
      mg->Add(ge);
    }
    ge->SetPoint(j, u, y);
    ge->SetPointError(j, 0, e);
  }
  /** 
   * Process a single step 
   * 
   * @param mode  Mode 
   * @param i     First index
   * @param j     Second index
   * @param n     Number of steps in each direction 
   * @param bt    Test result 
   */
  virtual void Step(UShort_t     mode, 
		    Int_t        i, 
		    Int_t        j, 
		    Int_t        n, 
		    TestF::Test* bt,
		    TNtuple*     nt) 
  { 
    Test* t = static_cast<Test*>(bt);
    Double_t xi    = t->Xi();
    Double_t sigma = t->Sigma();
    Double_t c     = t->fit->GetParameter(0);
    Double_t p     = t->fit->GetParameter(1);
    Double_t ec    = t->fit->GetParError(0);
    Double_t ep    = t->fit->GetParError(1);
    Int_t    idx   = i * n + j;
    if (mode == 0) {
      c2s->SetPoint(idx, xi, sigma, c);
      p2s->SetPoint(idx, xi, sigma, p);
      c2s->SetPointError(idx, 0, 0, ec);
      p2s->SetPointError(idx, 0, 0, ep);
      FillMG(mcs, i, j, n, sigma/xi, c, ec);
      FillMG(mps, i, j, n, sigma/xi, p, ep);
    }
    else { 
      cs->SetPoint(idx, (mode==1 ? sigma : xi), c);
      ps->SetPoint(idx, (mode==1 ? sigma : xi), p);
      cs->SetPointError(idx, 0, ec);
      ps->SetPointError(idx, 0, ep);
    }
    nt->Fill(1,xi,sigma,t->fit->GetChisquare()/t->fit->GetNDF(),c,ec,p,ep);
  }
  /** 
   * Called after end of loop 
   * 
   * @param mode  Mode of operation 
   * @param can   Canvas 
   * @param nt    NTuple 
   * @param out   Output file 
   */
  void PostLoop(UShort_t mode, TCanvas* can, TNtuple* /*nt*/, TFile* out) 
  { 
    can->Clear();
    if (mode == 0) { 
      can->Divide(2,2);
      can->cd(1); c2s->Draw("TRI2Z");
      can->cd(2); p2s->Draw("TRI2Z");
      can->cd(3); mcs->Draw("APL");
      can->cd(4); mps->Draw("APL");
      c2s->GetHistogram()->SetXTitle("#xi");
      p2s->GetHistogram()->SetXTitle("#xi");
      c2s->GetHistogram()->SetYTitle("#sigma");
      p2s->GetHistogram()->SetYTitle("#sigma");
      c2s->GetHistogram()->SetZTitle("c");
      p2s->GetHistogram()->SetZTitle("p");
      mcs->GetHistogram()->SetXTitle("#sigma/#xi");
      mps->GetHistogram()->SetXTitle("#sigma/#xi");
      mcs->GetHistogram()->SetYTitle("c");
      mps->GetHistogram()->SetYTitle("p");

      out->cd();
      c2s->Write();
      p2s->Write();
      can->Write();
      mcs->Write();
      mps->Write();
    } 
    else { 
      can->Divide(2,1);
      can->cd(1); cs->Draw("APL");
      can->cd(2); ps->Draw("APL");

      out->cd();
      can->Write();
      cs->Write();
      ps->Write();
    }
    PrintCanvas(can, 0, 0);
  }
};
#endif

//====================================================================
#ifdef TEST_FITTER
#include "AliLandauGausFitter.h"
#include <TRandom.h>

struct TestFit : public TestF
{
  struct Test : public TestF::Test 
  {
    TH1* dist;
    TF1* res;
    TH1* pars;
    Bool_t fDoShift;
    Test() : TestF::Test(), dist(0), res(0), pars(0), fDoShift(true) {}
    /** 
     * Re-implementation of TH1::FillRandom to accept a TF1 pointer 
     * 
     * @param h Histogram to fill 
     * @param f Function to use 
     * @param n Number of samples
     */
    static void FillRandom(TH1* h, TF1* f, Int_t n=1000000) { 
      TAxis* xAxis = h->GetXaxis();

      Int_t first  = xAxis->GetFirst();
      Int_t last   = xAxis->GetLast();
      Int_t nbinsx = last-first+1;

      TArrayD integral(nbinsx+1);
      integral[0] = 0;
      for (Int_t i = 1; i <= nbinsx; i++) {
	Double_t fint = f->Integral(xAxis->GetBinLowEdge(i+first-1),
				    xAxis->GetBinUpEdge(i+first-1));
	integral[i] = integral[i-1] + fint;
      }

      //   - Normalize integral to 1
      if (integral[nbinsx] == 0) {
	Error("FillRandom", "Integral = zero"); 
	return;
      }
      for (Int_t i = 1 ; i <= nbinsx; i++)  
	integral[i] /= integral[nbinsx];

      //   --------------Start main loop ntimes
      for (Int_t i = 0; i < n; i++) {
	Double_t r1 = gRandom->Rndm(i);
	Int_t    j  = TMath::BinarySearch(nbinsx,&integral[0],r1);
	Double_t x  = (xAxis->GetBinLowEdge(j+first) +
		       xAxis->GetBinWidth(j+first) * 
		       (r1-integral[j])/(integral[j+1] - integral[j]));
	h->Fill(x);
      }
    }
    /** 
     * Run the actual test
     */
    void DoRun() 
    {
      AliLandauGaus::EnableSigmaShift(1); // Enable sigma-shift
      Double_t xMin = XMin();
      Double_t xMax = XMax();
      Int_t    n    = N();
      src->SetLineColor(kBlue);
      src->SetLineStyle(2);
      dist = new TH1D("dist", "Landau-Gaus distribution",500,xMin,xMax);
      dist->SetXTitle("#Delta");
      dist->SetYTitle("dN/d#Delta");
      dist->SetFillStyle(3001);
      dist->SetFillColor(kMagenta+1);
      dist->Sumw2();
      dist->SetDirectory(0);
      FillRandom(dist, src);
      
      Int_t    peakBin = dist->GetMaximumBin();
      Double_t max     = dist->GetBinContent(peakBin);
      dist->Scale(1/max);

      AliLandauGaus::EnableSigmaShift(fDoShift); // Disable sigma-shift
      AliLandauGausFitter f(xMin, xMax, 10);
      // f.SetDebug(true);
      
      TF1* r = f.FitNParticle(dist, n);
      if (!r) return;
      r->SetName("fit");

      // Get the status for the covariance calculation: 
      // 0: Not calculatuted 
      // 1: Approximate 
      // 2: Forced positive definite 
      // 3: Accurate 
      Int_t covStatus = 
	static_cast<TFitResult*>(f.GetFitResults().At(n-1))->CovMatrixStatus();
      if (covStatus < 2) ::Warning("", "Fit may not be good");
      
      r->SetRange(xMin, xMax);
      res = static_cast<TF1*>(r->Clone());
      res->SetLineColor(kPink+10);
      res->SetLineStyle(1);
      res->SetLineWidth(3);
      dist->GetListOfFunctions()->Add(res);
      
      Double_t min     = 0.1*dist->GetMinimum();
      TF1*     old     = src; // Tuck away src
      src              = res; // Use fit now
      for (Int_t i = 1; i <= n; i++) { 
	Double_t rDelta  = Delta();
	Double_t rxi     = Xi();
	Double_t rsigma  = Sigma();
	TF1*     comp    = Fi(i);
	comp->SetLineWidth(2);
	comp->SetLineStyle(1);

	// Modifies arguments!
	AliLandauGaus::IPars(i, rDelta, rxi, rsigma);
	TLine* l = VerticalLine(1, comp->GetLineColor(), rDelta, min,
				1.1*old->Eval(rDelta));
	dist->GetListOfFunctions()->Add(comp);
	dist->GetListOfFunctions()->Add(l);
      }
      src = old; // restore 
      
    }
    /** 
     * Draw results in pad 
     * 
     * @param p Pad 
     */
    void DrawInPad(TVirtualPad* p)
    {
      // Double_t scale = src->GetMaximum();
      // src->SetParameter(0, 1/scale);
      p->Clear();
      p->SetTopMargin(0.0);
      p->SetRightMargin(0.0);
      p->SetBottomMargin(0.0);
      p->SetLeftMargin(0.0);
      p->SetFillColor(kGray);
      p->Divide(2,1,0.001,0.001);
      TVirtualPad* q = p->cd(1);
      q->SetLogy();
      q->SetTopMargin(0.01);
      q->SetRightMargin(0.02);
      q->SetFillColor(kWhite);
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(11111);
      gStyle->SetOptTitle(false);
      gStyle->SetStatY(1-q->GetTopMargin());

      dist->DrawCopy("hist");
      dist->DrawCopy("e same");
      AliLandauGaus::EnableSigmaShift(1); // Enable sigma-shift
      src->DrawClone("same");
      AliLandauGaus::EnableSigmaShift(fDoShift); // Disable sigma-shift
      if (!res) return;


      q = p->cd(2);
      q->SetTopMargin(0.01);
      q->SetRightMargin(0.02);
      q->SetFillColor(kWhite);
      Int_t nPar = src->GetNpar();
      pars = new TH1D("pars", "#chi^{2}/#nu & #Deltap_{i}/#deltap_{i}", 
		      nPar+1, 0, nPar+1);
      pars->SetDirectory(0);
      pars->SetFillColor(kCyan+2);
      pars->SetFillStyle(3001);
      pars->SetYTitle("#Deltap_{i}/#deltap_{i}");

      TPaveStats* stats = new TPaveStats(0,0,0,0);
      stats->SetX1NDC(gStyle->GetStatX()-gStyle->GetStatW());      
      stats->SetY1NDC(gStyle->GetStatY()-gStyle->GetStatH());
      stats->SetX2NDC(gStyle->GetStatX());
      stats->SetY2NDC(gStyle->GetStatY());
      stats->SetTextFont(42);
      
      stats->SetBorderSize(1);
      stats->SetFillColor(kWhite);
      stats->AddText("");
      pars->GetListOfFunctions()->Add(stats);
      TString fmt(Form("%%s = %%%s", stats->GetFitFormat()));
      for (Int_t i = 0; i < nPar; i++) { 
	Double_t low, high;
	res->GetParLimits(i, low, high);
	Double_t pSrc = src->GetParameter(i);
	Double_t pRes = res->GetParameter(i);
	Double_t eRes = res->GetParError(i);
	Bool_t   fix  = (low == high && ((pRes == 0 && low == 1) || low==pRes));
	Double_t diff = fix ? 0 : TMath::Abs(pSrc-pRes)/eRes;
	pars->GetXaxis()->SetBinLabel(i+1,res->GetParName(i));
	pars->SetBinContent(i+1, diff);
	stats->AddText(Form(fmt.Data(), res->GetParName(i), pSrc));
      }
      pars->GetXaxis()->SetBinLabel(nPar+1,"#chi^{2}/#nu");
      pars->SetBinContent(nPar+1,res->GetChisquare() / res->GetNDF());
      pars->DrawCopy();
      q->Modified();
      q->Update();

      p->Modified();
      p->Update();
      p->cd();
      
      PrintFs(false);
    }
    /** 
     * Print results 
     * 
     * @param full 
     */
    void PrintFs(Bool_t full=false) {
      if (full) {
	Int_t nPar = src->GetNpar();
	Printf("%-2s %-10s | %-8s | %8s %-9s | %9.5s", 
	       "#", "Name", "Source", "Fit", "Error", "delta");
	for (Int_t i = 0; i < nPar; i++) { 
	  Double_t low, high;
	  res->GetParLimits(i, low, high);
	  Double_t pSrc = src->GetParameter(i);
	  Double_t pRes = res->GetParameter(i);
	  Double_t eRes = res->GetParError(i);
	  Bool_t   fix  = (low == high && ((pRes == 0 && low==1) || low==pRes));
	  Double_t diff = fix ? 0 : TMath::Abs(pSrc-pRes)/eRes;
	  Printf("%2d %10s | %8.4f | %8.4f %9.5f | %9.5f %s", 
		 i, src->GetParName(i), pSrc, pRes, eRes, diff, 
		 (fix ? "fixed" : (diff < 2 ? "ok" : "bad")));
	}
	Printf("chi^2/nu = %8.4f/%-3d = %9.5f", 
	       res->GetChisquare(), res->GetNDF(), 
	       res->GetChisquare() / res->GetNDF());
      }
      else 
	Printf("%8.4f", res->GetChisquare() / res->GetNDF());
    }
    /** 
     * Write results to disk 
     * 
     * @param d Directory to write to 
     */
    virtual void WriteOut(TDirectory* d) 
    {
      d->cd();
      src->Write();
      dist->Write();
      if (res) res->Write();
      if (pars) pars->Write();
    }
    
  };
  //__________________________________________________________________
  TestF::Test* MakeTest() { return new Test(); } 
  const char* Prefix() const { return "fit"; }
  void PreLoop(UShort_t,Int_t) {
    Printf("%-8s | %-8s | %-8s | %-8s ",  "Delta_p", "xi", "sigma", "chi^2/nu");
  }
  void Step(UShort_t,Int_t,Int_t,Int_t,TestF::Test*,TNtuple*) {}
  void PostLoop(UShort_t,TCanvas* c,TNtuple*,TFile*) { PrintCanvas(c,0,0); }
};
#endif
// 
// EOF
//

  
