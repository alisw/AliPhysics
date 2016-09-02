/**
 * @file   ToyModeldNdTheta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Sep  1 10:58:01 2016
 * 
 * @brief A class to analyse output of the toy model for
 * @f$\mathrm{d}N/\mathrm{d}\vartheta@f$. 
 * 
 * @ingroup pwglf_forward_tracklets_toy
 */
#ifndef __CINT__
# include <TFile.h>
# include <TNtuple.h>
# include <TH1.h>
# include <THStack.h>
# include <TGraphErrors.h>
# include <TLatex.h>
# include <TCanvas.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TParameter.h>
# include <TError.h>
# include <TMath.h>
#else
class TFile;
class TNtuple;
class TH1;
class TGraph;
class TCanvas;
class THStack;
#endif

/**
 * A class to analyse output of the toy model for
 * @f$\mathrm{d}N/\mathrm{d}\vartheta@f$. 
 * 
 * The procedure requires two input files: 
 * 
 * - A file with data used for corrections - the simulation sample 
 * - A file with data used as the observation - the real sample 
 * 
 * The result is calculated as 
 *
 * @f[ 
 R = \frac{G'}{(1-\beta')M'}(1-\beta)M 
 @f] 
 * where 
 *
 * - @f$ G'@f$ is the distribution of generated tracks in simulation
 * - @f$ M'@f$ is the distribution of the measured tracklets in simulation
 * - @f$ (1-\beta')=C'/M'@f$ is the fraction of good tracklets in simulation 
 * - @f$ C'@f$ is the distribution of fake tracklets in simulation 
 * - @f$ M@f$  is the distribution of the measured tracklets in  real 
 * - @f$\beta=k\beta'@f$ is the estimated fraction of fakes tracklets in real 
 * - @f$ k@f$ is a scaling factor evaluated as 
 * @f[ 
 k = \frac{\int_{l}^{\infty}\mathrm{d}q Q}{\int_{l}^{\infty}\mathrm{d}q Q'}
 @f] 
 * where 
 * 
 * - @f$ Q@f$ is the distribution of tracklet quality (see ToyModel::Tracklet::Quality) in real 
 * - @f$ Q'@f$ is the distribution of tracklet quality in simulation 
 * - @f$ l@f$ is a cut-off 
 *
 * The final result is compared to 
 *
 * - @f$ G@f$ is the distribution of generated tracks in real
 *
 * We may way each and every tracklet by some weight.  See also
 * ToyModeldNdTheta::ThetaP.
 *
 * @ingroup pwglf_forward_tracklets_toy
 */
struct ToyModeldNdTheta
{
  ToyModeldNdTheta() : fNTheta(30), fNQual(1000), fOutput(0) {}
#ifndef __CINT__
  typedef TParameter<double> DParam_t;
#endif
  
  //==================================================================
  /** 
   * @{ 
   * @name Helper functions 
   */
  /** 
   * Open a file in read-only mode 
   * 
   * @param name Name of file to open
   * 
   * @return Pointer to file or null 
   */
  static TFile* OpenFile(const char* name)
  {
    TFile* file = TFile::Open(name, "READ");
    if (!file) {
      ::Error("OpenFile", "Failed to open %s", name);
      return 0;
    }
    return file;
  }
  /** 
   * Get an object from a directory 
   * 
   * @param d       Directory
   * @param name    Name of object 
   * @param cls     Optional class to match to 
   * 
   * @return Pointer to object or null
   */
  TObject* GetO(TDirectory* d, const char* name, TClass* cls=0)
  {
    if (!d) {
      ::Warning("GetO", "No directory to read %s from", name);
      return 0;
    }
    TObject* o = d->Get(name);
    if (!o) {
      ::Error("GetO", "Couldn't find %s in %s", name, d->GetName());
      return 0;
    }
    if (cls && !o->IsA()->InheritsFrom(cls)) {
      ::Error("GetO", "Object %s read from %s is not a %s, but a %s",
	      name, d->GetName(), cls->GetName(), o->ClassName());
      return 0;
    }
    return o;
  }
  /** 
   * Get an ntuple from a directory 
   * 
   * @param d     Directory 
   * @param name  Name of ntuple 
   * 
   * @return Pointer to ntuple or null
   */
  TNtuple* GetN(TDirectory* d, const char* name)
  {
    static_cast<TNtuple*>(GetO(d, name, TNtuple::Class()));
  }
  /** 
   * Get an 1-D histogram from a directory 
   * 
   * @param d     Directory 
   * @param name  Name of 1-D histogram 
   * 
   * @return Pointer to histogram or null
   */
  TH1* GetH1(TDirectory* d, const char* name)
  {
    static_cast<TH1*>(GetO(d, name, TH1::Class()));
  }
  /** 
   * Get an histogram stack from a directory 
   * 
   * @param d     Directory 
   * @param name  Name of histogram stack
   * 
   * @return Pointer to stack or null
   */
  THStack* GetHS(TDirectory* d, const char* name)
  {
    static_cast<THStack*>(GetO(d, name, THStack::Class()));
  }
  /** 
   * Get a graph from a directory 
   * 
   * @param d     Directory 
   * @param name  Name of graph
   * 
   * @return Pointer to graph or null
   */
  TGraph* GetG(TDirectory* d, const char* name)
  {
    static_cast<TGraph*>(GetO(d, name, TGraph::Class()));
  }
  /** 
   * Integrate histogram from low edge to infinity 
   * 
   * @param h    Histogram to integrate 
   * @param low  Low edge 
   * @param err  On return, the error on the integral 
   * 
   * @return The integral 
   */
  static Double_t Integrate(TH1* h, Double_t low, Double_t& err)
  {
    if (!h) {
      ::Warning("Integrate", "Nothing to integrate");
      return 0;
    }
    Int_t start = h->GetXaxis()->FindBin(low);
    Int_t end   = h->GetXaxis()->GetNbins();
    return h->IntegralAndError(start, end, err);
  }
  //____________________________________________________________________
  /** 
   * Take ratio of two numbers with errors 
   * 
   * @param n   Numerator value 
   * @param en  Numerator error 
   * @param d   Denominator value 
   * @param ed  Denominator error 
   * @param er  On return, ratio error 
   * 
   * @return Ratio 
   */
  static Double_t RatioE(Double_t n, Double_t en,
			 Double_t d, Double_t ed,
			 Double_t& er)
  {
    Double_t r = 0;
    er = 0;
    if (TMath::Abs(n) < 1e-16 || TMath::Abs(d) < 1e-16) return 0;
    r  = n/d;
    er = TMath::Sqrt(en*en/n/n + ed*ed/d/d);
    return r;
  }
  //____________________________________________________________________
  /** 
   * Scale a histogram by a constant with errors  
   * 
   * @param h   Histogram to scale 
   * @param x   Scalar
   * @param xe  Error on scalar 
   * 
   * @return The scaled histogram h 
   */
  static TH1* Scale(TH1* h, Double_t x, Double_t xe)
  {
    if (!h) return 0;
    Double_t rr    = xe/x;
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
      Double_t c  = h->GetBinContent(i);
      Double_t e  = h->GetBinError  (i);
      Double_t s  = (c > 0 ? e/c : 0);
      Double_t sc = x * c;
      Double_t se = sc*TMath::Sqrt(s*s+rr*rr);
      // Printf("Setting bin %2d,%2d to %f +/- %f", sc, se);
      h->SetBinContent(i,sc);
      h->SetBinError  (i,se);
    }
    return h;
  }
  /* @} */
  //====================================================================
  /** 
   * @{ 
   * @name Projections of @f$ \vartheta@f$ angles
   */
  /** 
   * Create a histogram of angles 
   * 
   * @param name   Name of histogram 
   * @param title  Title of histogram 
   * @param color  Color used 
   * @param style  Marker style 
   * 
   * @return Pointer to newly created histogram 
   */
  TH1* ThetaH(const char* name,
	      const char* title,
	      Color_t     color,
	      Style_t     style)
  {
    TH1* h = new TH1D(name,title,fNTheta,45,135);
    h->SetFillStyle(0);
    h->SetFillColor(color);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(style);
    h->SetXTitle("#vartheta [Degrees]");
    h->SetYTitle(title);
    h->SetStats(0);
    h->Sumw2();
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);
    return h;
  }
  //____________________________________________________________________
  /** 
   * Project an angle distribution into a histogram.  
   *
   * That is, we project @f$ a 360 / (2\pi)@f$ where @f$ a@f$ is the
   * angle of the stored tracklets.  If the passed ntuple has a branch
   * called @c fake, then we apply the weight according to the prescription
   * @f[ 
   w_{t} = \left\{\begin{array}{cc}w & \mbox{good tracklets}\\ 
   w^2 & \mbox{fake tracklets}\\ \end{array}\right.
   @f] 
   * 
   * We can impose additional cuts on the tracklets by passing a
   * condition in @a cut parameter.
   * 
   * @param t        Ntuple to project from 
   * @param nEv      Number of events
   * @param h        Histogram to project into 
   * @param weight   Weight to be applied to tracklets 
   * @param cut      Cut on tracklets 
   * 
   * @return The histogram project into, or null 
   */
  TH1* ThetaP(TNtuple*    t,
	      Int_t       nEv, 
	      TH1*        h,
	      Double_t    weight=1,
	      const char* cut="qual<1e-4")
  {
    if (!t) {
      Warning("ThetaP", "No tuple to project from");
      return 0;
    }
    if (!h) {
      Warning("ThetaP", "No histogram to project into");
      return 0;
    }
    const char* expr = "angle*TMath::RadToDeg()";
    TString var;  var.Form("%s>>%s", expr, h->GetName());
    TString con;
    if (t->GetBranch("fake")) con.Form("%f*(fake!=0?%f:1)", weight, weight);
    else                      con.Form("%f", weight);
    if (cut && cut[0] != '\0') {
      if (!con.IsNull()) con.Append("*");  
      con.Append(Form("(%s)", cut));
    }
    //Printf("Projecting %s with %s into %s",
    //       var.Data(), con.Data(), h->GetName());
    t->Draw(var, con);

    h->Scale(1./nEv, "width");
  
    return h;
  }
  //____________________________________________________________________
  /** 
   * Create a histogram of angles and project into it 
   * 
   * @param t        Ntuple to project from 
   * @param nEv      Number of events
   * @param name     Name of histogram 
   * @param title    Title of histogram 
   * @param color    Color used 
   * @param style    Marker style 
   * @param weight   Weight to be applied to tracklets 
   * @param cut      Cut on tracklets 
   * 
   * @return Pointer to newly allocated histogram or null
   */
  TH1* ThetaHP(TNtuple*    t,
	       Int_t       nEv, 
	       const char* name,
	       const char* title,
	       Color_t     color,
	       Style_t     style,
	       Double_t    weight=1,
	       const char* cut="qual<1e-4")
  {
    TH1* h = ThetaH(name, title, color, style);
    return ThetaP(t, nEv, h, weight, cut);
  }
  /* @} */
  //====================================================================
  /** 
   * @{ 
   * @name Projections of @f$ \vartheta@f$ angles
   */
  /** 
   * Create a histogram of qualities 
   * 
   * @param name   Name of histogram 
   * @param title  Title of histogram 
   * @param color  Color used 
   * @param style  Marker style 
   * 
   * @return Pointer to newly created histogram 
   */
  TH1* QualH(const char* name,
	     const char* title,
	     Color_t     color,
	     Style_t     style)
  {
    TH1* h = new TH1D(name,title,fNQual,0,0.004);
    // TH1* h = new TH1D("h","H",100, -1, 1);
    // h->SetDirectory(0);
    h->SetFillStyle(0);
    h->SetFillColor(color);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(style);
    h->SetXTitle("Quality");
    h->SetYTitle(title);
    h->SetStats(0);
    h->Sumw2();
    h->GetXaxis()->SetNdivisions(205);
    h->GetYaxis()->SetNdivisions(205);
    return h;
  }
  //____________________________________________________________________
  /** 
   * Project a quality distribution into a histogram 
   * 
   * @param t        Ntuple to project from 
   * @param nEv      Number of events
   * @param h        Histogram to project into 
   * @param weight   Weight to be applied to tracklets 
   * @param cut      Cut on tracklets 
   * 
   * @return The histogram project into, or null 
   */
  TH1* QualP(TNtuple*    t,
	     Int_t       nEv, 
	     TH1*        h,
	     Double_t    weight=1,
	     const char* cut="")
  {
    if (!t) {
      Warning("QualP", "No tuple to project from");
      return 0;
    }
    if (!h) {
      Warning("QualP", "No histogram to project into");
      return 0;
    }
    const char* expr = "qual";
    TString var;  var.Form("%s>>%s", expr, h->GetName());
    TString con;  con.Form("%f", weight);
    if (cut && cut[0] != '\0') con.Append(Form("*(%s)", cut));
    t->Draw(var, con);

    h->Scale(1./nEv, "width");
  
    return h;
  }
  //____________________________________________________________________
  /** 
   * Create a histogram of qualities and project into it 
   * 
   * @param t        Ntuple to project from 
   * @param nEv      Number of events
   * @param name     Name of histogram 
   * @param title    Title of histogram 
   * @param color    Color used 
   * @param style    Marker style 
   * @param weight   Weight to be applied to tracklets 
   * @param cut      Cut on tracklets 
   * 
   * @return Pointer to newly allocated histogram or null
   */
  TH1* QualHP(TNtuple*    t,
	      Int_t       nEv, 
	      const char* name,
	      const char* title,
	      Color_t     color,
	      Style_t     style,
	      Double_t    weight=1,
	      const char* cut="")
  {
    TH1* h = QualH(name, title, color, style);
    return QualP(t, nEv, h, weight, cut);
  }
  /* @} */

  //====================================================================
  /** 
   * Run the analysis 
   * 
   * @param realFileName The measured tracklets 
   * @param simFileName  The tracklets used for corrections
   * @param weight       Reweighing factor 
   */
  void Run(const char* realFileName="real.root",
	   const char* simFileName="reduced.root",
	   Double_t    weight=1)
  {
    // Constants for colors 
    const Color_t cM = kBlue+1;
    const Color_t cG = kBlack;
    const Color_t cC = kMagenta+1;
    const Color_t cS = kGreen+1;
    const Color_t cB = kPink+1;
    const Color_t cA = kCyan+1;
    const Color_t cR = kRed+1;

    // Cuts for fakes, signals, and both 
    const char*   fC = "fake!=0";
    const char*   sC = "qual<1e-4";
    TString       bC; bC.Form("%s&&%s", fC, sC);

    // Square weight 
    Double_t w2      = weight*weight;

    // Open files 
    TFile* realFile = OpenFile(realFileName);
    TFile* simFile  = OpenFile(simFileName);
    if (!realFile || !simFile) return;

    // Get the tuples 
    TNtuple* realTkl = GetN(realFile, "record");
    TNtuple* simTkl  = GetN(simFile,  "record");
    TNtuple* realTrk = GetN(realFile, "generated");
    TNtuple* simTrk  = GetN(simFile,  "generated");
    TNtuple* realSta = GetN(realFile, "stat");
    TNtuple* simSta  = GetN(simFile,  "stat");
    Int_t    nR      = realSta->GetEntries();
    Int_t    nS      = simSta ->GetEntries();

    // Output file
    if (fOutput) {
      fOutput->Close();
    }
    fOutput  = TFile::Open("result.root", "RECREATE");

    // Read data into histograms
    TH1*     realM   = ThetaHP(realTkl,nR,"realM","M", cM,21,1,     sC);
    TH1*     realG   = ThetaHP(realTrk,nR,"realG","G", cG,20,1,     "");
    TH1*     simM    = ThetaHP(simTkl, nS,"simM", "M'",cM,25,weight,sC);
    TH1*     simG    = ThetaHP(simTrk, nS,"simG", "G'",cG,24,weight,"");
    TH1*     simC    = ThetaHP(simTkl, nS,"simC", "C'",cC,27,weight,bC);

    // Read quality distributions 
    TH1*     realQM  = QualHP(realTkl,nR,"realQM","#Delta_{M}", cM,21);  
    TH1*     simQM   = QualHP(simTkl, nS,"simQM", "#Delta_{M'}",cM,25,weight);
    TH1*     simQC   = QualHP(simTkl, nS,"simQC","#Delta_{C'}",cC,32,weight,fC);

    // Declare derived histograms
    TH1*     simS    = ThetaH("simS",   "S'",      cS, 26);
    TH1*     simB    = ThetaH("simB",   "#beta'",  cB, 26);
    TH1*     sim1B   = ThetaH("sim1B",  "1-#beta'",cB, 26);
    TH1*     real1B  = ThetaH("real1B", "1-#beta", cB, 22);
    TH1*     simA    = ThetaH("simA",   "#alpha'", cA, 26);
    TH1*     realS   = ThetaH("realS",  "S",       cS, 22);
    TH1*     realC   = ThetaH("realC",  "C",       cC, 23);
    TH1*     realR   = ThetaH("realR",  "R",       cR, 20);
    TH1*     closure = ThetaH("closure","R/G",     cR, 20);
    // Make a histogram of ones 
    TH1*     one     = ThetaH("one",    "1",       kBlack, 0);
    for (Int_t i = 1; i <= one->GetNbinsX(); i++) {
      one->SetBinContent(i, 1);
      one->SetBinError(i, 0);
    }
    // Scale quality histograms
    Double_t realQIE, realQI = Integrate(realQM, 0, realQIE);
    Double_t simQIE,  simQI  = Integrate(simQM,  0, simQIE);  
    Scale(realQM, 1/realQI, 1/realQIE/realQIE);
    Scale(simQM,  1/simQI,  1/simQIE/simQIE);
    Scale(simQC,  1/simQI,  1/simQIE/simQIE);
  
    // Calculate scaling factor k
    realQI         = Integrate(realQM, 1e-4, realQIE);
    simQI          = Integrate(simQM,  1e-4, simQIE);  
    Double_t ek, k = RatioE(realQI, realQIE, simQI, simQIE, ek);
    Printf("Scalar is equal to %f +/- %f", k, ek);

    // Calculate derived histograms 
    simB   ->Divide(simC, simM);      // Calculate simulated beta   
    simS   ->Add(simM, simC, 1, -1);  // Calculate simulated signal   
    sim1B  ->Add(one, simB, 1, -1);   // Calculate simulated 1-beta   
    real1B ->Add(one, simB, 1, -k);   // Calculate real 1-beta     
    simA   ->Divide(simG, simS);      // Calculate simulated alpha
    realC  ->Add(simC, k);            // Calculate estimated background
    realS  ->Multiply(real1B, realM); // Calculate real signal   
    realR  ->Multiply(simA, realS);   // Calculate the result
    closure->Divide(realR, realG);    // Closure test
    closure->SetTitle("Closure");

    // Create stack of parts 
    THStack* parts   = new THStack("parts","");
    parts->Add(realG);
    parts->Add(simG);
    parts->Add(realM);
    parts->Add(simM);
    parts->Add(realC);
    parts->Add(simC);
    parts->Add(realS);
    parts->Add(simS);
    parts->Add(realR);

    // Create stack of quality histograms 
    THStack* deltas = new THStack("deltas", "");
    deltas->Add(realQM);
    deltas->Add(simQM);
    deltas->Add(simQC);

    // Create stack of other histograms 
    THStack* other = new THStack("other", "");
    other->Add(simA);
    other->Add(simB);
    other->Add(sim1B);
    other->Add(real1B);
    TGraphErrors* g = new TGraphErrors;
    g->SetName("k");
    g->SetTitle(Form("k=%5.3f#pm%5.3f",k,ek));
    g->SetFillStyle(3002);
    g->SetFillColor(cC);
    g->SetLineColor(cC);
    g->SetLineStyle(2);
    g->SetLineWidth(2);  
    g->SetPoint(0,simA->GetXaxis()->GetXmin(), k);
    g->SetPoint(1,simA->GetXaxis()->GetXmax(), k);
    g->SetPointError(0,0, ek);
    g->SetPointError(1,0, ek);

    fOutput->cd();
    TNamed*    realName = new TNamed("realFile",realFileName);
    TNamed*    simName  = new TNamed("simFile", simFileName);
    DParam_t*  pw       = new DParam_t("weight", weight);
    realName->Write();
    simName->Write();
    pw->Write();
    parts->Write();
    other->Write();
    deltas->Write();
    g->Write();
    fOutput->Write();
  }
  void Visualize()
  {
    // Make a canvas 
    TCanvas* c = new TCanvas("c","C", 800, 900);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.01);
    c->cd();

    TObject*  realFile = GetO(fOutput, "realFile");
    TObject*  simFile  = GetO(fOutput, "simFile");
    DParam_t* weight   = static_cast<DParam_t*>(GetO(fOutput, "weight"));
    TLatex* title = new TLatex(0.5,.99,
			       Form("Data: %s  Corr: %s  Weight: %5.3f",
				    realFile->GetTitle(),
				    simFile->GetTitle(),
				    weight->GetVal()));
    title->SetNDC();
    title->SetTextAlign(23);
    title->SetTextFont(42);
    title->SetTextSize(0.03);
    title->Draw();
  
    TPad* sub = new TPad("sub","sub", 0, 0, 1, .95);
    sub->SetTopMargin(0.10);
    sub->SetRightMargin(0.01);
    sub->Draw();
    sub->cd();

    TPad* left = new TPad("left","left",0,0,.5,1);
    left->SetTopMargin(0.10);
    left->SetRightMargin(0.01);
    left->Draw();
    left->cd();
    left->Divide(1,2,0,0);

    sub->cd();
    TPad* right = new TPad("right","right",.5,0,1,1);
    right->SetTopMargin(0.10);
    right->SetRightMargin(0.01);
    right->Draw();
    right->cd();
    right->Divide(1,2);
    
    TVirtualPad* p = 0;
    TLegend* l     = 0;
    // Draw our stack of parts 
    p = left->cd(1);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetBottomMargin(0);
    p->SetTicks();
    THStack* parts = GetHS(fOutput, "parts");
    if (parts) {
      parts->SetMaximum(1.5*parts->GetMaximum("nostack"));
      parts->Draw("nostack");
      parts->GetHistogram()->SetXTitle("#vartheta [Degrees]");
      parts->GetHistogram()->SetYTitle("R,G',G,M,M',C,C'");
      l = p->BuildLegend(p->GetLeftMargin(), .7,
			 1-p->GetRightMargin(), 1-p->GetTopMargin());
      l->SetNColumns(2);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      p->Modified();
    }

    // Draw our stack of qualities 
    p = right->cd(1);
    p->SetLogy();
    p->SetLogx();
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    THStack* deltas = GetHS(fOutput,"deltas");
    if (deltas) { 
      deltas->Draw("nostack");
      deltas->GetHistogram()->SetXTitle("#Delta=(#sigma_{1}^{2}+#sigma_{2}^{2}+"
					"#delta_{z}^{2})");
      deltas->GetHistogram()->SetYTitle("");
      l = p->BuildLegend(.7, .7,
			 1-p->GetRightMargin(), 1-p->GetTopMargin());
      l->SetNColumns(1);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      p->Modified();
    }
  

    // Draw closure test 
    p = left->cd(2);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetTicks();
    TH1* closure = GetH1(fOutput, "closure");
    if (closure) closure->Draw();

    // Draw other histograms 
    p = right->cd(2);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetTicks();
    THStack* other = GetHS(fOutput, "other");
    TGraph*  g     = GetG(fOutput, "k");
    if (other) {
      other->SetMaximum(1.5*other->GetMaximum("nostack"));
      other->Draw("nostack");
      other->GetHistogram()->GetListOfFunctions()->Add(g, "3");
      other->GetHistogram()->SetXTitle("#vartheta [Degrees]");
      other->GetHistogram()->SetYTitle("#alpha',#beta',(1-#beta'),"
				       "(1-k#beta'),k");
      if (g) g->Draw("3");
      l = p->BuildLegend(p->GetLeftMargin(), .7,
			 1-p->GetRightMargin(), 1-p->GetTopMargin());
      l->SetNColumns(2);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      p->Modified();
    }

    c->Modified();
    c->Update();
    c->cd();
  }

  Int_t  fNTheta; // Number of @f$\vartheta@f$ bins
  Int_t  fNQual;  // Number of quality bins
  TFile* fOutput; // Output file 
};

//
// EOF
//
