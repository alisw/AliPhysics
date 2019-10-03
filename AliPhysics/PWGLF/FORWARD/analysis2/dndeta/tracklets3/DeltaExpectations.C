/**
 * @file   DeltaExpectations.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 17:07:01 2016
 * 
 * @brief  Do back of the envelope calculation of effect of reweighing
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */

#ifndef __CINT__
# include <vector>
# include <TFile.h>
# include <TError.h>
# include <TCollection.h>
# include <TString.h>
# include <TH1.h>
# include <TH2.h>
# include <THStack.h>
# include <TClass.h>
# include <TBrowser.h>
# include <TCanvas.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TStyle.h>
#else
class TFile;
class TClass;
class TH1;
class TH2;
class TDirectory;
class THStack;
class TCollection;
class TBrowser;
class TCanvas;
class TVirtualPad;
#endif

/**
 * Print format 
 * 
 * @ingroup pwglf_forward_tracklets
 */
const char* parFmt = "%-10s | %6.3f |";
const char* rowFmt = " %5.2f | %5.2f | %5.2f | %5.2f |";
const char* eRowFmt = " %5.2f+/-%5.2f | %5.2f+/-%5.2f | "
  "%5.2f+/-%5.2f | %5.2f+/-%5.2f |";
struct Utilities
{
  /** 
   * Check class of object
   * 
   * @param o Object to check
   * @param c Optional class 
   * 
   * @return o if the class match, or null
   */
  static TObject* ChkC(TObject* o, TClass* c)
  {
    if (!o) return 0;
    if (!c) return o;
    if (!o->IsA()->InheritsFrom(c)) {
      ::Warning("ChkC", "Object %s is not a %s but a %s",
		o->GetName(), c->GetName(), o->ClassName());
      return 0;
    }
    return o;
  }
  /** 
   * Check that we got an object from a source 
   * 
   * @param src  Source 
   * @param o    Object
   * @param name Object name 
   * 
   * @return Object or null
   */
  static TObject* ChkO(TObject* src, TObject* o, const char* name)
  {
    if (!o) {
      ::Warning("ChkO", "Object %s not found in %s", name, src->GetName());
      src->ls();
      return 0;
    }
    return o;
  }
  /** 
   * Get an object from a source 
   * 
   * @param d     Directory 
   * @param name  Name 
   * @param cls   Optional class to check for 
   * 
   * @return Found object or null
   */
  static TObject* GetO(TDirectory* d, const char* name, TClass* cls)
  {
    if (!d) {
      ::Warning("GetO", "No directory passed for %s", name);
      return 0;
    }
    TObject* o = d->Get(name);
    return ChkC(ChkO(d, o, name), cls);    
  }
  /** 
   * Get an object from a source 
   * 
   * @param d     Collection 
   * @param name  Name 
   * @param cls   Optional class to check for 
   * 
   * @return Found object or null
   */
  static TObject* GetO(TCollection* d, const char* name, TClass* cls)
  {
    if (!d) {
      ::Warning("GetO", "No collection passed for %s", name);
      return 0;
    }
    TObject* o = d->FindObject(name);
    return ChkC(ChkO(d, o, name), cls);    
  }
  /** 
   * Get a collection from a source 
   * 
   * @param d    Source 
   * @param name Name of collection
   * 
   * @return Collection or null
   */
  static TCollection* GetC(TDirectory* d, const char* name)
  {
    return static_cast<TCollection*>(GetO(d, name, TCollection::Class()));
  }
  /** 
   * Get a collection from a source 
   * 
   * @param d    Source 
   * @param name Name of collection
   * 
   * @return Collection or null
   */
  static TCollection* GetC(TCollection* d, const char* name)
  {
    return static_cast<TCollection*>(GetO(d, name, TCollection::Class()));
  }
  /** 
   * Get a histogram from a source 
   * 
   * @param d    Source 
   * @param name Name of histogram
   * 
   * @return Histogram or null
   */
  static TH1* GetH1(TCollection* d, const char* name)
  {
    return static_cast<TH1*>(GetO(d, name, TH1::Class()));
  }
  /** 
   * Get a histogram from a source 
   * 
   * @param d    Source 
   * @param name Name of histogram
   * 
   * @return Histogram or null
   */
  static TH2* GetH2(TCollection* d, const char* name)
  {
    return static_cast<TH2*>(GetO(d, name, TH2::Class()));
  }
  /** 
   * Get a histogram stack from a source 
   * 
   * @param d    Source 
   * @param name Name of histogram stack
   * 
   * @return Histogram stack or null
   */
  static THStack* GetHS(TCollection* d, const char* name)
  {
    return static_cast<THStack*>(GetO(d, name, THStack::Class()));
  }
};

/**
 * Class to do the calculation 
 * 
 * @ingroup pwglf_forward_tracklets
 */
struct DeltaCalculations : public Utilities
{
  enum EColumn {
    kPs,
    kSs,
    kCs,
    kCb
  };
  enum ERow {
    kReduced,
    kExpectation,
    kReweighed,
    kAsIs
  };
  enum EKind {
    kK0S     = 1,
    kKpm,
    kLambda,
    kXi,
    kSigma,
    kOther
  };
  struct Row
  {
    Double_t Ps;
    Double_t Ss;
    Double_t Cs;
    Double_t Cb;
    Double_t ePs;
    Double_t eSs;
    Double_t eCs;
    Double_t eCb;
    /** 
     * Get a quantity
     * 
     * @param i number of parameter 
     * 
     * @return Parameter
     */
    Double_t Get(EColumn i) const
    {
      return (i == kPs ? Ps :
	      i == kSs ? Ss :
	      i == kCs ? Cs :
	      i == kCb ? Cb : 0);
    }
    /** 
     * Get the error on a quantity
     * 
     * @param i number of parameter 
     * 
     * @return Parameter error 
     */
    Double_t GetE(EColumn i) const
    {
      return (i == kPs ? ePs :
	      i == kSs ? eSs :
	      i == kCs ? eCs :
	      i == kCb ? eCb : 0);
    }
    /** 
     * Scale a row 
     * 
     * @param i Row number 
     * @param s Scalar 
     */
    void Scale(EColumn i, Double_t s)
    {
      switch (i) {
      case kPs: Ps *= s; ePs *= s; break;
      case kSs: Ss *= s; eSs *= s; break;
      case kCs: Cs *= s; eCs *= s; break;
      case kCb: Cb *= s; eCb *= s; break;
      }
    }
    /** 
     * Print a row 
     * 
     */
    void Print(Bool_t err=false) const
    {
      if (!err) printf(rowFmt, Ps, Ss, Cs, Cb);
      else      printf(eRowFmt, Ps, ePs, Ss, eSs, Cs, eCs, Cb, eCb);
    }
  };
  /** 
   * Per particle information
   * 
   * @ingroup pwglf_forward_tracklets
   */
  struct Particle
  {
    TString  n;
    Double_t w;
    Row      reduced;
    Row      expected;
    Row      reweighed;
    Row      asis;
    /** 
     * Get the row corresponding to @f$ P^s, S^s, C^s, C^b@f$
     * 
     * @param i Row number 
     * 
     * @return Reference to row 
     */
    Row& Get(ERow i)
    {
      return (i == kReduced     ? reduced :
	      i == kExpectation ? expected :
	      i == kReweighed   ? reweighed :
	      asis);
    }
    /** 
     * Print the particle 
     * 
     */
    void Print(Bool_t err=false)
    {
      printf(parFmt, n.Data(), w);
      reduced.  Print(err);
      expected. Print(err);
      reweighed.Print(err);
      asis.     Print(err);
      printf("\n");
    }
    /** 
     * Read the particle details from passed collections
     * 
     * @param row     Which row (i.e., particle)
     * @param kind    Which column (i.e., @f$ P^s, S^s, C^s, C^b@f$)
     * @param bin     The top-level (centrality) collection 
     * @param sub     Which sub-part to get (i.e., @c mid or @c fwd)
     * @param oth     The other sub-part 
     * @param abso    If true, store absolute numbers 
     * 
     * @return Ourselves 
     */
    Particle* Read(ERow           row,
		   EKind          kind,
		   TCollection*   bin,
		   const TString& sub,
		   const TString& oth,
		   Bool_t         abso=true)
    {
      TString      pn  = "primaries";
      TString      sn  = "secondaries";
      TString      cn  = "combinatorics";
      TCollection* pc  = GetC(GetC(GetC(bin,pn),"specie"),sub);
      TCollection* sc  = GetC(GetC(GetC(bin,sn),"specie"),sub);
      TCollection* cc  = GetC(GetC(GetC(bin,cn),"specie"),sub);
      TCollection* po  = GetC(GetC(GetC(bin,pn),"specie"),oth);
      TCollection* so  = GetC(GetC(GetC(bin,sn),"specie"),oth);
      TCollection* co  = GetC(GetC(GetC(bin,cn),"specie"),oth);
      THStack*     hp  = GetHS(pc,"rows");
      THStack*     hs  = GetHS(sc,"rows");
      THStack*     hc  = GetHS(cc,"rows");
      TH1*         rps = GetH1(hp->GetHists(), "rowSig");
      TH1*         rss = GetH1(hs->GetHists(), "rowSig");
      TH1*         rcs = GetH1(hc->GetHists(), "rowSig");
      TH1*         rcb = GetH1(hc->GetHists(), "rowBg");
      TH1*         pi  = GetH1(pc, "totalIntegrals");
      TH1*         si  = GetH1(sc, "totalIntegrals");
      TH1*         ci  = GetH1(cc, "totalIntegrals");
      TH1*         pj  = GetH1(po, "totalIntegrals");
      TH1*         sj  = GetH1(so, "totalIntegrals");
      TH1*         cj  = GetH1(co, "totalIntegrals");
      Double_t     ww  = (row == kExpectation ? w : 1);
      Int_t        bn  = kind;
      Row&         ret = Get(row);
      if (n.IsNull()) n = rps->GetXaxis()->GetBinLabel(bn);
      if (abso) {
	ret.Ps  = ww*pi->GetBinContent(2)*rps->GetBinContent(bn);
	ret.Ss  = ww*si->GetBinContent(2)*rss->GetBinContent(bn);
	ret.Cs  = ww*ci->GetBinContent(2)*rcs->GetBinContent(bn);
	ret.Cb  = ww*ci->GetBinContent(3)*rcb->GetBinContent(bn);
	ret.ePs = ww*pi->GetBinContent(2)*rps->GetBinError  (bn);
	ret.eSs = ww*si->GetBinContent(2)*rss->GetBinError  (bn);
	ret.eCs = ww*ci->GetBinContent(2)*rcs->GetBinError  (bn);
	ret.eCb = ww*ci->GetBinContent(3)*rcb->GetBinError  (bn);
      }
      else {
	ret.Ps  = ww*rps->GetBinContent(bn);
	ret.Ss  = ww*rss->GetBinContent(bn);
	ret.Cs  = ww*rcs->GetBinContent(bn);
	ret.Cb  = ww*rcb->GetBinContent(bn);
	ret.ePs = ww*rps->GetBinError  (bn);
	ret.eSs = ww*rss->GetBinError  (bn);
	ret.eCs = ww*rcs->GetBinError  (bn);
	ret.eCb = ww*rcb->GetBinError  (bn);
      }
      return this;
    }
  };
#ifndef __CINT__ 
  typedef std::vector<Particle*> ParticleList;
#else
  typedef std::vector<void*> ParticleList;
#endif 
  ParticleList particles;
  Double_t     fP;  // Ratio of primaries to all tracklets 	       
  Double_t     fS;  // Ratio of secondaries to all tracklets       
  Double_t     fC;  // Ratio of combinatorics to all tracklets     
  Double_t     fPs; // Ratio of signal primaries to all tracklets 	       
  Double_t     fSs; // Ratio of signal secondaries to all tracklets   
  Double_t     fCs; // Ratio of signal combinatorics to all tracklets
  /** 
   * Print out 
   * 
   */
  void Print(Bool_t err=false)
  {
    printf("%-18s |", "Sample");
    if (err) {
      printf(" %-29s |", "Reduced");
      printf(" %-29s |", "Expected");
      printf(" %-29s |", "Reweighed");
      printf(" %-29s |", "As-is");
    }
    else {
      printf(" %-61s |", "Reduced");
      printf(" %-61s |", "Expected");
      printf(" %-61s |", "Reweighed");
      printf(" %-61s |", "As-is");
    }
    printf("\n");
    TString head(parFmt);
    for (Int_t i = 0; i < 4; i++) head.Append(rowFmt);
    head.ReplaceAll("%6.3f", "%6s");
    head.ReplaceAll("%5.2f", err ? "%13s" : "%5s");
    TString lne(head);
    lne.ReplaceAll("|", "+");
    lne.ReplaceAll("%5s", "-----");
    lne.ReplaceAll("%6s", "------");
    lne.ReplaceAll("%10s", "----------");
    lne.ReplaceAll("%13s", "-------------");
    lne.ReplaceAll(" ", "-");
    Printf(head.Data(), "Type", "W",
	   "Ps", "Ss", "Cs", "Cb",
	   "Ps", "Ss", "Cs", "Cb",
	   "Ps", "Ss", "Cs", "Cb",
	   "Ps", "Ss", "Cs", "Cb");
    Printf(lne);
    for (ParticleList::const_iterator i = particles.begin();
	 i != particles.end(); ++i)
      (*i)->Print(err);
    Printf(lne);
  }
  /** 
   * Scale the expected numbers to 100% 
   * 
   */
  void ScaleExpected()
  {
    for (Int_t i = 0; i < 4; i++) {
      EColumn c = (i == 0 ? kPs :
		   i == 1 ? kSs :
		   i == 2 ? kCs : kCb);
      Double_t sum = 0;
      for (ParticleList::const_iterator j = particles.begin();
	   j != particles.end(); ++j)
	sum += (*j)->expected.Get(c);
      for (ParticleList::const_iterator j = particles.begin();
	   j != particles.end(); ++j)
	(*j)->expected.Scale(c,100/sum);
    }
  }
      
      
  /** 
   * Create a particle 
   * 
   * @param w           Weight of particle 
   * @param column      Column in histograms to read off 
   * @param reduced     The reduced set 
   * @param reweighed   The reweighed set
   * @param asis        The asis set 
   * @param bin         The bin name 
   * @param sub         The sub name 
   * @param oth         The other sub-name 
   * @param abso        If true, get absolute numbers 
   * 
   * @return Create Particle object
   */
  Particle* Create(Double_t       w,
		   EKind          kind,		   
		   TCollection*   reduced,
		   TCollection*   reweighed,
		   TCollection*   asis,
		   const TString& sub,
		   const TString& oth,
		   Bool_t         abso)
  {
    Particle* ret = new Particle;
    ret->w = w;
    ret->Read(kReduced,     kind, reduced,   sub, oth, abso);
    ret->Read(kExpectation, kind, reduced,   sub, oth, abso);
    ret->Read(kReweighed,   kind, reweighed, sub, oth, abso);
    ret->Read(kAsIs,        kind, asis,      sub, oth, abso);
    particles.push_back(ret);
    return ret;
  }
  /** 
   * Make a stack 
   * 
   * @param column Which column (Ps,Ss,Cs, or Cb) 
   * @param name   Name of stack 
   * @param title  Title of stack 
   * 
   * @return The created stack 
   */
  THStack* MakeStack(EColumn     column,
		     const char* name,
		     const char* title,
		     Int_t       mode)
  {
    const char*   names[]  = { "red", "exp", "rew", "asi" };
    const char*   titles[] = { "Reduced", "Expected", "Reweighed", "As-is" };
    const Color_t colors[] = { kRed+1, kMagenta+1, kBlue+1, kGreen+1 };
    Int_t         nPart    = particles.size();
    Double_t      bW       = 0.15;
    Double_t      bS       = 0.05;
    Double_t      bT       = 4*bW+3*bS;
    Double_t      bO       = (1-bT)/2;
    THStack* stack = new THStack(name, title);
    for (Int_t i = 0; i < 4; i++) {
      ERow q = (i == 0 ? kReduced     :
		i == 1 ? kExpectation :
		i == 2 ? kReweighed   : kAsIs);
      TH1* h = new TH1D(names[i],titles[i],nPart,0,nPart);
      h->SetFillColor(colors[i]);
      h->SetLineColor(colors[i]);
      h->SetMarkerColor(colors[i]);      
      h->SetBarOffset(bO+i*(bW+bS));
      h->SetMarkerSize(mode == 2 ? 3.5 : 2);
      h->SetBarWidth (bW);
      h->SetDirectory(0);
      h->SetYTitle("Fraction [%]");
      stack->Add(h, "bar0 text90");
      
      for (Int_t j = 0; j < nPart; j++) {
	Particle* p = particles[j];
	Row&      r = p->Get(q);
	h->GetXaxis()->SetBinLabel(j+1,p->n);
	h->SetBinContent          (j+1,r.Get (column));
	h->SetBinError            (j+1,r.GetE(column));
      }
    }
    return stack;
  }
  /** 
   * Draw a stack 
   * 
   * @param mother Mother pad  
   * @param column Column being drawn
   * @param stack  Stack to draw 
   * @param logy   If true, plot on logarithmic scale 
   * @param abso   If true, assume absolute numbers
   */
  void DrawStack(TVirtualPad* mother,
		 EColumn      column,
		 THStack*     stack,
		 Bool_t       logy,
		 Bool_t       abso,
		 Int_t        mode)
  {
    TVirtualPad* p = mother->cd(column+1);
    p->SetGridx();
    p->SetGridy();
    p->SetTicks();
    if (logy) p->SetLogy();
    if (mode==2 || (p->GetNumber() % (mode==1 ? 4 : 2)) == 0)
      p->SetRightMargin(0.01);
    if (!abso) {
      stack->SetMinimum(0.003);
      stack->SetMaximum(logy ? (mode == 2 ? 900 : 400) : 110);
    }
    gStyle->SetTitleFontSize(0.02/p->GetHNDC());
    stack->Draw("nostack hist");
    stack->GetHistogram()->GetYaxis()->SetNdivisions(210);
    stack->GetHistogram()->GetYaxis()->SetLabelSize(0.03/p->GetHNDC());
    stack->GetHistogram()->GetYaxis()->SetTitleSize(0.03/p->GetHNDC());
    stack->GetHistogram()->GetXaxis()->SetLabelSize(0.03/p->GetHNDC());
    TString tit = abso ? "#it{N}_{tracklet}" : "Fraction [%]";
    if ((mode == 2 || mode == 0) && p->GetNumber() != 1) tit="";
    stack->GetHistogram()->GetYaxis()->SetTitle(tit);
    stack->GetHistogram()->GetYaxis()->SetTitleOffset(mode == 2 ? 0.5 : 1);
  }
      
  /** 
   * Run it 
   * 
   * @param fileName Input file 
   * @param binName  Bin name 
   * @param mid      If true, for |eta|<1, otherwise |eta|>1
   */
  void Run(Double_t c1=0, Double_t c2=0, Bool_t mid=true)
  {
    Bool_t abso = false;
    Bool_t logy = true;
    Int_t  mode = 2; // 0: square, 1: landscape, 2: portrait
    
    TString bin("all");
    if (c2 > c1)
      bin.Form("cent%03dd%02d_%03dd%02d",
	       Int_t(c1), Int_t(100*c1) % 100, 
	       Int_t(c2), Int_t(100*c2) % 100);
    
    TFile*  redFile = TFile::Open("stk.root", "READ");
    TFile*  rweFile = TFile::Open("wstk.root", "READ");
    TFile*  trhFile = TFile::Open("hijing.root", "READ");
    if (!redFile || !rweFile || !trhFile) return;
    TString      sub    = (mid ? "mid" : "fwd");
    TString      oth    = (mid ? "fwd" : "mid");
    TCollection* redTop = GetC(GetC(redFile,"MidRapidityMCResults"),bin.Data());
    TCollection* rweTop = GetC(GetC(rweFile,"MidRapidityMCResults"),bin.Data());
    TCollection* trhTop = GetC(GetC(trhFile,"MidRapidityMCResults"),bin.Data());
    
    Create(1.52233, kK0S,    redTop, rweTop, trhTop, sub, oth, abso);
    Create(1.41178, kKpm,    redTop, rweTop, trhTop, sub, oth, abso);
    Create(2.75002, kLambda, redTop, rweTop, trhTop, sub, oth, abso);
    Create(3.24110, kXi,     redTop, rweTop, trhTop, sub, oth, abso);
    Create(2.75002, kSigma,  redTop, rweTop, trhTop, sub, oth, abso);
    Create(1,       kOther,  redTop, rweTop, trhTop, sub, oth, abso);

    if (!abso) ScaleExpected();
    Print(true);

    THStack* Ps = MakeStack(kPs, "Ps", "Signal primaries",   mode);
    THStack* Ss = MakeStack(kSs, "Ss", "Signal secondaries", mode);
    THStack* Cs = MakeStack(kCs, "Cs", "Signal fake",        mode);
    THStack* Cb = MakeStack(kCb, "Cb", "Background fake",    mode);

    gStyle->SetPaintTextFormat("5.2f%%");

    Int_t    cw = (mode == 1 ? 1600 : mode == 2 ? 800    : 1000);
    Int_t    ch = (mode == 1 ? cw/2 : mode == 2 ? 1.2*cw : cw);
    TCanvas* c  = new TCanvas(Form("deltaExpectation_%s_%s",
				   bin.Data(), mid?"mid":"fwd"),
			      "DeltaCanvas",cw,ch);
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.01);
    c->SetLeftMargin(mode == 1 ? 0.09 : 0.13);
    switch (mode) {
    case 1:  c->Divide(4,1,0,0); break;
    case 2:  c->Divide(1,4,0,0); break;
    default: c->Divide(2,2,0,0); break;
    }
    
    DrawStack(c, kPs, Ps, logy, abso, mode);
    DrawStack(c, kSs, Ss, logy, abso, mode);
    DrawStack(c, kCs, Cs, logy, abso, mode);
    DrawStack(c, kCb, Cb, logy, abso, mode);

    TVirtualPad* p = c->cd(1);
    TString  t = mid ? "|#eta|<1" : "|#eta|>1";
    if (c1 >= c2) t.Append(" MB");
    else          t.Append(Form(" %4.1f#minus%4.1f%%", c1, c2));
    TLegend* l = (logy
		  ? p->BuildLegend(0.45,(mode==2?0.3:0.5),0.95,0.9, t) 
		  : p->BuildLegend(0.16,0.4,0.6,0.8, t));
    l->SetFillStyle(0);
    l->SetBorderSize(0);

    TFile* output = TFile::Open(Form("%s.root", c->GetName()),"RECREATE");
    output->cd();
    THStack* stacks[] = { Ps, Ss, Cs, Cb, 0 };
    THStack** ptr = stacks;
    while (*ptr) {
      THStack*    stack = *ptr;
      TString     name  = stack->GetTitle(); name.ReplaceAll(" ","");
      TDirectory* dir   = output->mkdir(name);
      dir->cd();
      stack->Write("all");
      TIter next(stack->GetHists());
      TH1*  hist = 0;
      while ((hist = static_cast<TH1*>(next()))) {
	name = hist->GetTitle(); name.ReplaceAll("-", "");
	hist->Write(name);
      }
      ptr++;
      output->cd();
    }
    output->Write();
       
    
    c->SaveAs(Form("%s.png",c->GetName()));
  }
};

/** 
 * Do back of the envoloope calculation  
 * 
 * @param file 
 * @param c1 
 * @param c2 
 * @param mid 
 *
 * @ingroup pwglf_forward_tracklets
 */  
void DeltaExpectations(Double_t c1=0, Double_t c2=5)
{
  Bool_t mid=true;
  DeltaCalculations cm;
  cm.Run(c1, c2, mid);
  DeltaCalculations cf;
  cf.Run(c1, c2, !mid);
}
//
// EOF
// 
