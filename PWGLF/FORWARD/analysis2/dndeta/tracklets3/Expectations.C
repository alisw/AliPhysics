/**
 * @file   BackOfTheEnvelope.C
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
const char* parFmt = "%-10s | %5.3f |";
const char* rowFmt = " %5.2f | %5.2f | %5.2f | %5.2f |";
/**
 * Class to do the calculation 
 * 
 * @ingroup pwglf_forward_tracklets
 */
struct Calculation
{
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
    Double_t Get(Int_t i) const
    {
      return (i == 0 ? Ps :
	      i == 1 ? Ss :
	      i == 2 ? Cs :
	      i == 3 ? Cb : 0);
    }
    /** 
     * Get a quantity
     * 
     * @param i number of parameter 
     * 
     * @return Parameter
     */
    Double_t GetE(Int_t i) const
    {
      return (i == 0 ? ePs :
	      i == 1 ? eSs :
	      i == 2 ? eCs :
	      i == 3 ? eCb : 0);
    }
    void Scale(Int_t i, Double_t s)
    {
      switch (i) {
      case 0: Ps *= s; ePs *= s; break;
      case 1: Ss *= s; eSs *= s; break;
      case 2: Cs *= s; eCs *= s; break;
      case 3: Cb *= s; eCb *= s; break;
      }
    }
    void Print() const
    {
      printf(rowFmt, Ps, Ss, Cs, Cb);
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
    Row      truth;
    Row& Get(Int_t i)
    {
      return (i == 0 ? reduced :
	      i == 1 ? expected :
	      i == 2 ? reweighed :
	      truth);
    }
    void Print()
    {
      printf(parFmt, n.Data(), w);
      reduced.  Print();
      expected. Print();
      reweighed.Print();
      truth.    Print();
      printf("\n");
    }
    Particle* Read(Int_t          iRow,
		   Int_t          column,
		   TCollection*   bin,
		   const TString& sub)
    {
      TString  pn  = "primaries";
      TString  sn  = "secondaries";
      TString  cn  = "combinatorics";
      THStack* hp  = GetHS(GetC(GetC(GetC(bin,pn),"specie"),sub),"rows");
      THStack* hs  = GetHS(GetC(GetC(GetC(bin,sn),"specie"),sub),"rows");
      THStack* hc  = GetHS(GetC(GetC(GetC(bin,cn),"specie"),sub),"rows");
      TH1*     rps = GetH1(hp->GetHists(), "rowSig");
      TH1*     rss = GetH1(hs->GetHists(), "rowSig");
      TH1*     rcs = GetH1(hc->GetHists(), "rowSig");
      TH1*     rcb = GetH1(hc->GetHists(), "rowBg");
      Row&     ret = (iRow == 0 ? reduced :
		      iRow == 1 ? expected :
		      iRow == 2 ? reweighed : truth);
      Double_t ww  = (iRow == 1 ? w : 1);
      if (n.IsNull()) n = rps->GetXaxis()->GetBinLabel(column);
      ret.Ps     = ww*rps->GetBinContent(column);
      ret.Ss     = ww*rss->GetBinContent(column);
      ret.Cs     = ww*rcs->GetBinContent(column);
      ret.Cb     = ww*rcb->GetBinContent(column);
      ret.ePs    = ww*rps->GetBinError  (column);
      ret.eSs    = ww*rss->GetBinError  (column);
      ret.eCs    = ww*rcs->GetBinError  (column);
      ret.eCb    = ww*rcb->GetBinError  (column);
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
  void Print()
  {
    printf("%-18s |", "Sample");
    printf(" %-29s |", "Reduced");
    printf(" %-29s |", "Expected");
    printf(" %-29s |", "Reweighed");
    printf(" %-29s |", "Truth");
    printf("\n");
    TString head(parFmt);
    for (Int_t i = 0; i < 4; i++) head.Append(rowFmt);
    head.ReplaceAll("%5.3f", "%5s");
    head.ReplaceAll("%5.2f", "%5s");
    TString lne(head);
    lne.ReplaceAll("|", "+");
    lne.ReplaceAll("%5s", "-----");
    lne.ReplaceAll("%10s", "----------");
    lne.ReplaceAll(" ", "-");
    Printf(head.Data(), "Type", "W",
	   "Ps", "Ss", "Cs", "Cb",
	   "Ps", "Ss", "Cs", "Cb",
	   "Ps", "Ss", "Cs", "Cb",
	   "Ps", "Ss", "Cs", "Cb");
    Printf(lne);
    for (ParticleList::const_iterator i = particles.begin();
	 i != particles.end(); ++i)
      (*i)->Print();
    Printf(lne);
  }
  void ScaleExpected()
  {
    for (Int_t i = 0; i < 4; i++) {
      Double_t sum = 0;
      for (ParticleList::const_iterator j = particles.begin();
	   j != particles.end(); ++j)
	sum += (*j)->expected.Get(i);
      for (ParticleList::const_iterator j = particles.begin();
	   j != particles.end(); ++j)
	(*j)->expected.Scale(i,100/sum);
    }
  }
      
      
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
  /** 
   * Create a particle 
   * 
   * @param w           Weight of particle 
   * @param column      Column in histograms to read off 
   * @param reduced     The reduced set 
   * @param reweighed   The reweighed set
   * @param truth       The truth set 
   * @param bin         The bin name 
   * @param sub         The sub name 
   * 
   * @return Create Particle object
   */
  Particle* Create(Double_t       w,
		   Int_t          column,		   
		   TCollection*   reduced,
		   TCollection*   reweighed,
		   TCollection*   truth,
		   const TString& sub)
  {
    Particle* ret = new Particle;
    ret->w = w;
    ret->Read(0, column, reduced,   sub);
    ret->Read(1, column, reduced,   sub);
    ret->Read(2, column, reweighed, sub);
    ret->Read(3, column, truth,     sub);
    particles.push_back(ret);
    return ret;
  }
  THStack* MakeStack(Int_t column, const char* name, const char* title)
  {
    const char*   names[]  = { "red", "exp", "rew", "trh" };
    const char*   titles[] = { "Reduced", "Expected", "Reweighed", "Truth" };
    const Color_t colors[] = { kRed+1, kMagenta+1, kBlue+1, kGreen+1 };
    Int_t         nPart    = particles.size();
    THStack* stack = new THStack(name, title);
    for (Int_t i = 0; i < 4; i++) {
      TH1* h = new TH1D(names[i],titles[i],nPart,0,nPart);
      h->SetFillColor(colors[i]);
      h->SetLineColor(colors[i]);
      h->SetMarkerColor(colors[i]);
      h->SetBarOffset(i*0.225+.075);
      h->SetBarWidth (0.20);
      h->SetDirectory(0);
      stack->Add(h, "bar0");
      
      for (Int_t j = 0; j < nPart; j++) {
	Particle* p = particles[j];
	Row&      r = p->Get(i);
	h->GetXaxis()->SetBinLabel(j+1,p->n);
	h->SetBinContent          (j+1,r.Get (column));
	h->SetBinError            (j+1,r.GetE(column));
      }
    }
    return stack;
  }
  void DrawStack(TVirtualPad* mother, Int_t sub, THStack* stack, Bool_t logy)
  {
    TVirtualPad* p = mother->cd(sub);
    p->SetGridx();
    p->SetGridy();
    p->SetTicks();
    if (logy) p->SetLogy();
    if ((p->GetNumber() % 2) == 0) p->SetRightMargin(0.01);
    stack->SetMinimum(0.02);
    stack->SetMaximum(logy ? 150 : 110);
    stack->Draw("nostack hist");
    stack->GetHistogram()->GetYaxis()->SetNdivisions(210);
    stack->GetHistogram()->GetYaxis()->SetLabelSize(0.03/p->GetHNDC());
    stack->GetHistogram()->GetYaxis()->SetTitleSize(0.03/p->GetHNDC());
    stack->GetHistogram()->GetXaxis()->SetLabelSize(0.03/p->GetHNDC());
    stack->GetHistogram()->GetYaxis()->SetTitle("Fraction [%]");
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
    TCollection* redTop = GetC(GetC(redFile,"MidRapidityMCResults"),bin.Data());
    TCollection* rweTop = GetC(GetC(rweFile,"MidRapidityMCResults"),bin.Data());
    TCollection* trhTop = GetC(GetC(trhFile,"MidRapidityMCResults"),bin.Data());

    Create(1.52233, 1, redTop, rweTop, trhTop, sub);
    Create(1.41178, 2, redTop, rweTop, trhTop, sub);
    Create(2.75002, 3, redTop, rweTop, trhTop, sub);
    Create(3.24110, 4, redTop, rweTop, trhTop, sub);
    Create(2.75002, 5, redTop, rweTop, trhTop, sub);
    Create(1,       6, redTop, rweTop, trhTop, sub);

    ScaleExpected();
    Print();

    THStack* Ps = MakeStack(0, "Ps", "Signal primaries");
    THStack* Ss = MakeStack(1, "Ss", "Signal secondaries");
    THStack* Cs = MakeStack(2, "Cs", "Signal fake");
    THStack* Cb = MakeStack(3, "Cb", "Background fake");

    Bool_t logy = true;
    TCanvas* c = new TCanvas("c","c",1000,1000);
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.01);
    c->SetLeftMargin(0.13);
    c->Divide(2,2,0,0);

    DrawStack(c, 1, Ps, logy);
    DrawStack(c, 2, Ss, logy);
    DrawStack(c, 3, Cs, logy);
    DrawStack(c, 4, Cb, logy);

    TVirtualPad* p = c->cd(1);
    TString  t = mid ? "|#eta|<1" : "|#eta|>1";
    if (c1 >= c2) t.Append(" MB");
    else          t.Append(Form(" %4.1f#minus%4.1f%%", c1, c2));
    TLegend* l = (logy
		  ? p->BuildLegend(0.45,0.5,0.95,0.9, t) 
		  : p->BuildLegend(0.16,0.4,0.6,0.8, t));
    l->SetFillStyle(0);
    l->SetBorderSize(0);

    c->SaveAs(Form("expected_%s_%s.png", bin.Data(), mid ? "mid" : "fwd"));
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
void Expectations(Bool_t mid=true, Double_t c1=0, Double_t c2=5)
{
  Calculation c;
  c.Run(c1, c2, mid);
}
//
// EOF
// 
