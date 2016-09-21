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
#endif

/**
 * Print format 
 * 
 * @ingroup pwglf_forward_tracklets
 */
const char* fmt = "%-10s | %5.3f | %5.2f | %5.2f | %5.2f  %5.2f |";
/**
 * Class to do the calculation 
 * 
 * @ingroup pwglf_forward_tracklets
 */
struct Calculation
{
  /** 
   * Per particle information
   * 
   * @ingroup pwglf_forward_tracklets
   */
  struct Particle
  {
    TString  n;
    Double_t w;
    Double_t Ps;
    Double_t Ss;
    Double_t Cs;
    Double_t Cb;
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
    void Print() const
    {
      Printf(fmt, n.Data(), w, Ps, Ss, Cs, Cb);
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
    TString head(fmt);
    head.ReplaceAll("%5.3f", "%5s");
    head.ReplaceAll("%5.2f", "%5s");
    TString lne(head);
    lne.ReplaceAll("|", "+");
    lne.ReplaceAll("%5s", "-----");
    lne.ReplaceAll("%10s", "----------");
    lne.ReplaceAll(" ", "-");
    Printf(head.Data(), "Type", "W", "Ps", "Ss", "Cs", "Cb");
    Printf(lne);
    Printf(fmt, "Full",   0., fP,  fS,  fC,  0.);
    Printf(fmt, "Signal", 0., fPs, fSs, fCs, 0.);
    Printf(lne);
    for (ParticleList::const_iterator i = particles.begin();
	 i != particles.end(); ++i)
      (*i)->Print();
    Printf(lne);
  }
  /** 
   * Add a particle 
   * 
   * @param p 
   */
  void AddParticle(Particle* p)
  {
    particles.push_back(p);
  }
  /** 
   * Calculate the change 
   * 
   * @param i 
   * @param w 
   * @param verb 
   * 
   * @return 
   */
  Double_t R(Int_t i, const char* w="", Bool_t verb=true) const
  {
    Double_t ret = 1;
    if (verb) printf("%-10s: %6.4f", w, ret);
    for (ParticleList::const_iterator j = particles.begin();
	 j != particles.end(); ++j) {
      Particle* p =  *j;
      ret         += (p->w-1)*(p->Get(i)/100);
      if (verb) printf("+(%4.2f-1)*%6.4f", (p->w), (p->Get(i)/100));
    }
    if (verb) printf("=%6.4f\n", ret);
    return ret;
  }  
  /** 
   * Calculate the full change 
   * 
   * @param verb 
   * 
   * @return 
   */
  Double_t Calc(Bool_t verb=true)
  {
    Double_t rPs      = R(0, "rPs", verb);
    Double_t rSs      = R(1, "rSs", verb);
    Double_t rCs      = R(2, "rCs", verb);
    Double_t rCb      = R(3, "rCb", verb);
    Double_t rMsim    = (fP*rPs+fS*rSs+fC*rCs)/(fP+fS+fC);
    Printf("%-10s: (%5.3f*%6.4f+%5.3f*%6.4f+%5.3f*%6.4f)/"
	   "(%5.3f+%5.3f+%5.3f)=%6.4f",
	   "rM'", fP, rPs, fS, rSs, fC, rCs, fP, fS, fC, rMsim);
    Double_t beta     = (fCs)/(fPs+fSs+fCs);
    Printf("%-10s: %5.3f/(%5.3f+%5.3f+%5.3f)=%6.4f", "beta",fCs,fPs,fSs,fCs);
    Double_t rBetasim = (1-rCs/rMsim*beta)/(1-beta);
    Printf("%-10s: (1-%6.4f/%6.4f*%6.4f)/(1-%6.4f)=%6.4f",
	   "rbeta'", rCs, rMsim, beta, beta);
    Double_t rScale   = rCs/rCb;
    Printf("%-10s: %6.4f/%6.4f=%6.4f", "rk", rCs, rCb, rScale);
    Double_t rBetarea = (1-rScale*beta)/(1-beta);
    Printf("%-10s: (1-%6.4f*%6.4f)/(1-%6.4f)=%6.4f",
	   "rbeta", rScale, beta, beta);   
    Double_t rR       = rPs * rBetarea / rBetasim / rMsim;
    Printf("%-10s: %6.4f*%6.4f/%6.4f/%6.4f=%6.4f",
	   "rR", rPs, rBetarea, rBetasim, rMsim);
    return rR;
  }
  TObject* ChkC(TObject* o, TClass* c) const
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
  TObject* ChkO(TObject* src, TObject* o, const char* name) const
  {
    if (!o) {
      ::Warning("ChkO", "Object %s not found in %s", name, src->GetName());
      src->ls();
      return 0;
    }
    return o;
  }
  TObject* GetO(TDirectory* d, const char* name, TClass* cls) const
  {
    if (!d) {
      ::Warning("GetO", "No directory passed for %s", name);
      return 0;
    }
    TObject* o = d->Get(name);
    return ChkC(ChkO(d, o, name), cls);    
  }
  TObject* GetO(TCollection* d, const char* name, TClass* cls) const
  {
    if (!d) {
      ::Warning("GetO", "No collection passed for %s", name);
      return 0;
    }
    TObject* o = d->FindObject(name);
    return ChkC(ChkO(d, o, name), cls);    
  }
  TCollection* GetC(TDirectory* d, const char* name)
  {
    return static_cast<TCollection*>(GetO(d, name, TCollection::Class()));
  }
  TCollection* GetC(TCollection* d, const char* name)
  {
    return static_cast<TCollection*>(GetO(d, name, TCollection::Class()));
  }
  TH1* GetH1(TCollection* d, const char* name)
  {
    return static_cast<TH1*>(GetO(d, name, TH1::Class()));
  }
  TH2* GetH2(TCollection* d, const char* name)
  {
    return static_cast<TH2*>(GetO(d, name, TH2::Class()));
  }
  THStack* GetHS(TCollection* d, const char* name)
  {
    return static_cast<THStack*>(GetO(d, name, THStack::Class()));
  }
  /** 
   * Create a particle 
   * 
   * @param w    Weight 
   * @param bin  Bin to look in 
   * @param hp   Stack for primaries 
   * @param hs   Stack for secondaries 
   * @param hc   Stack for combinatorics
   * 
   * @return The particle 
   */
  Particle* Create(Double_t w, Int_t bin,
		   THStack* hp, THStack* hs, THStack* hc)
  {
    TH1* rps = GetH1(hp->GetHists(), "rowSig");
    TH1* rss = GetH1(hs->GetHists(), "rowSig");
    TH1* rcs = GetH1(hc->GetHists(), "rowSig");
    TH1* rcb = GetH1(hc->GetHists(), "rowBg");
    Particle* ret = new Particle;
    ret->n   = rps->GetXaxis()->GetBinLabel(bin);
    ret->w   = w;
    ret->Ps  = rps->GetBinContent(bin);
    ret->Ss  = rss->GetBinContent(bin);
    ret->Cs  = rcs->GetBinContent(bin);
    ret->Cb  = rcb->GetBinContent(bin);
    AddParticle(ret);
    return ret;
  }  
  /** 
   * Run it 
   * 
   * @param fileName Input file 
   * @param binName  Bin name 
   * @param mid      If true, for |eta|<1, otherwise |eta|>1
   */
  void Run(const char* fileName, const char* binName, Bool_t mid=true)
  {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file) return;
    TString      sub  = (mid ? "mid" : "fwd");
    TString      oth  = (mid ? "fwd" : "mid");
    TCollection* top  = GetC(file, "MidRapidityMCResults");
    TCollection* bin  = GetC(top,  binName);
    TCollection* prim = GetC(GetC(GetC(bin,"primaries"),    "specie"),sub);
    TCollection* seco = GetC(GetC(GetC(bin,"secondaries"),  "specie"),sub);
    TCollection* comb = GetC(GetC(GetC(bin,"combinatorics"),"specie"),sub);
    TCollection* opri = GetC(GetC(GetC(bin,"primaries"),    "specie"),oth);
    TCollection* osec = GetC(GetC(GetC(bin,"secondaries"),  "specie"),oth);
    TCollection* ocom = GetC(GetC(GetC(bin,"combinatorics"),"specie"),oth);
    TH1*         ip   = GetH1(prim, "totalIntegrals");
    TH1*         is   = GetH1(seco, "totalIntegrals");
    TH1*         ic   = GetH1(comb, "totalIntegrals");
    TH1*         op   = GetH1(opri, "totalIntegrals");
    TH1*         os   = GetH1(osec, "totalIntegrals");
    TH1*         oc   = GetH1(ocom, "totalIntegrals");
    THStack*     hp   = GetHS(prim, "rows");
    THStack*     hs   = GetHS(seco, "rows");
    THStack*     hc   = GetHS(comb, "rows");
    Double_t     tot  = (ip->GetBinContent(1)+
			 is->GetBinContent(1)+
			 ic->GetBinContent(1)+
			 op->GetBinContent(1)+
			 os->GetBinContent(1)+
			 oc->GetBinContent(1));
    ip->Scale(1./tot);
    is->Scale(1./tot);
    ic->Scale(1./tot);
    fP                = ip->GetBinContent(1)*100;
    fS                = is->GetBinContent(1)*100;
    fC                = ic->GetBinContent(1)*100;
    fPs               = ip->GetBinContent(2)*100;
    fSs               = is->GetBinContent(2)*100;
    fCs               = ic->GetBinContent(2)*100;

    Create(1.52233, 1, hp, hs, hc);
    Create(1.41178, 2, hp, hs, hc);
    Create(2.75002, 3, hp, hs, hc);
    Create(3.24110, 4, hp, hs, hc);
    Create(2.75002, 5, hp, hs, hc);
    Create(1,       6, hp, hs, hc);

    TCanvas* c = new TCanvas("c","c");
    c->Divide(2,3);
    c->cd(1); ip->Draw("hist text30");
    c->cd(3); is->Draw("hist text30");
    c->cd(5); ic->Draw("hist text30");
    c->cd(2); hp->Draw("nostack");
    c->cd(4); hs->Draw("nostack");
    c->cd(6); hc->Draw("nostack");

    Print();
    
    Double_t r = Calc();
    Printf("R = %5.3f", r);

    // file->Close();
  }
  void Run(const char* fileName, Double_t c1=0, Double_t c2=0, Bool_t mid=true)
  {
    TString bin("all");
    if (c2 > c1)
      bin.Format("cent%03dd%02d_%03dd%02d",
		 Int_t(c1), Int_t(100*c1) % 100, 
		 Int_t(c2), Int_t(100*c2) % 100);
    Run(fileName, bin, mid);
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
void BackOfTheEnvelope(const char* file, Double_t c1=0, Double_t c2=0,
		       Bool_t mid=true)
{
  Calculation c;
  c.Run(file, c1, c2, mid);
}
//
// EOF
// 
