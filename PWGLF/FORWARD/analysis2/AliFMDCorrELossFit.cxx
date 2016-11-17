// Object holding the Energy loss fit 'correction' 
// 
// These are generated from Monte-Carlo or real ESDs. 
//
#include "AliFMDCorrELossFit.h"
#include "AliForwardUtil.h"
#include "AliLandauGaus.h"
#include <TF1.h>
#include <TGraph.h>
#include <TBrowser.h>
#include <TVirtualPad.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TH1D.h>
#include <AliLog.h>
#include <TMath.h>
#include <TList.h>
#include <iostream>
#include <iomanip>
namespace { 
  Int_t fDebug = 1;
}

Double_t AliFMDCorrELossFit::ELossFit::fgMaxRelError = .25;
Double_t AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-7;
Double_t AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 20;

//____________________________________________________________________
AliFMDCorrELossFit::ELossFit::ELossFit()
  : fN(0),
    fNu(0),
    fChi2(0),
    fC(0),
    fDelta(0),
    fXi(0),
    fSigma(0),
    fSigmaN(0),
    fA(0),
    fEC(0),
    fEDelta(0),
    fEXi(0),
    fESigma(0),
    fESigmaN(0),
    fEA(0),
    fQuality(0), 
    fDet(0), 
    fRing('\0'),
    fBin(0),
    fMaxWeight(0)
{
  //
  // Default constructor 
  // 
  //
}
//____________________________________________________________________
AliFMDCorrELossFit::ELossFit::ELossFit(Int_t quality, const TF1& f)
  : fN(f.GetNpar() > AliLandauGaus::kN ? 
       Int_t(f.GetParameter(AliLandauGaus::kN)) : 
       1),
    fNu(f.GetNDF()),
    fChi2(f.GetChisquare()),
    fC(f.GetParameter(AliLandauGaus::kC)),
    fDelta(f.GetParameter(AliLandauGaus::kDelta)),
    fXi(f.GetParameter(AliLandauGaus::kXi)),
    fSigma(f.GetParameter(AliLandauGaus::kSigma)),
    fSigmaN(f.GetParameter(AliLandauGaus::kSigmaN)),
    fA(0),
    fEC(f.GetParError(AliLandauGaus::kC)),
    fEDelta(f.GetParError(AliLandauGaus::kDelta)),
    fEXi(f.GetParError(AliLandauGaus::kXi)),
    fESigma(f.GetParError(AliLandauGaus::kSigma)),
    fESigmaN(f.GetParError(AliLandauGaus::kSigmaN)),
    fEA(0),
    fQuality(quality),
    fDet(0), 
    fRing('\0'),
    fBin(0),
    fMaxWeight(0)
{
  // 
  // Construct from a function
  // 
  // Parameters:
  //    quality Quality flag
  //    f       Function
  //
  if (fN <= 0) return;
  fA  = new Double_t[fN];
  fEA = new Double_t[fN];
  for (Int_t i = 0; i < fN-1; i++) { 
    fA[i]  = f.GetParameter(AliLandauGaus::kA+i);
    fEA[i] = f.GetParError(AliLandauGaus::kA+i);
  }
  fA[fN-1]  = -9999;
  fEA[fN-1] = -9999;
}

//____________________________________________________________________
AliFMDCorrELossFit::ELossFit::ELossFit(Int_t     quality,UShort_t  n, 
				       Double_t  chi2,   UShort_t  nu, 
				       Double_t  c,      Double_t  ec, 
				       Double_t  delta,  Double_t  edelta, 
				       Double_t  xi,     Double_t  exi,
				       Double_t  sigma,  Double_t  esigma, 
				       Double_t  sigman, Double_t  esigman, 
				       const Double_t* a,const Double_t* ea)
  : fN(n),
    fNu(nu),
    fChi2(chi2),
    fC(c),
    fDelta(delta),
    fXi(xi),
    fSigma(sigma),
    fSigmaN(sigman),
    fA(0),
    fEC(ec),
    fEDelta(edelta),
    fEXi(exi),
    fESigma(esigma),
    fESigmaN(esigman),
    fEA(0),
    fQuality(quality),
    fDet(0), 
    fRing('\0'),
    fBin(0),
    fMaxWeight(0)
{
  // 
  // Constructor with full parameter set
  // 
  // Parameters:
  //    quality   Quality flag
  //    n         @f$ N@f$ - Number of fitted peaks
  //    chi2      @f$ \chi^2 @f$
  //    nu        @f$ \nu @f$ - number degrees of freedom
  //    c         @f$ C@f$ - scale constant
  //    ec        @f$ \delta C@f$ - error on @f$ C@f$ 
  //    delta     @f$ \Delta@f$ - Most probable value		  
  //    edelta    @f$ \delta\Delta@f$ - error on @f$\Delta@f$ 
  //    xi        @f$ \xi@f$ - width  
  //    exi       @f$ \delta\xi@f$ - error on @f$\xi@f$ 
  //    sigma     @f$ \sigma@f$ - Width of Gaussian		   
  //    esigma    @f$ \delta\sigma@f$ - error on @f$\sigma@f$ 
  //    sigman    @f$ \sigma_n@f$ - Noise width		  
  //    esigman   @f$ \delta\sigma_n@f$ - error on @f$\sigma_n@f$ 
  //    a         Array of @f$ N-1@f$ weights @f$ a_i@f$ for 
  //                  @f$ i=2,\ldots@f$ 
  //    ea        Array of @f$ N-1@f$ error on the weights @f$ a_i@f$ for 
  //                  @f$ i=2,\ldots@f$ 
  //
  if (fN <= 0) return;
  fA  = new Double_t[fN];
  fEA = new Double_t[fN];
  for (Int_t i = 0; i < fN-1; i++) { 
    fA[i]  = a[i];
    fEA[i] = ea[i];
  }
  fA[fN-1]  = -9999;
  fEA[fN-1] = -9999;
}
//____________________________________________________________________
AliFMDCorrELossFit::ELossFit::ELossFit(const ELossFit& o)
  : TObject(o), 
    fN(o.fN),
    fNu(o.fNu),
    fChi2(o.fChi2),
    fC(o.fC),
    fDelta(o.fDelta),
    fXi(o.fXi),
    fSigma(o.fSigma),
    fSigmaN(o.fSigmaN),
    fA(0),
    fEC(o.fEC),
    fEDelta(o.fEDelta),
    fEXi(o.fEXi),
    fESigma(o.fESigma),
    fESigmaN(o.fESigmaN),
    fEA(0),
    fQuality(o.fQuality),
    fDet(o.fDet), 
    fRing(o.fRing),
    fBin(o.fBin),
    fMaxWeight(o.fMaxWeight)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  if (fN <= 0) return;
  fA  = new Double_t[fN];
  fEA = new Double_t[fN];
  for (Int_t i = 0; i < fN-1; i++) { 
    fA[i]  = o.fA[i];
    fEA[i] = o.fEA[i];
  }
  fA[fN-1]  = -9999;
  fEA[fN-1] = -9999;
}

//____________________________________________________________________
AliFMDCorrELossFit::ELossFit&
AliFMDCorrELossFit::ELossFit::operator=(const ELossFit& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
  if (&o == this) return *this; 
  fN	     = o.fN;
  fNu	     = o.fNu;
  fChi2	     = o.fChi2;
  fC	     = o.fC;
  fDelta     = o.fDelta;
  fXi	     = o.fXi;
  fSigma     = o.fSigma;
  fSigmaN    = o.fSigmaN;
  fEC	     = o.fEC;
  fEDelta    = o.fEDelta;
  fEXi	     = o.fEXi;
  fESigma    = o.fESigma;
  fESigmaN   = o.fESigmaN;
  fQuality   = o.fQuality;
  fDet       = o.fDet; 
  fRing      = o.fRing;
  fBin       = o.fBin;
  fMaxWeight = o.fMaxWeight;
  if (fA)  delete [] fA;
  if (fEA) delete [] fEA; 
  fA  = 0;
  fEA = 0;

  if (fN <= 0) return *this;
  fA  = new Double_t[fN];
  fEA = new Double_t[fN];
  for (Int_t i = 0; i < fN; i++) { 
    fA[i]  = o.fA[i];
    fEA[i] = o.fEA[i];
  }

  return *this;
}

//____________________________________________________________________
AliFMDCorrELossFit::ELossFit::~ELossFit()
{
  if (fA)  delete[] fA;
  if (fEA) delete[] fEA;
}


//____________________________________________________________________
Int_t
AliFMDCorrELossFit::ELossFit::FindMaxWeight(Double_t maxRelError, 
					    Double_t leastWeight, 
					    UShort_t maxN) const
{
  // 
  // Find the maximum weight to use.  The maximum weight is the
  // largest i for which 
  // 
  // - @f$ i \leq \max{N}@f$ 
  // - @f$ a_i > \min{a}@f$ 
  // - @f$ \delta a_i/a_i > \delta_{max}@f$ 
  // 
  // Parameters:
  //    maxRelError @f$ \min{a}@f$ 
  //    leastWeight @f$ \delta_{max}@f$ 
  //    maxN        @f$ \max{N}@f$      
  // 
  // Return:
  //    The largest index @f$ i@f$ for which the above
  // conditions hold.  Will never return less than 1. 
  //
  if (fMaxWeight > 0) return fMaxWeight;
  // Info("FindMaxWeight", "Maximum weight (fN=%d, maxN=%d)", fN, maxN);
  Int_t n = TMath::Min(maxN, UShort_t(fN-1));
  Int_t m = 1;
  // fN is one larger than we have data 
  for (Int_t i = 0; i < n; i++, m++) {
    // Info("FindMaxWeight", "Investigating fA[%d/%d]: %f +/- %f (>%g,%f/%f)",
    //      i, n, fA[i], fEA[i], leastWeight, fEA[i] / fA[i], maxRelError);
    if (fA[i] < leastWeight)  break;
    if (fEA[i] / fA[i] > maxRelError) break;
  }
  // Info("FindMaxWeight","Found %d",m);
  return fMaxWeight = m;
}

//____________________________________________________________________
Double_t 
AliFMDCorrELossFit::ELossFit::Evaluate(Double_t x, 
				       UShort_t maxN) const
{
  // 
  // Evaluate 
  // @f[ 
  //  f_N(x;\Delta,\xi,\sigma') = 
  //     \sum_{i=1}^{n} a_i f(x;\Delta_i,\xi_i,\sigma_i')
  // @f] 
  //
  // (see AliLandauGaus::NLandauGaus) for the maximum @f$ N @f$
  // that fulfills the requirements 
  // 
  // Parameters:
  //    x           Where to evaluate 
  //    maxN 	  @f$ \max{N}@f$    
  // 
  // Return:
  //    @f$ f_N(x;\Delta,\xi,\sigma')@f$ 
  //
  return AliLandauGaus::Fn(x, fDelta, fXi, fSigma, fSigmaN, 
			   TMath::Min(maxN, UShort_t(fN)), fA);
}

//____________________________________________________________________
Double_t 
AliFMDCorrELossFit::ELossFit::EvaluateWeighted(Double_t x, 
					       UShort_t maxN) const
{									
  // 
  // Evaluate 
  // @f[ 
  //   f_W(x;\Delta,\xi,\sigma') = 
  //   \frac{\sum_{i=1}^{n} i a_i f_i(x;\Delta,\xi,\sigma')}{
  //     f_N(x;\Delta,\xi,\sigma')} = 
  //   \frac{\sum_{i=1}^{n} i a_i f(x;\Delta_i,\xi_i,\sigma_i')}{
  //     \sum_{i=1}^{n} a_i f(x;\Delta_i,\xi_i,\sigma_i')}
  // @f] 
  // where @f$ n@f$ fulfills the requirements (see FindMaxWeight). 
  //
  // If the denominator is zero, then 1 is returned. 
  //
  // See also AliLandauGaus::Fi and AliLandauGaus::NLandauGaus
  // for more information on the evaluated functions. 
  // 
  // Parameters:
  //    x           Where to evaluate 
  //    maxN 	  @f$ \max{N}@f$      
  // 
  // Return:
  //    @f$ f_W(x;\Delta,\xi,\sigma')@f$.  
  //
  UShort_t n   = TMath::Min(maxN, UShort_t(fN-1));
  Double_t num = 0;
  Double_t den = 0;
  for (Int_t i = 1; i <= n; i++) {
    Double_t a = (i == 1 ? 1 : fA[i-1]);
    if (fA[i-1] < 0) break;
    Double_t f = AliLandauGaus::Fi(x,fDelta,fXi,fSigma,fSigmaN,i);
    num += i * a * f;
    den += a * f;
  }
  if (den <= 0) return 1;
  return num / den;
}


#define OUTPAR(N,V,E) 			\
  std::setprecision(9) <<               \
  std::setw(12) << N << ": " << 	\
  std::setw(14) << V << " +/- " <<  	\
  std::setw(14) << E << " (" <<  	\
  std::setprecision(-1) <<              \
  std::setw(5)  << 100*(V>0?E/V:1) << "%)\n"
  

//____________________________________________________________________
Int_t
AliFMDCorrELossFit::ELossFit::Compare(const TObject* o) const
{
  // 
  // Compare to another ELossFit object. 
  // 
  // - +1, if this quality is better (larger) than other objects quality
  // - -1, if this quality is worse (smaller) than other objects quality
  // - +1, if this @f$|\chi^2/\nu-1|@f$ is smaller than the same for other
  // - -1, if this @f$|\chi^2/\nu-1|@f$ is larger than the same for other
  // - +1, if this @f$ N_{peak}@f$ is larger than the the same for other
  // - -1, if this @f$ N_{peak}@f$ is smaller than the the same for other
  // - 0 otherwise 
  // 
  // Parameters:
  //    o Other object to compare to 
  //
  const ELossFit* other = static_cast<const ELossFit*>(o);
  if (this->fQuality < other->fQuality) return -1;
  if (this->fQuality > other->fQuality) return +1;
  Double_t chi2nu  = (fNu == 0 ? 1000 : fChi2 / fNu);
  Double_t oChi2nu = (other->fNu == 0 ? 1000 : other->fChi2 / other->fNu);
  if (TMath::Abs(chi2nu-1) < TMath::Abs(oChi2nu-1)) return +1;
  if (TMath::Abs(chi2nu-1) > TMath::Abs(oChi2nu-1)) return -1;
  if (fN > other->fN)                               return +1;
  if (fN < other->fN)                               return -1;
  return 0;
}

//____________________________________________________________________
void
AliFMDCorrELossFit::ELossFit::Print(Option_t* option) const
{
  // 
  // Information to standard output 
  // 
  // Parameters:
  //    option Not used 
  //
  TString o(option);
  if (o.Contains("S", TString::kIgnoreCase)) {
    Printf("%15s: q=%2d n=%1d chi2/nu=%6.3f",
	   GetName(), fQuality, fN, (fNu <= 0 ? 999 : fChi2 / fNu));
    return;
  }
  
  std::cout << GetName() << ":\n"
	    << " chi^2/nu = " << fChi2 << "/" << fNu << " = " 
	    << (fNu == 0 ? 999 : fChi2 / fNu) << "\n"
	    << "     Quality:   " << fQuality << "\n" 
	    << "  NParticles:   " << fN << "  (" << FindMaxWeight() << ")\n"
	    << OUTPAR("Delta",   fDelta,  fEDelta) 
	    << OUTPAR("xi",      fXi,     fEXi)
	    << OUTPAR("sigma",   fSigma,  fESigma)
	    << OUTPAR("sigma_n", fSigmaN, fESigmaN);
  for (Int_t i = 0; i < fN-1; i++) 
    std::cout << OUTPAR(Form("a%d", i+2), fA[i], fEA[i]);
  std::cout << std::flush;
}
//____________________________________________________________________
TF1*
AliFMDCorrELossFit::ELossFit::GetF1(Int_t i, Double_t max) const
{
  const Double_t lowX = 0.001;
  const Double_t upX  = (max < 0 ? 10 : max);
  Int_t          maxW = FindMaxWeight();
  TF1*           ret  = 0;
  
  if (i <= 0)
    ret = AliLandauGaus::MakeFn(fC * 1, fDelta, fXi, 
				fSigma, fSigmaN, maxW/*fN*/, fA,  lowX, upX);
  else if (i == 1) 
    ret = AliLandauGaus::MakeF1(fC, fDelta, fXi, fSigma, fSigmaN, lowX, upX);
  else if (i <= maxW) 
    ret = AliLandauGaus::MakeFi(fC*(i == 1 ? 1 : fA[i-2]), 
				fDelta, fXi, fSigma, fSigmaN, i, lowX, upX);
  
  return ret;
}

//____________________________________________________________________
const Char_t*
AliFMDCorrELossFit::ELossFit::GetName() const 
{
  // 
  // Get the name of this object 
  // 
  // 
  // Return:
  //    
  //
  return Form("FMD%d%c_etabin%03d", fDet, fRing, fBin);
}

//____________________________________________________________________
void
AliFMDCorrELossFit::ELossFit::Browse(TBrowser* b)
{
  // 
  // Browse this object 
  // 
  // Parameters:
  //    b Browser
  //
  Draw(b ? b->GetDrawOption() : "comp values");
  gPad->SetLogy();
  gPad->Update();
}
     
//____________________________________________________________________
void
AliFMDCorrELossFit::ELossFit::Draw(Option_t* option)
{
  // 
  // Draw this fit 
  // 
  // Parameters:
  //    option Options 
  //  - COMP  Draw components too 
  //
  TString opt(option);
  opt.ToUpper();
  bool comp = false;
  bool good = false;
  bool vals = false;
  bool legd = false;
  bool peak = false;
  if (opt.Contains("COMP")) { 
    opt.ReplaceAll("COMP","");
    comp = true;
  }
  if (opt.Contains("GOOD")) { 
    opt.ReplaceAll("GOOD","");
    good = true;
  }
  if (opt.Contains("VALUES")) { 
    opt.ReplaceAll("VALUES","");
    vals = true;
  }
  if (opt.Contains("LEGEND")) { 
    opt.ReplaceAll("LEGEND","");
    legd = comp;
  }
  if (!opt.Contains("SAME")) { 
    gPad->Clear();
  }
  if (opt.Contains("PEAK")) { 
    peak = true;
  }
  TLegend* l = 0;
  if (legd) { 
    l = new TLegend(.3, .5, .59, .94);
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetFillStyle(0);
  }
  TObjArray cleanup;
  Int_t maxW = FindMaxWeight();
  TF1* tot = AliLandauGaus::MakeFn(fC * 1, fDelta, fXi, fSigma, fSigmaN, 
				   maxW/*fN*/,     fA,  0.01,   10);
  tot->SetLineColor(kBlack);
  tot->SetLineWidth(2);
  tot->SetLineStyle(1);
  tot->SetTitle(GetName());
  if (l) l->AddEntry(tot, "Total", "l");
  Double_t max = tot->GetMaximum();

  
  if (!opt.Contains("SAME")) {
    TH1* frame = new TH1F(GetName(), 
		     Form("FMD%d%c, eta bin %d",fDet,fRing,fBin),
		     100, 0, 10);
    frame->SetMinimum(max/10000);
    frame->SetMaximum(max*1.4);
    frame->SetXTitle("#Delta / #Delta_{mip}");
    frame->Draw();
    opt.Append(" SAME");
  }
  TF1* cpy = tot->DrawCopy(opt.Data());
  cleanup.Add(tot);

  if (vals) { 
    Double_t x1 = .72;
    Double_t x2 = .73;
    Double_t y  = .90;
    Double_t dy = .05;
    TLatex* ltx1 = new TLatex(x1, y, "");
    TLatex* ltx2 = new TLatex(x2, y, "");
    ltx1->SetNDC();
    ltx1->SetTextAlign(33);
    ltx1->SetTextFont(132);
    ltx1->SetTextSize(dy-.01);
    ltx2->SetNDC();
    ltx2->SetTextAlign(13);
    ltx2->SetTextFont(132);
    ltx2->SetTextSize(dy-.01);

    ltx1->DrawLatex(x1, y, "Quality");
    ltx2->DrawLatex(x2, y, Form("%d", fQuality));
    y -= dy;

    ltx1->DrawLatex(x1, y, "#chi^{2}/#nu");
    ltx2->DrawLatex(x2, y, Form("%7.3f", (fNu > 0 ? fChi2 / fNu : -1)));
    y -= dy;
    
    const Char_t* pn[] = { "C", "#Delta", "#xi", "#sigma" };
    Double_t      pv[] = { fC,  fDelta,  fXi,  fSigma };
    Double_t      pe[] = { fEC, fEDelta, fEXi, fESigma };
    for (Int_t i = 0; i < 4; i++) { 
      ltx1->DrawLatex(x1, y, pn[i]);
      ltx2->DrawLatex(x2, y, Form("%6.4f#pm%6.4f", pv[i], pe[i]));
      y -= dy;
    }
    for (Int_t i=2; i <= fN; i++) { 
      if (i > maxW) {
	ltx1->SetTextColor(kRed+3);
	ltx2->SetTextColor(kRed+3);
      }
      ltx1->DrawLatex(x1, y, Form("a_{%d}", i));
      ltx2->DrawLatex(x2, y, Form("%6.4f#pm%6.4f", fA[i-2], fEA[i-2]));
      y -= dy;
    }

  }

  if (peak) { 
    TLine* pl = new TLine(fDelta, 0.01*max, fDelta, 1.5*max);
    pl->SetLineStyle(2);
    pl->SetLineColor(kBlack);
    pl->Draw();
  }
  if (!comp) { 
    gPad->cd();
    return;
  }

  Double_t min = max;
  opt.Append(" same");
  for (Int_t i=1; i <= fN; i++) { 
    if (good && i > maxW) break;
    TF1* f = AliLandauGaus::MakeFi(fC*(i == 1 ? 1 : fA[i-2]), 
					     fDelta, fXi, 
					     fSigma, fSigmaN, 
					     i,      0.01, 10);
    f->SetLineWidth(2);
    f->SetLineStyle(i > maxW ? 2 : 1);
    min = TMath::Min(f->GetMaximum(), min);
    f->DrawCopy(opt.Data());
    if (l) l->AddEntry(f, Form("%d MIP%s", i, (i>1 ? "s" : "")), "l");

    cleanup.Add(f);
  }
  min /= 100;
  if (max <= 0) max = 0.1;
  if (min <= 0) min = 1e-4;
  cpy->SetMaximum(max);
  cpy->SetMinimum(min);
  cpy->GetHistogram()->SetMaximum(max);
  cpy->GetHistogram()->SetMinimum(min);
  cpy->GetHistogram()->GetYaxis()->SetRangeUser(min, max);
  if (l) l->Draw();

  gPad->cd();
}

						 
//____________________________________________________________________
#define CHECKPAR(V,E,T) ((V > 0) && (E / V < T))

//____________________________________________________________________
Double_t
AliFMDCorrELossFit::ELossFit::GetMpvCut(Double_t f) const
{
  return f * fDelta;
}
//____________________________________________________________________
Double_t
AliFMDCorrELossFit::ELossFit::GetXiCut(Double_t f) const
{
  return fDelta - f * fXi;
}
//____________________________________________________________________
Double_t
AliFMDCorrELossFit::ELossFit::GetXiSigmaCut(Double_t f) const
{
  return fDelta - f * (fXi+fSigma);
}

//____________________________________________________________________
Double_t
AliFMDCorrELossFit::ELossFit::GetAvgXiSigmaCut(Double_t f) const
{
  Double_t sum = fXi / (fEXi*fEXi) + fSigma / (fESigma+fESigma);
  Double_t nor = 1   / (fEXi*fEXi) + 1      / (fESigma+fESigma);
  return fDelta - f * sum / nor;
}

//____________________________________________________________________
Double_t
AliFMDCorrELossFit::ELossFit::FindProbabilityCut(Double_t low) const
{
  Double_t ret = 1000;
  TF1*     f   = 0;
  TGraph*  g   = 0;
  try {
    if (!(f = GetF1(1,5))) // First landau up to Delta/Delta_{mip}=5
      throw TString("Didn't TF1 object");
    if (!(g = new TGraph(f, "i")))
      throw TString("Failed to integrate function");

    Int_t    n     = g->GetN();
    Double_t total = g->GetY()[n-1];
    if (total <= 0) 
      throw TString::Format("Invalid integral: %lf", total);

    for (Int_t i = 0; i < n; i++) { 
      Double_t prob = g->GetY()[i] / total;
      if (prob > low) {
	ret = g->GetX()[i];
	break;
      }
    }
    if (ret >= 1000) 
      throw TString::Format("Couldn't find x value for cut %lf", low);
  }
  catch (const TString& str) {
    AliWarningF("%s: %s", GetName(), str.Data());
  }
  if (f) delete f;
  if (g) delete g;
  return ret;
}

//____________________________________________________________________
void 
AliFMDCorrELossFit::ELossFit::CalculateQuality(Double_t maxChi2nu, 
					       Double_t maxRelError, 
					       Double_t leastWeight)
{
  // 
  // Calculate the quality 
  //
  Double_t decline = maxChi2nu;
  Int_t qual = 0;
  if (fNu > 0) {
    Double_t red = fChi2 / fNu;
    if (red < maxChi2nu) qual += 4;
    else {
      Int_t q = Int_t((maxChi2nu+decline - red) / decline * 4);
      if (q > 0) qual += q;
    }
  }
  if (CHECKPAR(fDelta,  fEDelta,  maxRelError)) qual++;
  if (CHECKPAR(fXi,     fEXi,     maxRelError)) qual++;
  if (CHECKPAR(fSigma,  fESigma,  maxRelError)) qual++;
  if (CHECKPAR(fSigmaN, fESigmaN, maxRelError)) qual++;
  // Large penalty for large sigma - often a bad fit. (factor was 10)
  if (fSigma > 4*fXi)                          qual -= 4;
  if (fXi    < 0.01)                           qual -= 2;
  if (fSigma < 0.01)                           qual -= 2;
  //Info("CalculateQuality", "FMD%d%c [%3d] sigma=%f xi=%f sigma/xi=%f qual=%d",
  //     fDet, fRing, fBin, fSigma, fXi, fSigma/fXi, qual);
  qual += FindMaxWeight(2*maxRelError, leastWeight, fN);
  fQuality = qual;
}

//____________________________________________________________________
AliFMDCorrELossFit::AliFMDCorrELossFit()
  : TObject(), 
    fRings(), 
    fEtaAxis(0,0,0), 
    fLowCut(0),
    fCache(0)
{
  // 
  // Default constructor 
  //
  fRings.SetOwner(kTRUE);
  fEtaAxis.SetTitle("#eta");
  fEtaAxis.SetName("etaAxis");
  fRings.SetName("rings");
  SetBit(kIsGoodAsserted,false);
}

//____________________________________________________________________
AliFMDCorrELossFit::AliFMDCorrELossFit(const AliFMDCorrELossFit& o)
  : TObject(o), 
    fRings(o.fRings),
    fEtaAxis(o.fEtaAxis.GetNbins(),o.fEtaAxis.GetXmin(),o.fEtaAxis.GetXmax()), 
    fLowCut(0),
    fCache(0)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  fEtaAxis.SetTitle("#eta");
  fEtaAxis.SetName("etaAxis");
  SetBit(kIsGoodAsserted,false);
}

//____________________________________________________________________
AliFMDCorrELossFit::~AliFMDCorrELossFit()
{
  // 
  // Destructor 
  //
  fRings.Clear();
}

//____________________________________________________________________
AliFMDCorrELossFit&
AliFMDCorrELossFit::operator=(const AliFMDCorrELossFit& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
  if (&o == this) return *this; 
  fRings = o.fRings;
  fLowCut = o.fLowCut;
  fCache  = o.fCache;
  SetEtaAxis(o.fEtaAxis.GetNbins(), o.fEtaAxis.GetXmin(), o.fEtaAxis.GetXmax());

  return *this;
}
#define CACHE(BIN,IDX,OFFSET) fCache[IDX*OFFSET+BIN-1]
#define CACHEDR(BIN,D,R,OFFSET) \
  CACHE(BIN,(D == 1 ? 0 : (D - 2) * 2 + 1 + (R=='I' || R=='i' ? 0 : 1)),OFFSET)

//____________________________________________________________________
void
AliFMDCorrELossFit::CacheBins(UShort_t minQuality) const
{
  AliLandauGaus::EnableSigmaShift(TestBit(kHasShift));
  if (fCache.GetSize() > 0) return;

  Int_t nRings = fRings.GetEntriesFast();
  Int_t offset = fEtaAxis.GetNbins();

  fCache.Set(nRings*offset);
  fCache.Reset(-1);
  
  for (Int_t i = 0; i < nRings; i++) { 
    TObjArray* ringArray  = static_cast<TObjArray*>(fRings.At(i));

    // First loop to find where we actually have fits
    Int_t nFits      = 0;
    Int_t nGood      = 0;
    Int_t minBin     = offset+1;
    Int_t maxBin     = -1;
    Int_t realMinBin = offset+1;
    Int_t realMaxBin = -1;
    for (Int_t j = 1; j < ringArray->GetEntriesFast(); j++) {
      ELossFit* fit = static_cast<ELossFit*>(ringArray->At(j));
      if (!fit) continue;
      nFits++;

      // Update our range off possible fits 
      realMinBin = TMath::Min(j, realMinBin);
      realMaxBin = TMath::Max(j, realMaxBin);
      
      // Check the quality of the fit 
      fit->CalculateQuality(AliFMDCorrELossFit::ELossFit::fgMaxChi2nu, 
			    AliFMDCorrELossFit::ELossFit::fgMaxRelError, 
			    AliFMDCorrELossFit::ELossFit::fgLeastWeight);
      if (minQuality > 0 && fit->fQuality < minQuality) continue;
      nGood++;
      
      // Check this bin 
      CACHE(j,i,offset) = j;
      minBin            = TMath::Min(j, minBin);
      maxBin            = TMath::Max(j, maxBin);
    }
    AliInfoF("Out of %d bins, %d had fits, of which %d are good (%5.1f%%)", 
	     offset, nFits, nGood, (nFits > 0 ? 100*float(nGood)/nFits : 0));
    
    // Now loop and set neighbors 
    realMinBin = TMath::Max(1,      realMinBin-1); // Include one more 
    realMaxBin = TMath::Min(offset, realMaxBin+1); // Include one more 
    for (Int_t j = realMinBin; j <= realMaxBin; j++) {
      if (CACHE(j,i,offset) > 0) continue;
      
      Int_t nK    = TMath::Max(realMaxBin - j, j - realMinBin);
      Int_t found = -1;
      for (Int_t k = 1; k <= nK; k++) {
	Int_t left  = j - k;
	Int_t right = j + k;
	if      (left  > realMinBin && 
		 CACHE(left,i,offset)  == left) found = left;
	else if (right < realMaxBin && 
		 CACHE(right,i,offset) == right) found = right;
	if (found > 0) break;
      }
      // Now check that we found something 
      if (found) CACHE(j,i,offset) = CACHE(found,i,offset);
      else AliWarningF("No fit found for etabin=%d in ring=%d", j, i);
    }
  }
}


//____________________________________________________________________
Int_t
AliFMDCorrELossFit::FindEtaBin(Double_t eta) const
{
  // 
  // Find the eta bin corresponding to the given eta 
  // 
  // Parameters:
  //    eta  Eta value 
  // 
  // Return:
  //    Bin (in @f$[1,N_{bins}]@f$) corresponding to the given
  // eta, or 0 if out of range.
  //
  if (TMath::Abs(fEtaAxis.GetXmin() - fEtaAxis.GetXmax()) < 1e-6 
      || fEtaAxis.GetNbins() == 0) {
    AliWarning("No eta axis defined");
    return -1;
  }
  Int_t bin = const_cast<TAxis&>(fEtaAxis).FindBin(eta);
  if (bin <= 0 || bin > fEtaAxis.GetNbins()) return 0;
  return bin;
}

//____________________________________________________________________
Bool_t
AliFMDCorrELossFit::SetFit(UShort_t d, Char_t r, Int_t etaBin, ELossFit* fit)
{
  // 
  // Set the fit parameters from a function 
  // 
  // Parameters:
  //    d       Detector
  //    r       Ring 
  //    etaBin  Eta (bin number, 1->nBins)
  //    f       ELoss fit result - note, the object will take ownership
  //  
  TObjArray* ringArray = GetOrMakeRingArray(d, r);
  if (!ringArray) { 
    AliError(Form("Failed to make ring array for FMD%d%c", d, r));
    return kFALSE;
  }
  if (etaBin <= 0 || etaBin >= fEtaAxis.GetNbins()+1) { 
    AliError(Form("bin=%d is out of range [%d,%d]", 
		  etaBin, 1, fEtaAxis.GetNbins()));
    return kFALSE;
  }
  // AliInfo(Form("Adding fit %p at %3d", fit, etaBin));
  ringArray->AddAtAndExpand(fit, etaBin);
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDCorrELossFit::SetFit(UShort_t d, Char_t r, Double_t eta, ELossFit* fit)
{
  // 
  // Set the fit parameters from a function 
  // 
  // Parameters:
  //    d    Detector
  //    r    Ring 
  //    eta  Eta 
  //    f    ELoss fit result - note, the object will take ownership
  //  
  Int_t bin = FindEtaBin(eta);
  if (bin <= 0) { 
    AliError(Form("eta=%f is out of range [%f,%f]", 
		  eta, fEtaAxis.GetXmin(), fEtaAxis.GetXmax()));
    return kFALSE;
  }

  return SetFit(d, r, bin, fit);
}
//____________________________________________________________________
Bool_t
AliFMDCorrELossFit::SetFit(UShort_t  d,      Char_t    r, 
			   Double_t  eta, 
			   Int_t     quality,UShort_t  n, 
			   Double_t  chi2,   UShort_t  nu, 
			   Double_t  c,      Double_t  ec, 
			   Double_t  delta,  Double_t  edelta, 
			   Double_t  xi,     Double_t  exi,
			   Double_t  sigma,  Double_t  esigma, 
			   Double_t  sigman, Double_t  esigman, 
			   Double_t* a,      Double_t* ea)
{
  // 
  // Set the fit parameters from a function 
  // 
  // Parameters:
  //    d         Detector number
  //    r         Ring identifier 
  //    eta       Eta value
  //    quality   Quality flag
  //    n         @f$ N@f$ - Number of fitted peaks
  //    chi2      @f$ \chi^2 @f$
  //    nu        @f$ \nu @f$ - number degrees of freedom
  //    c         @f$ C@f$ - scale constant
  //    ec        @f$ \delta C@f$ - error on @f$ C@f$ 
  //    delta     @f$ \Delta@f$ - most probable value
  //    edelta    @f$ \delta\Delta@f$ - error on @f$\Delta@f$ 
  //    xi        @f$ \xi@f$ - Landau width		  
  //    exi       @f$ \delta\xi@f$ - error on @f$\xi@f$ 
  //    sigma     @f$ \sigma@f$ - Gaussian width
  //    esigma    @f$ \delta\sigma@f$ - error on @f$\sigma@f$ 
  //    sigman    @f$ \sigma_n@f$ - Noise width		  
  //    esigman   @f$ \delta\sigma_n@f$ - error on @f$\sigma_n@f$ 
  //    a         Array of @f$ N-1@f$ weights @f$ a_i@f$ for 
  //                  @f$ i=2,\ldots@f$ 
  //    ea        Array of @f$ N-1@f$ errors on weights @f$ a_i@f$ for 
  //                  @f$ i=2,\ldots@f$ 
  //
  ELossFit* e = new ELossFit(quality, n, 
			     chi2,    nu,
			     c,       ec,
			     delta,   edelta,
			     xi,      exi,
			     sigma,   esigma,
			     sigman,  esigman,
			     a,       ea);
  if (!SetFit(d, r, eta, e)) { 
    delete e;
    return kFALSE;
  }
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDCorrELossFit::SetFit(UShort_t  d, Char_t r, Double_t eta, 
			   Int_t quality, const TF1& f)
{
  // 
  // Set the fit parameters from a function 
  // 
  // Parameters:
  //    d        Detector
  //    r        Ring 
  //    eta      Eta 
  //    quality  Quality flag
  //    f        Function from fit 
  //  
  ELossFit* e = new ELossFit(quality, f);
  if (!SetFit(d, r, eta, e)) { 
    delete e;
    return kFALSE;
  }
  return kTRUE;
}
//____________________________________________________________________
AliFMDCorrELossFit::ELossFit*
AliFMDCorrELossFit::GetFit(UShort_t  d, Char_t r, Int_t etabin) const
{
  // 
  // Get the fit corresponding to the specified parameters 
  // 
  // Parameters:
  //    d      Detector 
  //    r      Ring 
  //    etabin Eta bin (1 based)
  // 
  // Return:
  //    Fit parameters or null in case of problems 
  //
  TObjArray* ringArray = GetRingArray(d, r);
  if (!ringArray)                                   return 0;
  if (etabin <= 0 || etabin >= fEtaAxis.GetNbins()) return 0;
  if      (etabin > ringArray->GetEntriesFast())    return 0;
  else if (etabin >= ringArray->GetEntriesFast())   etabin--;
  else if (!ringArray->At(etabin))                  etabin++;
  return static_cast<ELossFit*>(ringArray->At(etabin));
}
//____________________________________________________________________
AliFMDCorrELossFit::ELossFit*
AliFMDCorrELossFit::GetFit(UShort_t  d, Char_t r, Double_t eta) const
{
  // 
  // Find the fit corresponding to the specified parameters 
  // 
  // Parameters:
  //    d   Detector 
  //    r   Ring 
  //    eta Eta value 
  // 
  // Return:
  //    Fit parameters or null in case of problems 
  //
  Int_t etabin = FindEtaBin(eta);
  return GetFit(d, r, etabin);
}

//____________________________________________________________________
AliFMDCorrELossFit::ELossFit*
AliFMDCorrELossFit::FindFit(UShort_t  d, Char_t r, Int_t etabin,
			    UShort_t minQ) const
{
  // 
  // Find the fit corresponding to the specified parameters 
  // 
  // Parameters:
  //    d      Detector 
  //    r      Ring 
  //    etabin Eta bin (1 based)
  // 
  // Return:
  //    Fit parameters or null in case of problems 
  //
  if (etabin <= 0 || etabin >= fEtaAxis.GetNbins()) { 
    // AliError(Form("Eta bin=%3d out of bounds [%d,%d] for FMD%d%c", 
    //  	     etabin, 1, fEtaAxis.GetNbins(), d, r));
    return 0;
  }

  TObjArray* ringArray = GetRingArray(d, r);
  if (!ringArray) { 
    AliError(Form("Failed to make ring array for FMD%d%c", d, r));
    return 0;
  }
  DMSG(fDebug, 10, "Got ringArray %s for FMD%d%c", ringArray->GetName(), d, r);
  if (fCache.GetSize() <= 0) CacheBins(minQ);
#if 0
  Int_t idx = (d == 1 ? 0 : 
	       (d - 2) * 2 + 1 + (r=='I' || r=='i' ? 0 : 1));
#endif
  Int_t bin = CACHEDR(etabin, d, r, fEtaAxis.GetNbins());
  
  if (bin < 0 || bin > ringArray->GetEntriesFast()) return 0;
  
  return static_cast<ELossFit*>(ringArray->At(bin));
}
//____________________________________________________________________
AliFMDCorrELossFit::ELossFit*
AliFMDCorrELossFit::FindFit(UShort_t  d, Char_t r, Double_t eta,
			    UShort_t minQ) const
{
  // 
  // Find the fit corresponding to the specified parameters 
  // 
  // Parameters:
  //    d   Detector 
  //    r   Ring 
  //    eta Eta value 
  // 
  // Return:
  //    Fit parameters or null in case of problems 
  //
  Int_t etabin = FindEtaBin(eta);
  return FindFit(d, r, etabin, minQ);
}
//____________________________________________________________________
Int_t
AliFMDCorrELossFit::GetRingIndex(UShort_t d, Char_t r) const
{
  switch (d) { 
  case 1:   return 0;
  case 2:   return (r == 'i' || r == 'I') ? 1 : 2;
  case 3:   return (r == 'i' || r == 'I') ? 3 : 4;
  }
  return -1;
}
  
//____________________________________________________________________
TObjArray*
AliFMDCorrELossFit::GetRingArray(UShort_t d, Char_t r) const
{
  // 
  // Get the ring array corresponding to the specified ring
  // 
  // Parameters:
  //    d Detector 
  //    r Ring 
  // 
  // Return:
  //    Pointer to ring array, or null in case of problems
  //
  Int_t idx = GetRingIndex(d, r);
  if (idx < 0) return 0;
  return static_cast<TObjArray*>(fRings.At(idx));
}
//____________________________________________________________________
TObjArray*
AliFMDCorrELossFit::GetOrMakeRingArray(UShort_t d, Char_t r)
{
  // 
  // Get the ring array corresponding to the specified ring
  // 
  // Parameters:
  //    d Detector 
  //    r Ring 
  // 
  // Return:
  //    Pointer to ring array, or newly created container 
  //
  Int_t idx = GetRingIndex(d, r);
  if (idx < 0) return 0;
  if (idx >= fRings.GetEntriesFast() || !fRings.At(idx)) {
    TObjArray* a = new TObjArray(0);
    // TOrdCollection* a = new TOrdCollection(fEtaAxis.GetNbins());
    a->SetName(Form("FMD%d%c", d, r));
    a->SetOwner();
    fRings.AddAtAndExpand(a, idx);
  }
  return static_cast<TObjArray*>(fRings.At(idx));
}

//____________________________________________________________________
Double_t
AliFMDCorrELossFit::GetMpvCut(UShort_t  d,
			     Char_t    r,
			     Int_t     etabin,
			     Double_t  f) const
{
  ELossFit* fit = FindFit(d, r, etabin, 20);
  if (!fit) return -1024;
  return fit->GetMpvCut(f);
}

//____________________________________________________________________
Double_t
AliFMDCorrELossFit::GetXiCut(UShort_t  d,
			     Char_t    r,
			     Int_t     etabin,
			     Double_t  f) const
{
  ELossFit* fit = FindFit(d, r, etabin, 20);
  if (!fit) return -1024;
  return fit->GetXiCut(f);
}
//____________________________________________________________________
Double_t
AliFMDCorrELossFit::GetXiSigmaCut(UShort_t  d,
				  Char_t    r,
				  Int_t     etabin,
				  Double_t  f) const
{
  ELossFit* fit = FindFit(d, r, etabin, 20);
  if (!fit) return -1024;
  return fit->GetXiSigmaCut(f);
}
//____________________________________________________________________
Double_t
AliFMDCorrELossFit::GetAvgXiSigmaCut(UShort_t  d,
				     Char_t    r,
				     Int_t     etabin,
				     Double_t  f) const
{
  ELossFit* fit = FindFit(d, r, etabin, 20);
  if (!fit) return -1024;
  return fit->GetAvgXiSigmaCut(f);
}
//____________________________________________________________________
Double_t
AliFMDCorrELossFit::GetProbabilityCut(UShort_t  d,
				      Char_t    r,
				      Int_t     etabin,
				      Double_t  f) const
{
  ELossFit* fit = FindFit(d, r, etabin, 20);
  if (!fit) return -1024;
  return fit->FindProbabilityCut(f);
}

//____________________________________________________________________
Bool_t
AliFMDCorrELossFit::IsGood(Bool_t   verbose,
			   Double_t minRate, 
			   Int_t    maxGap,			   
			   Int_t    minInner,
			   Int_t    minOuter,
			   Int_t    minQuality)
{
  if (TestBit(kIsGoodAsserted)) return TestBit(kIsGood);

  Bool_t     isGood    = true; // Assume success  

  if (fRings.GetEntries() <= 0) isGood = false;
  TIter      nextRing(&fRings);
  TObjArray* ringArray = 0;
  while ((ringArray = static_cast<TObjArray*>(nextRing()))) {
    Char_t r = ringArray->GetName()[4];

    Int_t gap  = 0;
    Int_t good = 0;
    Int_t max  = 0;
    Int_t len  = 0;
    Bool_t first = true;
    for (Int_t iEta = 1; iEta <= fEtaAxis.GetNbins(); iEta++) {
      if (iEta > ringArray->GetEntriesFast()) {
	max = TMath::Max(gap, max);
	break;
      }
      TObject* o = ringArray->At(iEta);
      if (!o) {
	if (!first) { gap++; len++; }
	continue;
      }

      len++;
      first = false;
      ELossFit* fit = static_cast<ELossFit*>(o);
      if (fit->GetQuality() < minQuality) {
	gap++;
	continue;
      }
      good++;
      max = TMath::Max(max, gap);
      gap = 0;
    }
    Bool_t   thisGood = true; // Local result, assume success 
    Int_t    min      = (r == 'I' ? minInner : minOuter);
    Double_t rate     = len > 1 ? Double_t(good)/(len-1) : 0;
    if (rate < minRate) thisGood = false;
    if (len  < min)     thisGood = false;
    if (max  > maxGap)  thisGood = false;
    if (!thisGood)      isGood   = false; // Global result 
    if (verbose) 
      Printf("%s: %2d/%2d=%5.1f%%%-2s%5.1f%% good (%d) fits (%-2s %2d) "
	     "max gap %2d (%-2s %2d) -> %s",
	     ringArray->GetName(), good, len, 100*rate,
	     (rate < minRate ? "<" : ">="), 100*minRate, 	     
	     minQuality,
	     (len <  min    ? "<" : ">="), min, max,
	     (max >  maxGap ? ">" : "<="), maxGap,
	     (thisGood ? "good" : "bad"));
	
  }
  SetBit(kIsGoodAsserted,true);
  SetBit(kIsGood,        isGood);
  
  return isGood;
}

				  
//____________________________________________________________________
namespace { 
  TH1D* MakeHist(const TAxis& axis, const char* name, const char* title, 
		 Int_t color)
  {
    TH1D* h = new TH1D(Form("%s_%s", name, title), 
		       Form("%s %s", name, title), 
		       axis.GetNbins(), axis.GetXmin(), axis.GetXmax());
    h->SetDirectory(0);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(color);
    h->SetMarkerSize(0.5);
    h->SetFillColor(color);
    h->SetFillStyle(3001);
    h->SetLineColor(color);
    return h;
  }
}

  

#define IDX2RING(I) (i == 0 || i == 1 || i == 3 ? 'I' : 'O')
#define IDX2DET(I)  (i == 0 ? 1 : (i == 1 || i == 2 ? 2 : 3))
//____________________________________________________________________
TList*
AliFMDCorrELossFit::GetStacks(Bool_t   err, 
			      Bool_t   rel, 
			      Bool_t   good, 
			      UShort_t maxN) const
{
  // Get a list of THStacks 
  Int_t nRings = fRings.GetEntriesFast();
  // Int_t nPad   = 6+maxN-1; // 7 regular params, and maxN-1 weights

  enum { 
    kChi2nu = 0, 
    kC      = 1, 
    kDelta  = 2, 
    kXi     = 3, 
    kSigma  = 4, 
    kN      = 5 
  };
  
  THStack* sChi2nu;
  THStack* sC;
  THStack* sDelta;
  THStack* sXi;
  THStack* sSigma;
  // THStack* sigman;
  THStack* n;
  TList* stacks = new TList;
  stacks->AddAt(sChi2nu= new THStack("chi2",   "#chi^{2}/#nu"), kChi2nu);
  stacks->AddAt(sC     = new THStack("c",       "C"),           kC);
  stacks->AddAt(sDelta = new THStack("delta",  "#Delta_{mp}"),  kDelta);
  stacks->AddAt(sXi    = new THStack("xi",     "#xi"),          kXi);
  stacks->AddAt(sSigma = new THStack("sigma",  "#sigma"),       kSigma);
  //stacks->AddAt(sigman= new THStack("sigman", "#sigma_{n}"),   5);
  stacks->AddAt(n      = new THStack("n",      "N"),            kN);
  if (rel) { 
    sChi2nu->SetName("qual");
    sChi2nu->SetTitle("Quality");
    n->SetName("good");
    n->SetTitle("Bin map");
  }
  for (Int_t i = 1; i <= maxN; i++) {
    stacks->AddAt(new THStack(Form("a_%02d", i+1), Form("a_{%d}", i+1)), kN+i);
  }
  
  // TArrayD min(nPad);
  // TArrayD max(nPad);
  // min.Reset(100000);
  // max.Reset(-100000);

  for (Int_t i = 0; i < nRings; i++) { 
    if (!fRings.At(i)) continue;
    TObjArray* a     = static_cast<TObjArray*>(fRings.At(i));
    Int_t      nFits = a->GetEntriesFast();
    UShort_t   d     = IDX2DET(i);
    Char_t     r     = IDX2RING(i);
    Int_t      color = AliForwardUtil::RingColor(d, r);

    TH1* hChi    = MakeHist(fEtaAxis,a->GetName(), "chi2",   color);
    TH1* hC      = MakeHist(fEtaAxis,a->GetName(), "c",      color);
    TH1* hDelta  = MakeHist(fEtaAxis,a->GetName(), "delta",  color);
    TH1* hXi     = MakeHist(fEtaAxis,a->GetName(), "xi",     color);
    TH1* hSigma  = MakeHist(fEtaAxis,a->GetName(), "sigma",  color);
    // TH1D* hSigmaN = MakeHist(fEtaAxis,a->GetName(), "sigman", color);
    TH1* hN      = MakeHist(fEtaAxis,a->GetName(), "n",      color);
    TH1* hA[maxN];
    for (Int_t j = 0; j < maxN; j++) hA[j] = 0;
    const char* ho = (rel || !err ? "hist" : "e");
    sChi2nu->Add(hChi,   "hist"); // 0
    sC     ->Add(hC,     ho);     // 1
    sDelta ->Add(hDelta, ho);     // 2
    sXi    ->Add(hXi,    ho);     // 3
    sSigma ->Add(hSigma, ho);     // 4
    // sigman->Add(hSigmaN,ho);     // 5
    n     ->Add(hN,     "hist"); // 5
    hChi->SetFillColor(color);
    hChi->SetFillStyle(3001);
    hN->SetFillColor(color);
    hN->SetFillStyle(3001);

    for (Int_t k = 1; k <= maxN; k++) { 
      hA[k-1] = MakeHist(fEtaAxis,a->GetName(), Form("a%02d",k+1), color);
      static_cast<THStack*>(stacks->At(kN+k))->Add(hA[k-1], ho);
    }
    
    if (good) {
      for (Int_t j = 1; j <= fEtaAxis.GetNbins(); j++) {
	ELossFit* f = FindFit(d, r, j, 20);
	if (!f) continue;
	
	UpdateStackHist(f, rel, j, hChi, hN, hC, hDelta, hXi, hSigma, maxN, hA);
      }
    }
    else {
      for (Int_t j = 0; j < nFits; j++) {
	ELossFit* f = static_cast<ELossFit*>(a->At(j));
	if (!f) continue;
	
	UpdateStackHist(f, rel, CACHE(j,i,fEtaAxis.GetNbins()), 
			hChi, hN, hC, hDelta, hXi, hSigma, maxN, hA);
      }
    }
  }
  return stacks;
}
//____________________________________________________________________
void
AliFMDCorrELossFit::UpdateStackHist(ELossFit* f,    Bool_t rel, 
				    Int_t     used, 
				    TH1*      hChi, TH1*   hN, 
				    TH1*      hC,   TH1*   hDelta, 
				    TH1*      hXi,  TH1*   hSigma, 
				    Int_t     maxN, TH1**  hA) const
{
  Int_t    b      =f->fBin;
  Double_t chi2nu =(rel ? f->fQuality : (f->fNu <= 0 ? 0 : f->fChi2 / f->fNu));
  Double_t c      =(rel ? (f->fC     >0 ?f->fEC     /f->fC     :0) : f->fC);
  Double_t delta  =(rel ? (f->fDelta >0 ?f->fEDelta /f->fDelta :0) : f->fDelta);
  Double_t xi     =(rel ? (f->fXi    >0 ?f->fEXi    /f->fXi    :0) : f->fXi);
  Double_t sigma  =(rel ? (f->fSigma >0 ?f->fESigma /f->fSigma :0) : f->fSigma);
  Int_t    w      =(rel ? used : f->FindMaxWeight());
  // Double_t  sigman = (rel ? (f->fSigmaN>0 ?f->fESigmaN/f->fSigmaN :0) 
  //                     : f->SigmaN); 
  hChi   ->SetBinContent(b, chi2nu);
  hN     ->SetBinContent(b, w);
  hC     ->SetBinContent(b, c);
  hDelta ->SetBinContent(b, delta);
  hXi    ->SetBinContent(b, xi);
  hSigma ->SetBinContent(b, sigma);

  if (!rel) {
    hC     ->SetBinError(b, f->fEC);
    hDelta ->SetBinError(b, f->fEDelta);
    hXi    ->SetBinError(b, f->fEXi);
    hSigma ->SetBinError(b, f->fESigma);
    // hSigmaN->SetBinError(b, f->fESigmaN);
  }
  for (Int_t k = 0; k < f->fN-1 && k < maxN; k++) { 
    Double_t a = (rel ? (f->fA[k] > 0 ? f->fEA[k] / f->fA[k] : 0) : f->fA[k]);
    hA[k]->SetBinContent(b, a);
    if (!rel) hA[k]->SetBinError(b, f->fEA[k]);
  }
}

//____________________________________________________________________
void
AliFMDCorrELossFit::Draw(Option_t* option)
{
  // 
  // Draw this object 
  // 
  // Parameters:
  //    option Options.  Possible values are 
  //  - err Plot error bars 
  //
  TString opt(Form("nostack %s", option));
  opt.ToLower();
  Bool_t  rel = (opt.Contains("relative"));
  Bool_t  err = (opt.Contains("error"));
  Bool_t  clr = (opt.Contains("clear"));
  Bool_t  gdd = (opt.Contains("good"));
  if (rel) opt.ReplaceAll("relative","");
  if (err) opt.ReplaceAll("error","");
  if (clr) opt.ReplaceAll("clear", "");
  if (gdd) opt.ReplaceAll("good", "");

  UShort_t maxN   = 0;
  Int_t nRings = fRings.GetEntriesFast();
  for (Int_t i = 0; i < nRings; i++) { 
    if (!fRings.At(i)) continue;
    TObjArray* a     = static_cast<TObjArray*>(fRings.At(i));
    Int_t      nFits = a->GetEntriesFast();

    for (Int_t j = 0; j < nFits; j++) {
      ELossFit* fit = static_cast<ELossFit*>(a->At(j));
      if (!fit) continue;
      maxN          = TMath::Max(maxN, UShort_t(fit->fN));
    }
  }
  // AliInfo(Form("Maximum N is %d", maxN));
  Int_t nPad = 6+maxN-1; // 7 regular params, and maxN-1 weights
  TVirtualPad* pad = gPad;
  if (clr) { 
    pad->Clear();
    pad->SetTopMargin(0.02);
    pad->SetRightMargin(0.02);
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(0.10);
  }
  pad->Divide(2, (nPad+1)/2, 0.1, 0, 0);

  TList* stacks = GetStacks(err, rel, gdd, maxN);

  Int_t nPad2 = (nPad+1) / 2;
  for (Int_t i = 0; i < nPad; i++) {
    Int_t iPad = 1 + i/nPad2 + 2 * (i % nPad2);
    TVirtualPad* p = pad->cd(iPad);
    p->SetLeftMargin(.15);
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetGridx();
    p->SetGridy();
    if (rel && i != 0 && i != 6 && i != 5 && i != 4) p->SetLogy();


    THStack* stack = static_cast<THStack*>(stacks->At(i));
    if (!stack->GetHists() || stack->GetHists()->GetEntries() <= 0) { 
      AliWarningF("No histograms in %s", stack->GetName());
      continue;
    }
    // Double_t powMax = TMath::Log10(max[i]);
    // Double_t powMin = min[i] <= 0 ? powMax : TMath::Log10(min[i]);
    // if (powMax-powMin > 2. && min[i] != 0) p->SetLogy();

    // stack->SetMinimum(min[i]);
    // stack->SetMaximum(max[i]);
    stack->Draw(opt.Data());

    TString tit(stack->GetTitle());
    if (rel && i != 0 && i != 5)
      tit = Form("#delta %s/%s", tit.Data(), tit.Data());
    TH1*   hist  = stack->GetHistogram();
    TAxis* yaxis = hist->GetYaxis();
    yaxis->SetTitle(tit.Data());
    yaxis->SetTitleSize(0.15);
    yaxis->SetLabelSize(0.08);
    yaxis->SetTitleOffset(0.35);
    yaxis->SetTitleFont(132);
    yaxis->SetLabelFont(132);
    yaxis->SetNdivisions(5);


    TAxis* xaxis = stack->GetHistogram()->GetXaxis();
    xaxis->SetTitle("#eta");
    xaxis->SetTitleSize(0.15);
    xaxis->SetLabelSize(0.08);
    xaxis->SetTitleOffset(0.35);
    xaxis->SetTitleFont(132);
    xaxis->SetLabelFont(132);
    xaxis->SetNdivisions(10);

    stack->Draw(opt.Data());
  }
  pad->cd();      
}

//____________________________________________________________________
void
AliFMDCorrELossFit::Print(Option_t* option) const
{
  // 
  // Print this object.  
  // 
  // Parameters:
  //    option Options 
  //   - R   Print recursive  
  //
  //
  TString opt(option);
  opt.ToUpper();
  Int_t nRings  = fRings.GetEntriesFast();
  bool  recurse = opt.Contains("R");
  bool  cache   = opt.Contains("C") && fCache.GetSize() > 0;
  Int_t nBins   = fEtaAxis.GetNbins();

  std::cout << "Low cut in fit range: " << fLowCut << "\n"
	    << "Eta axis:             " << nBins
	    << " bins, range [" << fEtaAxis.GetXmin() << "," 
	    << fEtaAxis.GetXmax() << "]" << std::endl;
  
  for (Int_t i = 0; i < nRings; i++) { 
    if (!fRings.At(i)) continue;
    TObjArray* a     = static_cast<TObjArray*>(fRings.At(i));
    Int_t      nFits = a->GetEntriesFast();

    std::cout << a->GetName() << " [" << nFits << " entries]" 
	      << (recurse ? ":\n" : "\t");
    Int_t min = fEtaAxis.GetNbins()+1;
    Int_t max = 0;
    for (Int_t j = 0; j < nFits; j++) {
      if (!a->At(j)) continue;
      
      min = TMath::Min(j, min);
      max = TMath::Max(j, max);

      if (recurse) {
	std::cout << "Bin # " << j << "\t";
	ELossFit* fit = static_cast<ELossFit*>(a->At(j));
	fit->Print(option);
      }
    }
    if (!recurse) 
      std::cout << " bin range: " << std::setw(3) << min 
		<< "-" << std::setw(3) << max << " " << std::setw(3) 
		<< (max-min+1) << " bins" << std::endl;
  }

  if (!cache) return;

  std::cout << " eta bin           |           Fit bin              \n"
	    << " #       range     | FMD1i  FMD2i  FMD2o  FMD3i  FMD3o"
    // << "----+-----+++------+-----------------------------------"
	    << std::endl;
  size_t oldPrec = std::cout.precision();
  std::cout.precision(3);
  for (Int_t i = 1; i <= nBins; i++) { 
    std::cout << std::setw(4) << i << " " 
	      << std::setw(5) << std::showpos << fEtaAxis.GetBinLowEdge(i)
	      << " - " << std::setw(5) << fEtaAxis.GetBinUpEdge(i) 
	      << std::noshowpos << " | ";
    for (Int_t j = 0; j < 5; j++) {
      Int_t bin = CACHE(i,j,nBins);
      if (bin <= 0) std::cout << "       ";
      else          std::cout << std::setw(5) << bin 
			      << (bin == i ? ' ' : '*') << ' ';
    }
    std::cout << std::endl;
  }
  std::cout.precision(oldPrec);
}

//____________________________________________________________________
void
AliFMDCorrELossFit::Browse(TBrowser* b)
{
  // 
  // Browse this object 
  // 
  // Parameters:
  //    b 
  //
  b->Add(&fRings);
  b->Add(&fEtaAxis);
}



//____________________________________________________________________
//
// EOF
//
