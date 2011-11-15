// Object holding the Energy loss fit 'correction' 
// 
// These are generated from Monte-Carlo or real ESDs. 
//
#include "AliFMDCorrELossFit.h"
#include "AliForwardUtil.h"
#include <TF1.h>
#include <TBrowser.h>
#include <TVirtualPad.h>
#include <THStack.h>
#include <TH1D.h>
#include <AliLog.h>
#include <TMath.h>
#include <TList.h>
#include <iostream>
#include <iomanip>

Double_t AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
Double_t AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
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
    fBin(0)
{
  //
  // Default constructor 
  // 
  //
}
//____________________________________________________________________
AliFMDCorrELossFit::ELossFit::ELossFit(Int_t quality, const TF1& f)
  : fN(f.GetNpar() > AliForwardUtil::ELossFitter::kN ? 
       Int_t(f.GetParameter(AliForwardUtil::ELossFitter::kN)) : 
       1),
    fNu(f.GetNDF()),
    fChi2(f.GetChisquare()),
    fC(f.GetParameter(AliForwardUtil::ELossFitter::kC)),
    fDelta(f.GetParameter(AliForwardUtil::ELossFitter::kDelta)),
    fXi(f.GetParameter(AliForwardUtil::ELossFitter::kXi)),
    fSigma(f.GetParameter(AliForwardUtil::ELossFitter::kSigma)),
    fSigmaN(f.GetParameter(AliForwardUtil::ELossFitter::kSigmaN)),
    fA(0),
    fEC(f.GetParError(AliForwardUtil::ELossFitter::kC)),
    fEDelta(f.GetParError(AliForwardUtil::ELossFitter::kDelta)),
    fEXi(f.GetParError(AliForwardUtil::ELossFitter::kXi)),
    fESigma(f.GetParError(AliForwardUtil::ELossFitter::kSigma)),
    fESigmaN(f.GetParError(AliForwardUtil::ELossFitter::kSigmaN)),
    fEA(0),
    fQuality(quality),
    fDet(0), 
    fRing('\0'),
    fBin(0)
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
    fA[i]  = f.GetParameter(AliForwardUtil::ELossFitter::kA+i);
    fEA[i] = f.GetParError(AliForwardUtil::ELossFitter::kA+i);
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
    fBin(0)
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
    fBin(o.fBin)
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
  fN	   = o.fN;
  fNu	   = o.fNu;
  fChi2	   = o.fChi2;
  fC	   = o.fC;
  fDelta   = o.fDelta;
  fXi	   = o.fXi;
  fSigma   = o.fSigma;
  fSigmaN  = o.fSigmaN;
  fEC	   = o.fEC;
  fEDelta  = o.fEDelta;
  fEXi	   = o.fEXi;
  fESigma  = o.fESigma;
  fESigmaN = o.fESigmaN;
  fQuality = o.fQuality;
  fDet     = o.fDet; 
  fRing    = o.fRing;
  fBin     = o.fBin;
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
					    UShort_t  maxN) const
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
  Int_t n = TMath::Min(maxN, UShort_t(fN-1));
  Int_t m = 1;
  // fN is one larger than we have data 
  for (Int_t i = 0; i < n-1; i++, m++) { 
    if (fA[i] < leastWeight)  break;
    if (fEA[i] / fA[i] > maxRelError) break;
  }
  return m;
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
  // (see AliForwardUtil::NLandauGaus) for the maximum @f$ N @f$
  // that fulfills the requirements 
  // 
  // Parameters:
  //    x           Where to evaluate 
  //    maxN 	  @f$ \max{N}@f$    
  // 
  // Return:
  //    @f$ f_N(x;\Delta,\xi,\sigma')@f$ 
  //
  return AliForwardUtil::NLandauGaus(x, fDelta, fXi, fSigma, fSigmaN, 
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
  // See also AliForwardUtil::ILandauGaus and AliForwardUtil::NLandauGaus
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
    Double_t f = AliForwardUtil::ILandauGaus(x,fDelta,fXi,fSigma,fSigmaN,i);
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
  return 0;
}

//____________________________________________________________________
void
AliFMDCorrELossFit::ELossFit::Print(Option_t*) const
{
  // 
  // Information to standard output 
  // 
  // Parameters:
  //    option Not used 
  //
  std::cout << GetName() << ":\n"
	    << " chi^2/nu = " << fChi2 << "/" << fNu << " = " 
	    << (fNu == 0 ? 999 : fChi2 / fNu) << "\n"
	    << "     Quality:   " << fQuality << "\n" 
	    << "  NParticles:   " << fN << "  (" << FindMaxWeight() << ")\n"
	    << OUTPAR("Delta", fDelta, fEDelta) 
	    << OUTPAR("xi", fXi, fEXi)
	    << OUTPAR("sigma", fSigma, fESigma)
	    << OUTPAR("sigma_n", fSigmaN, fESigmaN);
  for (Int_t i = 0; i < fN-1; i++) 
    std::cout << OUTPAR(Form("a%d", i+2), fA[i], fEA[i]);
  std::cout << std::flush;
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
  Draw(b ? b->GetDrawOption() : "comp");
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
  if (opt.Contains("COMP")) { 
    opt.ReplaceAll("COMP","");
    comp = true;
  }
  if (!opt.Contains("SAME")) { 
    gPad->Clear();
  }

  TObjArray cleanup;
  TF1* tot = AliForwardUtil::MakeNLandauGaus(1, 
					     fDelta, fXi, 
					     fSigma, fSigmaN, 
					     fN,     fA, 
					     0.01,   10);
  tot->SetLineColor(kBlack);
  tot->SetLineWidth(2);
  tot->SetLineStyle(1);
  tot->SetTitle(GetName());
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
  tot->DrawCopy(opt.Data());
  cleanup.Add(tot);

  if (!comp) { 
    gPad->cd();
    return;
  }

  Double_t min = max;
  opt.Append(" same");
  Int_t maxW = FindMaxWeight();
  for (Int_t i=1; i <= fN; i++) { 
    TF1* f = AliForwardUtil::MakeILandauGaus((i == 1 ? 1 : fA[i-2]), 
					     fDelta, fXi, 
					     fSigma, fSigmaN, 
					     i,      0.01, 10);
    f->SetLineWidth(2);
    f->SetLineStyle(i > maxW ? 2 : 1);
    min = TMath::Min(f->GetMaximum(), min);
    f->DrawCopy(opt.Data());
    cleanup.Add(f);
  }
  min /= 100;
  tot->GetHistogram()->SetMaximum(max);
  tot->GetHistogram()->SetMinimum(min);
  tot->GetHistogram()->GetYaxis()->SetRangeUser(min, max);

  gPad->cd();
}

						 
//____________________________________________________________________
#define CHECKPAR(V,E,T) ((V > 0) && (E / V < T))

//____________________________________________________________________
Double_t
AliFMDCorrELossFit::ELossFit::GetLowerBound(Double_t f) const
{
  // 
  // Return 
  //    Delta * f
  return f * fDelta;
}
//____________________________________________________________________
Double_t
AliFMDCorrELossFit::ELossFit::GetLowerBound(Double_t f, 
					    Bool_t includeSigma) const
{
  // 
  // Return 
  //    Delta - f * (xi + sigma)
  return fDelta - f * (fXi + (includeSigma ? fSigma : 0));
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
  Int_t qual = 0;
  if (fNu > 0 && fChi2 / fNu < maxChi2nu) qual += 4;;
  if (CHECKPAR(fDelta,  fEDelta,  maxRelError)) qual++;
  if (CHECKPAR(fXi,     fEXi,     maxRelError)) qual++;
  if (CHECKPAR(fSigma,  fESigma,  maxRelError)) qual++;
  if (CHECKPAR(fSigmaN, fESigmaN, maxRelError)) qual++;
  qual += FindMaxWeight(1.5*maxRelError, leastWeight, fN);
  fQuality = qual;
}

//____________________________________________________________________
AliFMDCorrELossFit::AliFMDCorrELossFit()
  : TObject(), 
    fRings(), 
    fEtaAxis(0,0,0), 
    fLowCut(0)
{
  // 
  // Default constructor 
  //
  fRings.SetOwner(kTRUE);
  fEtaAxis.SetTitle("#eta");
  fEtaAxis.SetName("etaAxis");
  fRings.SetName("rings");
}

//____________________________________________________________________
AliFMDCorrELossFit::AliFMDCorrELossFit(const AliFMDCorrELossFit& o)
  : TObject(o), 
    fRings(o.fRings),
    fEtaAxis(o.fEtaAxis.GetNbins(),o.fEtaAxis.GetXmin(),o.fEtaAxis.GetXmax()), 
    fLowCut(0)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  fEtaAxis.SetTitle("#eta");
  fEtaAxis.SetName("etaAxis");
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
  SetEtaAxis(o.fEtaAxis.GetNbins(), o.fEtaAxis.GetXmin(), o.fEtaAxis.GetXmax());

  return *this;
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
AliFMDCorrELossFit::FindFit(UShort_t  d, Char_t r, Int_t etabin) const
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
  TObjArray* ringArray = GetRingArray(d, r);
  if (!ringArray) { 
    AliError(Form("Failed to make ring array for FMD%d%c", d, r));
    return 0;
  }
  if (etabin <= 0 || etabin >= fEtaAxis.GetNbins()) { 
    // AliError(Form("Eta bin=%3d out of bounds [%d,%d] for FMD%d%c", 
    //  	     etabin, 1, fEtaAxis.GetNbins(), d, r));
    return 0;
  }
  if (etabin > ringArray->GetEntriesFast()) { 
    // AliError(Form("Eta bin=%3d out of bounds [%d,%d] for FMD%d%c", 
    // 		     etabin, 1, ringArray->GetEntriesFast(), d, r));
    return 0;
  }
  else if (etabin >= ringArray->GetEntriesFast()) { 
    // AliWarning(Form("Eta bin=%3d out of bounds by +1 [%d,%d] for FMD%d%c, " 
    //		    "trying %3d", etabin, 1, ringArray->GetEntriesFast(), d, r,
    //		    etabin-1));
    etabin--;
  }
  else if (!ringArray->At(etabin)) { 
    // AliWarning(Form("Eta bin=%d has no fit for FMD%d%c, trying %03d", 
    // 		    etabin, d, r, etabin+1));
    etabin++;
  }
  return static_cast<ELossFit*>(ringArray->At(etabin));
}
//____________________________________________________________________
AliFMDCorrELossFit::ELossFit*
AliFMDCorrELossFit::FindFit(UShort_t  d, Char_t r, Double_t eta) const
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
  return FindFit(d, r, etabin);
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
  Int_t idx = -1;
  switch (d) { 
  case 1:   idx = 0; break;
  case 2:   idx = (r == 'i' || r == 'I') ? 1 : 2; break;
  case 3:   idx = (r == 'o' || r == 'I') ? 3 : 4; break;
  }
  if (idx < 0 || idx >= fRings.GetEntriesFast()) return 0;
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
  Int_t idx = -1;
  switch (d) { 
  case 1:   idx = 0; break;
  case 2:   idx = (r == 'i' || r == 'I') ? 1 : 2; break;
  case 3:   idx = (r == 'o' || r == 'I') ? 3 : 4; break;
  }
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
AliFMDCorrELossFit::GetLowerBound(UShort_t  d, Char_t r, Int_t etabin,
				  Double_t f) const
{
  ELossFit* fit = GetFit(d, r, etabin);
  if (!fit) return -1024;
  return fit->GetLowerBound(f);
}
//____________________________________________________________________
Double_t
AliFMDCorrELossFit::GetLowerBound(UShort_t  d, Char_t r, Double_t eta,
				  Double_t f) const
{
  Int_t bin = FindEtaBin(eta);
  if (bin <= 0) return -1024;
  return GetLowerBound(d, r, Int_t(bin), f);
}
//____________________________________________________________________
Double_t
AliFMDCorrELossFit::GetLowerBound(UShort_t  d, Char_t r, Int_t etabin,
				  Double_t f, Bool_t showErrors, 
				  Bool_t includeSigma) const
{
  ELossFit* fit = GetFit(d, r, etabin);
  if (!fit) { 
    if (showErrors) {
      AliWarning(Form("No fit for FMD%d%c @ etabin=%d", d, r, etabin));
    }
    return -1024;
  }
  return fit->GetLowerBound(f, includeSigma);
}

//____________________________________________________________________
Double_t
AliFMDCorrELossFit::GetLowerBound(UShort_t  d, Char_t r, Double_t eta,
				  Double_t f, Bool_t showErrors, 
				  Bool_t includeSigma) const
{
  Int_t bin = FindEtaBin(eta);
  if (bin <= 0) { 
    if (showErrors)
      AliError(Form("eta=%f out of bounds for FMD%d%c", eta, d, r));
    return -1024;
  }
  return GetLowerBound(d, r, bin, f, showErrors, includeSigma);
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
AliFMDCorrELossFit::GetStacks(Bool_t err, Bool_t rel, UShort_t maxN) const
{
  // Get a list of THStacks 
  Int_t nRings = fRings.GetEntriesFast();
  Int_t nPad   = 6+maxN-1; // 7 regular params, and maxN-1 weights

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
  stacks->AddAt(n     = new THStack("n",      "N"),            kN);
  for (Int_t i = 1; i <= maxN; i++) {
    stacks->AddAt(new THStack(Form("a_%02d", i+1), Form("a_{%d}", i+1)), kN+i);
  }
  
  TArrayD min(nPad);
  TArrayD max(nPad);
  min.Reset(100000);
  max.Reset(-100000);

  for (Int_t i = 0; i < nRings; i++) { 
    if (!fRings.At(i)) continue;
    TObjArray* a     = static_cast<TObjArray*>(fRings.At(i));
    Int_t      nFits = a->GetEntriesFast();
    Int_t      color = AliForwardUtil::RingColor(IDX2DET(i), IDX2RING(i));

    TH1D* hChi    = MakeHist(fEtaAxis,a->GetName(), "chi2",   color);
    TH1D* hC      = MakeHist(fEtaAxis,a->GetName(), "c",      color);
    TH1D* hDelta  = MakeHist(fEtaAxis,a->GetName(), "delta",  color);
    TH1D* hXi     = MakeHist(fEtaAxis,a->GetName(), "xi",     color);
    TH1D* hSigma  = MakeHist(fEtaAxis,a->GetName(), "sigma",  color);
    // TH1D* hSigmaN = MakeHist(fEtaAxis,a->GetName(), "sigman", color);
    TH1D* hN      = MakeHist(fEtaAxis,a->GetName(), "n",      color);
    TH1D* hA[maxN];
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
			  
    for (Int_t j = 0; j < nFits; j++) {
      ELossFit* f = static_cast<ELossFit*>(a->At(j));
      if (!f) continue;

      Int_t     b      = f->fBin;
      Int_t     nW     = f->FindMaxWeight();
      Double_t  vChi2nu = (f->fNu <= 0 ? 0 : f->fChi2 / f->fNu);
      Double_t  vC      = (rel ? (f->fC     >0 ?f->fEC     /f->fC      :0) 
			   : f->fC);
      Double_t  vDelta  = (rel ? (f->fDelta >0 ?f->fEDelta /f->fDelta  :0) 
			   : f->fDelta);
      Double_t  vXi     = (rel ? (f->fXi    >0 ?f->fEXi    /f->fXi     :0) 
			   : f->fXi);
      Double_t  vSigma  = (rel ? (f->fSigma >0 ?f->fESigma /f->fSigma  :0) 
			   : f->fSigma);
      // Double_t  sigman = (rel ? (f->fSigmaN>0 ?f->fESigmaN/f->fSigmaN :0) 
      //                     : f->SigmaN); 
      hChi   ->SetBinContent(b, vChi2nu);
      hN     ->SetBinContent(b, nW);
      hC     ->SetBinContent(b, vC);
      hDelta ->SetBinContent(b, vDelta);
      hXi    ->SetBinContent(b, vXi);
      hSigma ->SetBinContent(b, vSigma);
      if (vChi2nu > 1e-12) {
	min[kChi2nu] = TMath::Min(min[kChi2nu], vChi2nu);
	max[kChi2nu] = TMath::Max(max[kChi2nu], vChi2nu);
      }
      if (vC > 1e-12) {
	min[kC] = TMath::Min(min[kC], vC);
	max[kC] = TMath::Max(max[kC], vC);
      }
      if (vDelta > 1e-12) {
	min[kDelta] = TMath::Min(min[kDelta], vDelta);
	max[kDelta] = TMath::Max(max[kDelta], vDelta);
      }
      if (vXi > 1e-12) {
	min[kXi] = TMath::Min(min[kXi], vXi);
	max[kXi] = TMath::Max(max[kXi], vXi);
      }
      if (vSigma > 1e-12) {
	min[kSigma] = TMath::Min(min[kSigma], vSigma);
	max[kSigma] = TMath::Max(max[kSigma], vSigma);
      }
      if (nW > 1e-12) { 
	min[kN] = TMath::Min(min[kN], Double_t(nW));
	max[kN] = TMath::Max(max[kN], Double_t(nW));
      }
      // hSigmaN->SetBinContent(b,  sigman);
      if (!rel) {
	hC     ->SetBinError(b, f->fEC);
	hDelta ->SetBinError(b, f->fEDelta);
	hXi    ->SetBinError(b, f->fEXi);
	hSigma ->SetBinError(b, f->fESigma);
	// hSigmaN->SetBinError(b, f->fESigmaN);
      }
      for (Int_t k = 0; k < f->fN-1 && k < maxN; k++) { 
	Double_t vA = (rel ? (f->fA[k] > 0 ? f->fEA[k] / f->fA[k] : 0) 
		       : f->fA[k]);
	hA[k]->SetBinContent(b, vA);
	if (!rel) hA[k]->SetBinError(b, f->fEA[k]);
	if (vA > 1e-12) {
	  min[kN+1+k] = TMath::Min(min[kN+1+k], vA);
	  max[kN+1+k] = TMath::Max(max[kN+1+k], vA);
	}
      }
    }
  }
  return stacks;
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
  if (rel) opt.ReplaceAll("relative","");
  if (err) opt.ReplaceAll("error","");

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
  pad->Divide(2, (nPad+1)/2, 0.1, 0, 0);

  TList* stacks = GetStacks(err, rel, maxN);

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

    // Double_t powMax = TMath::Log10(max[i]);
    // Double_t powMin = min[i] <= 0 ? powMax : TMath::Log10(min[i]);
    // if (powMax-powMin > 2. && min[i] != 0) p->SetLogy();

    // stack->SetMinimum(min[i]);
    // stack->SetMaximum(max[i]);
    stack->Draw(opt.Data());

    TString tit(stack->GetTitle());
    if (rel && i != 0 && i != 6)
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
  Int_t nRings = fRings.GetEntriesFast();
  bool recurse = opt.Contains("R");

  std::cout << "Low cut in fit range: " << fLowCut << "\n"
	    << "Eta axis:             " << fEtaAxis.GetNbins() 
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
