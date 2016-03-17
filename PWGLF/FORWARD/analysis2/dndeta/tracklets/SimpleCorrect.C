#ifndef __CINT__
# include <TCanvas.h>
# include <TList.h>
# include <TH1.h>
# include <TH2.h>
# include <THStack.h>
# include <TFile.h>
# include <TError.h>
# include <TLatex.h>
# include <TFitResult.h>
# include <TF1.h>
# include <TMath.h>
# include <THashList.h>
#else
class TList;
class TSeqCollection;
class TH1;
class TH2;
class THStack;
class TAxis;
class TDirectory;
#endif

class TCanvas; // Auto load 

const Bool_t kSimpleCorrectLoaded = true;

/**
 * dN/dy axis title 
 */
const char* dndetaTitle = "\\mathrm{d}N_{\\mathrm{ch}}/\\mathrm{d}\\eta";
/** 
 * Color scale 
 */
const Color_t cc[] = { kMagenta+2,
		       kBlue+2,
		       kAzure-1, // 10,
		       kCyan+2,
		       kGreen+1,
		       kSpring+5,//+10,
		       kYellow+1,
		       kOrange+5,//+10,
		       kRed+1,
		       kPink+5,//+10,
		       kBlack };

//====================================================================
TObject* GetO(TSeqCollection* l, const char* name, TClass* cls)
{
  if (!l) {
    Warning("GetO", "No list passed");
    return 0;
  }
  TObject* o = l->FindObject(name);
  if (!o) {
    Warning("GetO", "No object named %s found in list %s",
	    name, l->GetName());
    return 0;
  }
  if (!cls) return o;

  if (!o->IsA()->InheritsFrom(cls)) {
    Warning("GetO", "Object %s read from %s is not a %s but a %s",
	    name, l->GetName(), cls->GetName(), o->ClassName());
    return 0;
  }
  return o;
}
//____________________________________________________________________
TH1* GetH1(TSeqCollection* l, const char* name)
{
  return static_cast<TH1*>(GetO(l,name,TH1::Class()));
}
//____________________________________________________________________
TH2* GetH2(TSeqCollection* l, const char* name)
{
  return static_cast<TH2*>(GetO(l,name,TH2::Class()));
}

//====================================================================
Bool_t CheckAxisNBins(const char*  which,
		      const TAxis* a1,
		      const TAxis* a2)
{
  if (a1->GetNbins() != a2->GetNbins()) {
    ::Warning("CheckAxisNBins", "Incompatible number %s bins: %d vs %d",
	      which, a1->GetNbins(), a2->GetNbins());
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t CheckAxisLimits(const char*  which,
		       const TAxis* a1,
		       const TAxis* a2)
{
  if (!TMath::AreEqualRel(a1->GetXmin(), a2->GetXmin(),1.E-12) ||
      !TMath::AreEqualRel(a1->GetXmax(), a2->GetXmax(),1.E-12)) {
    Warning("CheckAxisLimits",
	    "Limits of %s axis incompatible [%f,%f] vs [%f,%f]", which,
	    a1->GetXmin(), a1->GetXmax(), a2->GetXmin(), a2->GetXmax());
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t CheckAxisBins(const char*  which,
		     const TAxis* a1,
		     const TAxis* a2)
{
  const TArrayD * h1Array = a1->GetXbins();
  const TArrayD * h2Array = a2->GetXbins();
  Int_t fN = h1Array->fN;
  if ( fN == 0 ) return true;
  if (h2Array->fN != fN) {
    // Redundant?
    Warning("CheckAxisBins", "Not equal number of %s bin limits: %d vs %d",
	    which, fN, h2Array->fN);
    return false;
  }
  else {
    for (int i = 0; i < fN; ++i) {
      if (!TMath::AreEqualRel(h1Array->GetAt(i),h2Array->GetAt(i),1E-10)) {
	Warning("CheckAxisBins",
		"%s limit # %3d incompatible: %f vs %f",
		which, i, h1Array->GetAt(i),h2Array->GetAt(i));
	return false;
      }
    }
  }
  return true;
}
//____________________________________________________________________
Bool_t CheckAxisLabels(const char*  which,
		       const TAxis* a1,
		       const TAxis* a2)
{
  // check that axis have same labels
  THashList *l1 = (const_cast<TAxis*>(a1))->GetLabels();
  THashList *l2 = (const_cast<TAxis*>(a2))->GetLabels();
  
  if (!l1 && !l2) return true;
  if (!l1 ||  !l2) {
    Warning("CheckAxisLabels", "Difference in %s labels: %p vs %p",
	    which, l1, l2);
    return false;
  }
  // check now labels sizes  are the same
  if (l1->GetSize() != l2->GetSize()) {
    Warning("CheckAxisLabels", "Different number of %s labels: %d vs %d",
	    which, l1->GetSize(), l2->GetSize());
    return false;
  }
  for (int i = 1; i <= a1->GetNbins(); ++i) {
    TString label1 = a1->GetBinLabel(i);
    TString label2 = a2->GetBinLabel(i);
    if (label1 != label2) {
      Warning("CheckAxisLabels", "%s label # %d not the same: '%s' vs '%s'",
	      which, i, label1.Data(), label2.Data());
      return false;
    }
  }
    
  return true;
}
//____________________________________________________________________
Bool_t CheckAxis(const char*  which, 
		 const TAxis* a1,
		 const TAxis* a2,
		 Bool_t       alsoLbls)
{
  if (!CheckAxisNBins (which, a1, a2)) return false;
  if (!CheckAxisLimits(which, a1, a2)) return false;
  if (!CheckAxisBins  (which, a1, a2)) return false;
  if (alsoLbls && !CheckAxisLabels(which, a1, a2)) return false;
  return true;
}
//____________________________________________________________________
Bool_t CheckConsistency(const TH1* h1, const TH1* h2)
{
  // Check histogram compatibility
  if (h1 == h2) return true;

  if (h1->GetDimension() != h2->GetDimension()) {
    Warning("CheckConsistency",
	    "%s and %s have different dimensions %d vs %d",
	    h1->GetName(), h2->GetName(),
	    h1->GetDimension(), h2->GetDimension());
    return false;
  }
  Int_t dim = h1->GetDimension(); 
    
  Bool_t ret = true;
  Bool_t alsoLbls = (h1->GetEntries() != 0 && h2->GetEntries() != 0);
  if (!CheckAxis("X", h1->GetXaxis(), h2->GetXaxis(), alsoLbls)) ret = false;
  if (dim > 1 &&
      !CheckAxis("Y", h1->GetYaxis(), h2->GetYaxis(), alsoLbls)) ret = false;
  if (dim > 2 &&
      !CheckAxis("Z", h1->GetZaxis(), h2->GetZaxis(), alsoLbls)) ret = false;
    
  return ret;
}

//____________________________________________________________________
TH1* Coerce(const TH1* templ, TH1* src)
{
  if (CheckConsistency(templ, src)) return src;
  TH1* ret = static_cast<TH1*>(templ->Clone(src->GetName()));
  ret->SetDirectory(0);
  ret->Reset();
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++) {
    Double_t x = ret->GetXaxis()->GetBinCenter(i);
    Int_t    b = src->GetXaxis()->FindBin(x);
    if (b <= 0 || b >= src->GetNbinsX()) continue;
    Double_t c = src->GetBinContent(b);
    Double_t e = src->GetBinError(b);
    ret->SetBinContent(i, c);
    ret->SetBinError  (i, e);
  }
  return ret;
}

//____________________________________________________________________
TH2* Coerce(const TH2* templ, TH2* src)
{
  if (CheckConsistency(templ, src)) return src;
  TH2* ret = static_cast<TH2*>(templ->Clone(src->GetName()));
  ret->SetDirectory(0);
  ret->Reset();
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++) {
    Double_t x  = ret->GetXaxis()->GetBinCenter(i);
    Int_t    bx = src->GetXaxis()->FindBin(x);
    if (bx <= 0 || bx >= src->GetNbinsX()) continue;
    for (Int_t j = 1; j <= ret->GetNbinsY(); j++) {
      Double_t y  = ret->GetYaxis()->GetBinCenter(j);
      Int_t    by = src->GetYaxis()->FindBin(y);
      Double_t c  = src->GetBinContent(bx,by);
      Double_t e  = src->GetBinError(bx,by);
      ret->SetBinContent(i, j, c);
      ret->SetBinError  (i, j, e);
    }
  }
  return ret;
}

//____________________________________________________________________
TH1* SetAttr(TH1* h,
	     Color_t  color,
	     Style_t  marker=20,
	     Double_t size=1.,
	     Style_t  fill=0,
	     Style_t  line=1,
	     Width_t  width=1)
{
  if (!h) return 0;
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(size);
  h->SetFillColor(color);
  h->SetFillStyle(fill);
  h->SetLineColor(color);
  h->SetLineStyle(line);
  h->SetLineWidth(width);
  h->GetXaxis()->SetNdivisions(210);
  h->GetYaxis()->SetNdivisions(210);
}

//____________________________________________________________________
void PrintH(TH2* h, Int_t prec=2)
{
  Printf("Content of %s - %s", h->GetName(), h->GetTitle());
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    printf("%3d: ", i);
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      Double_t c = h->GetBinContent(i,j);
      Double_t e = h->GetBinError  (i,j);
      if (TMath::IsNaN(c) || TMath::IsNaN(e))
	printf("*** NAN ***");
      else 
	printf("%.*f+/-%.*f ", prec, c, prec, e);
    }
    printf("\n");
  }
}
//____________________________________________________________________
void PrintH(TH1* h, Int_t prec=2)
{
  Printf("Content of %s - %s", h->GetName(), h->GetTitle());
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    if (h->GetBinContent(i) <= 1e-6) continue;
    printf("%3d (%+6.3f): %.*f+/-%.*f\n", i,
	   h->GetXaxis()->GetBinCenter(i),
	   prec, h->GetBinContent(i),
	   prec, h->GetBinError(i));
  }
}

//====================================================================
Double_t Integrate(TH1*      h,
		   Double_t  min,
		   Double_t  max,
		   Double_t& err)
{
  const Double_t epsilon = 1e-6;
  Int_t bMin = h->FindBin(min+epsilon);
  Int_t bMax = h->FindBin(max-epsilon);
  if (bMin < 1) bMin = 0; // Include underflow bin
  Double_t val = h->IntegralAndError(bMin, bMax, err);
  // For a-symmetric histograms, return 
  if (TMath::Abs(h->GetXaxis()->GetXmin()+h->GetXaxis()->GetXmax())>=epsilon)
    return val;

  // Otherwise, also integrate negative side
  Double_t err2;
  bMin =  h->FindBin(-min+epsilon);
  bMax =  h->FindBin(-max-epsilon);
  val  += h->IntegralAndError(bMin, bMax, err2);
  err  =  TMath::Sqrt(err*err+err2*err2);
  return val;
}
//____________________________________________________________________
Double_t RatioE(Double_t n, Double_t en,
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


//====================================================================
TH2* CorrectForIPz(TH2* in, TH1* ipz, Bool_t full=true)
{
  // Printf("Scaling to IP of %s", in->GetName());
  for (Int_t iy = 1; iy <= in->GetNbinsY(); iy++) {
    Double_t z   = in->GetYaxis()->GetBinCenter(iy);
    Int_t    bin = ipz->GetXaxis()->FindBin(z);
    Double_t nEv = ipz->GetBinContent(bin);
    Double_t eEv = ipz->GetBinError  (bin);
    Double_t esc = (nEv > 0 ? 1./nEv : 0);
    Double_t rE2 = esc*esc*eEv*eEv;
    // Printf("z=%f -> %f +/- %f events", z, nEv, eEv);
    for (Int_t ix = 1; ix <= in->GetNbinsX(); ix++) {
      Double_t c   = in->GetBinContent(ix,iy);
      Double_t e   = in->GetBinError  (ix,iy);
      Double_t r   = (c > 0 ? e/c : 0);
      // Scale by number of events, and error propagate 
      Double_t sc  = c * esc;
      Double_t se  = 0;
      if (full)se  = sc * TMath::Sqrt(r*r+rE2);
      else     se  = e * esc;
      Double_t scl = 1 / in->GetXaxis()->GetBinWidth(ix);
      // Printf("Setting bin (%3d,%3d) to %f +/- %f (%f %5.1f%% %5.1f%%)",
      //        ix,iy,scl*sc, scl*se, scl, 100*r, 100*esc*eEv);
      in->SetBinContent(ix, iy, scl*sc);
      in->SetBinError  (ix, iy, scl*se);
    }
  }
  return in;
}

//____________________________________________________________________
TH1* Avg(TH2*        h,
	 TH1*        ipz,
	 UShort_t    mode, 
	 const char* name,
	 TH2*        other=0) 
{
  if (!h) return 0;
  TH2* mask = other ? other : h;
  // Printf("Averaging %s - %s (%p) %s", h->GetName(), h->GetTitle(),
  //        ipz, (mode  ? "errors only" : "full scale"));
  Int_t nIPz = h->GetNbinsY();
  Int_t nEta = h->GetNbinsX();
  TH1*  p    = h->ProjectionX(name,1,nIPz,"e");
  p->SetDirectory(0);
  p->SetFillColor(0);
  p->SetFillStyle(0);
  p->SetYTitle(Form("\\langle(%s)\\rangle", h->GetYaxis()->GetTitle()));
  p->Reset();
  for (Int_t etaBin = 1; etaBin <= nEta; etaBin++) {
    TArrayD hv(nIPz);
    TArrayD he(nIPz);
    TArrayD hr(nIPz);
    TArrayI hb(nIPz);
    Int_t   j = 0;
    for (Int_t ipzBin = 1; ipzBin <= nIPz; ipzBin++) {
      Double_t bc = mask->GetBinContent(etaBin, ipzBin);
      if (bc < 1e-9) continue; // Low value
      Double_t be = mask->GetBinError  (etaBin, ipzBin);
      if (TMath::IsNaN(be)) continue; // Bad error value
      hv[j] = bc;
      he[j] = be;
      hr[j] = be/bc;
      hb[j] = ipzBin;
      j++;		
    }
    // Now we have all interesting values.  Sort the relative error
    // array to get the most significant first.
    TArrayI idx(nIPz);
    TMath::Sort(j, hr.fArray, idx.fArray, false);
    Double_t nsum  = 0; // Running weighted sum
    Double_t nsume = 0; // Running weighted sum error
    Double_t dsum  = 0;
    Double_t dsume = 0;
    Int_t    n     = 0;
    Double_t rat   = 1e99;
    
    Int_t k = 0;
    for (k = 0; k < j; k++) {
      Int_t    l      = idx[k]; // Sorted index
      Int_t    ipzBin = hb[l];
      Double_t hvv    = hv[l];      
      Double_t hee    = he[l];
      Double_t hrr    = hr[l];
      Double_t x      = TMath::Sqrt(nsume+hee*hee)/(nsum+hvv);
      if (x > rat) {
	continue; // Ignore - does not help
      }
      Double_t by = h->GetYaxis()->GetBinCenter(ipzBin);
      Int_t    ib = ipz ? ipz->FindBin(by) : 0;
      rat   =  x;
      nsum  += h->GetBinContent(etaBin, ipzBin);
      nsume += TMath::Power(h->GetBinError(etaBin, ipzBin), 2);
      // If we do not have the vertex distribution, then just count
      // number of observations. 
      dsum  += !ipz ? 1 : ipz->GetBinContent(ib);
      dsume += !ipz ? 0 : TMath::Power(ipz->GetBinError(ib),2);
      n++;
    }
    if (k == 0 || n == 0) {
      /// ::Warning("Average", "Eta bin # %3d has no data",etaBin);
      continue; // This eta bin empty!
    }

    Double_t norm = (mode > 0 ? n : dsum);
    Double_t rne  = nsume/nsum/nsum;
    Double_t rde  = dsume/dsum/dsum;
    Double_t av   = nsum/norm;
    Double_t ave  = 0;
    if (mode==2) ave = TMath::Sqrt(nsume)/n;
    else         ave = av*TMath::Sqrt(rne+rde);
#if 0
    Printf("%10s - bin %3d (%+6.3f) count=%2d n=%9.3f+/-%9.3f d=%9.3f+/-%9.3f "
	   "norm=%9.3f -> %9.3f +/- %9.3f",
	   h->GetName(), etaBin, h->GetXaxis()->GetBinCenter(etaBin), n, 
	   nsum, TMath::Sqrt(nsume)/(mode == 2 ? n : nsum),
	   dsum, TMath::Sqrt(dsume)/(mode == 2 ? n : dsum),
	   norm, av, ave);
#endif 
    //Info("AverageSim", "Setting eta bin # %3d to %f +/- %f", etaBin,av,ave);
    p->SetBinContent(etaBin, av);
    p->SetBinError  (etaBin, ave);
  }
  if (mode == 0) p->Scale(1, "width");
  return p;		   
}

//____________________________________________________________________
TH1* Avg2(TH2* h, TH1* ipz, UShort_t mode, const char* name)
{
  TH2* tmp = static_cast<TH2*>(h->Clone("tmp"));
  tmp->SetDirectory(0);
  CorrectForIPz(tmp, ipz, true);
  
  TH1* ret = Avg(tmp, 0, 2, name);
  if (tmp) delete tmp;

  return ret;
}

//____________________________________________________________________
Double_t GetDeltaScale(TSeqCollection* list,
		       Double_t&       eR,
		       Int_t           b,
		       Double_t        nEv=1)
{
  TH1*     deltaRec = GetH1(list,Form("b%d_TrData_WDist", b));
  TH1*     deltaBg  = GetH1(list,Form("b%d_TrInj_WDist", b));
  deltaRec->Scale(1./nEv);
  deltaBg ->Scale(1./nEv);
  Double_t top      = deltaRec->GetXaxis()->GetXmax();
  Double_t eRec, eBg;
  Double_t iRec     = Integrate(deltaRec, 5, top, eRec);
  Double_t iBg      = Integrate(deltaBg,  5, top, eBg);
  //Printf("Integral of reconstructed Delta: %14.1f +/- %14.2f", iRec, eRec);
  //Printf("Integral of injected Delta:      %14.1f +/- %14.2f", iBg,  eBg);
  Double_t r = RatioE(iRec, eRec, iBg, eBg, eR);
  Printf("%20s: data=%9g+/-%9g  inj=%9g+/-%9g  ->  %9g+/-%9g",
	 list->GetName(), iRec, eRec, iBg, eBg, r, eR);
  return r;
}

TH1* GetDeltaDist(TSeqCollection* list, Int_t b)
{
  TH2*     deltaRec = GetH2(list,Form("b%d_TrData_WDvsEta", b));
  TH2*     deltaBg  = GetH2(list,Form("b%d_TrInj_WDvsEta", b));
  Double_t top      = deltaRec->GetYaxis()->GetXmax();
  TH1*     ret      = deltaRec->ProjectionX("deltaInt");
  ret->Reset();
  ret->SetYTitle("\\delta\\hbox{ scale}");
  ret->SetXTitle("\\eta");
  ret->SetDirectory(0);
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++) {
    TH1* tmpRec = deltaRec->ProjectionY("tmpRec",i,i);
    TH1* tmpBg  = deltaBg ->ProjectionY("tmpBg", i,i);
    Double_t eRec, eBg, eR;
    Double_t iRec     = Integrate(tmpRec, 5, top, eRec);
    Double_t iBg      = Integrate(tmpBg,  5, top, eBg);
    Double_t r        =  RatioE(iRec, eRec, iBg, eBg, eR);
    ret->SetBinContent(i, r);
    ret->SetBinError  (i, eR);
    delete tmpRec;
    delete tmpBg;
  }
  return ret;
}

//____________________________________________________________________
TH2* Scale(TH2* h, Double_t x, Double_t xe)
{
  Double_t rr    = xe/x;
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      Double_t c  = h->GetBinContent(i,j);
      Double_t e  = h->GetBinError  (i,j);
      Double_t s  = (c > 0 ? e/c : 0);
      Double_t sc = x * c;
      Double_t se = sc*TMath::Sqrt(s*s+rr*rr);
      // Printf("Setting bin %2d,%2d to %f +/- %f", sc, se);
      h->SetBinContent(i,j,sc);
      h->SetBinError  (i,j,se);
    }
  }
  return h;
}
//____________________________________________________________________
TH2* Scale(TH2* h, TH1* g)
{
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t dc = g->GetBinContent(i);
    Double_t de = g->GetBinError  (i);
    Double_t dr = (dc > 1e-6 ? 1/dc : 0);
    Double_t dq = dr * de;
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      Double_t nc  = h->GetBinContent(i,j);
      Double_t ne  = h->GetBinError  (i,j);
      Double_t ns  = (nc > 0 ? ne/nc : 0);
      Double_t sc  = dr * nc;
      Double_t se  = sc*TMath::Sqrt(ns*ns+dq*dq);
      // Printf("Setting bin %2d,%2d to %f +/- %f", sc, se);
      h->SetBinContent(i,j,sc);
      h->SetBinError  (i,j,se);
    }
  }
  return h;
}


//____________________________________________________________________
TH2* GetBg(Bool_t          comb,
	   TSeqCollection* list,
	   Int_t           b)
{
  Double_t eR = 0;
  TH2* in     = GetH2(list,Form("b%d_Tr%s_ZvEtaCutT",b,comb ? "Comb" : "Inj"));
  TH2* ret    = static_cast<TH2*>(in->Clone("bg"));
  in ->SetDirectory(0);
  ret->SetDirectory(0);
  ret->SetStats(0);
  if (comb) return ret;

#if 1
  // We should do the below, but that doesn't work as it creates
  // negative signals!
  Double_t r     = GetDeltaScale(list, eR, b);
  // Printf("Ratio of integrals:              %14.3f +/- %f14.3", r, eR);
  return Scale(ret, r, eR);
#else
  TH1*    r      = GetDeltaDist(list, b);
  Scale(ret, r);
  // PrintH(r);
  delete r;
  // PrintH(ret);
  return ret;
#endif 
}

/** 
 * Correct measured distribution for background.  If we do not use
 * MC-labels then we simply substract the background.  
 *
 * If we use MC labels, then multiply the measures signal by
 * @f$ 1-\beta@f$, where @f$ \beta@f$ is 
 *
 * @f[
 *   C/M\prime
 * @f]
 *
 * where @f$ C@f$ is the combinatorial background, and @f$ M\prime@f$
 * is the measured distribution in simulations.
 * 
 * @param comb     true if we use combinatorial background 
 * @param meas     The measured distribution to correct (in-place)
 * @param bg       The background (or combinatorial background)
 * @param simMeas  The measured distribution (in MC)
 * 
 * @return Pointer to @f$ meas@f$ after correcting 
 */
TH2* CorrectSignal(Bool_t comb, TH2* meas, TH2* bg, TH2* simMeas)
{
  if (!comb) {
    meas->Add(bg, -1);
    return meas;
  }
  TH2* beta = static_cast<TH2*>(bg->Clone("beta"));
  beta->Divide(simMeas);
  beta->SetDirectory(0);
  TH2* one = static_cast<TH2*>(bg->Clone("one"));
  one->Reset();
  one->SetDirectory(0);
  for (Int_t i = 1; i <= one->GetNbinsX(); i++) 
    for (Int_t j = 1; j <= one->GetNbinsY(); j++)
      one->SetBinContent(i,j,1);
  one->Add(beta,-1);
  meas->Multiply(one);
  delete one;
  delete beta;    
  return meas;
}
/** 
 * Calculates the result as 
 *
 * @f[
 *  R &=& \frac{P}{(1-\beta\prime)M\prime} (1-\beta)M
 * @f] 
 *
 * where 
 * 
 * - @f$ P@f$ is the input primary particle distribution 
 * - @f$ \beta\prime@f$ is background correction for simulated data 
 * - @f$ M\prime@f$ is the observed distribution in simulations 
 * - @f$ \beta@f$ is background correction for real data 
 * - @f$ M@f$ is the observed distribution in real data 
 * 
 * In case we use MC-labels for the background correction in
 * simulations, we have that
 *
 * @f[
 *   (1-\beta\prime)M\prime = (1-C\prime/M\prime)M\prime = M\prime-C\prime
 * @f] 
 * 
 * where 
 *
 * - @f$ C\prime@f$ is the observed number of tracklets for which the 
 *   constituent clusters have different parent simulated traks.
 *
 * In case we use MC-labels for the background correction for real
 * data, we have that
 *
 * @f[ 
 *  (1-\beta)M = (1-k\beta\prime)M = (1-kC/M\prime)M
 * @f] 
 *
 * where @f$ k@f$ is a constant scaling factor (1.3). 
 *
 * If we use injected clusters as an estimator for the background
 * correction in simulations, we have 
 *
 * @f[
 *  (1-\beta\prime)M\prime = (1-B\prime/M\prime)M\prime = M\prime-B\prime
 * @f] 
 *
 * where 
 *
 * - @f$ B\prime@f$ is the observed number of tracklets when injecting
 *   fake clusters, scaled so that the corresponding @f$\Delta@f$
 *   distribution match the measured (in simulations) @f$\Delta@f$
 *   distribution in the tails (@f$\Delta\in[5,20]@f$)
 *
 * The same expression holds for the case when we use injected
 * clusters to estimate the background in real data.
 *
 *  (1-\beta)M = (1-B/M)M = M-B
 *
 * where @f$ B@f$ has the same meaning as @f$ B\prime@f$ except we use
 * real data. 
 *
 * Note, if we use MC-labels for both the simulated and real data
 * backgrounds, and @f$ k=1@f$ we find that 
 *
 * @f[
 *  R = \frac{P}{M\prime}M
 * @f] 
 *
 * @param stack     Stack to add results to 
 * @param dir       Directory to write to 
 * @param b         Bin number 
 * @param realList  List of real data histograms 
 * @param simList   List of simulated data histograms 
 * @param realComb  Whether to use MC-labels for real data 
 * @param simComb   Whether to use MC-labels for simulated data 
 * @param preScale  Whether to scale all histograms to the vertex
 *                  distribution fully before doing any calculations. 
 * @param full      If true, pre-scale fully (including errors)
 */
void CorrectOneBin(THStack*        stack,
		   TDirectory*     dir,
		   Int_t           b,
		   TSeqCollection* realList,
		   TSeqCollection* simList,
		   Bool_t          realComb,
		   Bool_t          simComb, 
		   Bool_t          preScale=false,
		   Bool_t          full=false,
		   Bool_t          showOne=false,
		   UShort_t        scaleMode=2)
{
  const Double_t k = (realList == simList ? 1 : 1.3);
  realList->SetName("realList");
  simList ->SetName("simList");
  TH2* realMeas = GetH2(realList, Form("b%d_TrData_ZvEtaCutT", b));
  TH2* realIPzC = GetH2(realList, "zv");
  TH1* realIPz  = realIPzC->ProjectionX("realIPz", b+1,b+1, "e");
  TH2* simMeas  = Coerce(realMeas,
			 GetH2(simList,  Form("b%d_TrData_ZvEtaCutT", b)));
  TH2* simIPzC  = Coerce(realIPzC, GetH2(simList,  "zv"));
  TH1* simIPz   = simIPzC->ProjectionX("simIPz", b+1,b+1, "e");
  TH2* trueGen  = Coerce(realMeas,
			 GetH2(simList,  Form("b%d_zvEtaPrimMC",      b)));
  TH2* trueIPzC = Coerce(realIPzC, GetH2(simList,  "zvMCNoPS"));
  TH1* trueIPz  = simIPzC->ProjectionX("trueIPz", b+1,b+1, "e");
  TH2* simBg    = Coerce(realMeas, GetBg(simComb,  simList, b));
  TH2* realBg   = Coerce(realMeas, GetBg(realComb, (realComb ?
						    simList :
						    realList), b));
  Double_t scale  = 1;
  Double_t scaleE = 0;
  TH1*     histK  = 0;
  // PrintH(simBg);
  // PrintH(realBg);
  if (realComb) {
    if      (scaleMode == 0) {
      // Fixed scale
      scale = k;
      Scale(realBg, scale, 0);
    }
    else if (scaleMode == 1) {
      // Eta indepedent scale 
      Double_t simRe, realRe;
      Double_t realR  = GetDeltaScale(realList, realRe, b,
				      realIPz->GetEntries());
      Double_t simR   = GetDeltaScale(simList,  simRe,  b,
				      simIPz->GetEntries());
      scale           = RatioE(realR, realRe, simR, simRe, scaleE);
      Scale(realBg, scale, scaleE);
    }
    else {
      // Eta dependent scale 
      TH1* tmpReal = GetDeltaDist(realList, b);
      TH1* tmpSim  = GetDeltaDist(simList,  b);
      tmpReal->Divide(tmpSim);
      tmpReal->SetName("k");
      tmpReal->SetTitle("Background scale");
      tmpReal->SetMarkerColor(tmpReal->GetMarkerColor());
      tmpReal->SetMarkerStyle(tmpReal->GetMarkerStyle()+4);
      tmpReal->SetLineColor(tmpReal->GetLineColor());
      histK = tmpReal;
      Scale(realBg, tmpReal);
    }
  }
  if (!CheckConsistency(realMeas, simMeas)) {
    Warning("CorrectOneBin", "Real data (%s) and simulated data (%s) done "
	    "with different binning", realList->GetName(), simList->GetName());
    return;
  }
  realMeas->SetName("realMeas");
  simMeas ->SetName("simMeas");
  realBg  ->SetName("realBg");
  simBg   ->SetName("simBg");
  trueGen ->SetName("trueGen");
  realMeas->SetTitle("Measured (real)");
  realBg  ->SetTitle("Background (real)");
  simMeas ->SetTitle("Measured (simulated)");
  simBg   ->SetTitle("Background (simulated)");
  trueGen ->SetTitle("Generated");
  realMeas->SetDirectory(0);
  simMeas ->SetDirectory(0);
  trueGen ->SetDirectory(0);
  realIPz ->SetDirectory(0);
  simIPz  ->SetDirectory(0);
  trueIPz ->SetDirectory(0);
  
  // Scale IPz histograms
  Double_t dummy;
  Double_t ipzMin  = realMeas->GetYaxis()->GetXmin();
  Double_t ipzMax  = realMeas->GetYaxis()->GetXmax();
  Double_t realNev = Integrate(realIPz,ipzMin,ipzMax,dummy);
  Double_t simNev  = Integrate(simIPz, ipzMin,ipzMax,dummy);
  Double_t trueNev = Integrate(trueIPz,ipzMin,ipzMax,dummy);
  realIPz->SetMarkerColor(kRed+2);  realIPz->SetMarkerStyle(20);
  simIPz ->SetMarkerColor(kBlue+2); simIPz ->SetMarkerStyle(21);
  trueIPz->SetMarkerColor(kBlack);  trueIPz->SetMarkerStyle(24);
#if 0
  // This is done in Rubens code 
  for (Int_t i = 1; i <= realIPz->GetNbinsX(); i++) {
    if (simIPz->GetBinContent(i) < 1e-6) {
      realIPz->SetBinContent(i,0);
      realIPz->SetBinError  (i,0);
    }
  }
#endif
  
  // Either scale fully here - up front 
  if (preScale) {
    CorrectForIPz(realMeas, realIPz,                       full);
    CorrectForIPz(realBg,   (realComb ? simIPz : realIPz), full);
    CorrectForIPz(simMeas,  simIPz,                        full);
    CorrectForIPz(simBg,    simIPz,                        full);
    CorrectForIPz(trueGen,  trueIPz,                       full);
  }
  // Or first scale by number of events in each sample 
  else {
    // Printf("Number of real events: %f", realNev);
    realMeas->Scale(1/realNev);
    realBg  ->Scale(1/(realComb ? simNev : realNev));
    realIPz ->Scale(1/realNev);
    simMeas ->Scale(1/simNev);
    simBg   ->Scale(1/simNev);
    simIPz  ->Scale(1/simNev);
    trueGen ->Scale(1/trueNev);
    trueIPz ->Scale(1/trueNev);
  }
  
  // Calculate real signal as measured minus background
  TH2* realSig = static_cast<TH2*>(realMeas->Clone("realSig"));
  realSig->SetTitle("Signal (real)");
  realSig->SetDirectory(0);
  CorrectSignal(realComb, realSig, realBg, simMeas);
  
  // Calculate simulated signal as measured minus background
  TH2* simSig = static_cast<TH2*>(simMeas->Clone("simSig"));
  simSig->SetTitle("Signal (simulated)");
  simSig->SetDirectory(0);
  // CorrectSignal(simComb, simSig, simBg, simMeas);
  simSig->Add(simBg, -1);

  // Calculate alpha as primary over simulated signal
  TH2* alpha = static_cast<TH2*>(trueGen->Clone("alpha"));
  alpha->SetTitle("\\alpha");
  alpha->SetDirectory(0);
  alpha->Divide(simSig);
  // Create fiducial cut as cut on alpha 
  TH2* fiducial = static_cast<TH2*>(alpha->Clone("fiducial"));
  fiducial->SetTitle("Fiducial cut");
  fiducial->SetDirectory(0);
  for (Int_t i = 1; i <= fiducial->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= fiducial->GetNbinsY(); j++) {
      Double_t a = fiducial->GetBinContent(i,j);
      fiducial->SetBinError(i,j,0);
      fiducial->SetBinContent(i,j, (a > 0 && a <= 2.5));
    }
  }
  alpha   ->Multiply(fiducial);

  // Create result as alpha times real signal 
  TH2* result = static_cast<TH2*>(realSig->Clone("result"));
  result->SetTitle("Result");
  result->SetDirectory(0);
  result->Multiply(alpha);  

  // The average mode
  UShort_t mode = (preScale ? (full ? 2 : 1) : 0);
  // Calculate the primary dN/deta
  TH1* truth  = Avg(trueGen, trueIPz, mode, "truth");
  truth->SetYTitle(dndetaTitle);
  SetAttr(truth, cc[b%11], 24, 1.5, 0, 3, 2);

  // Calculate the real dN/deta 
  TH1* dndeta = Avg(result, realIPz, mode, "dndeta");
  dndeta->SetYTitle(dndetaTitle);
  SetAttr(dndeta, cc[b%11], 20, 1.2, 0, 1, 1);
  if (full) {
    for (Int_t i = 1; i <= dndeta->GetNbinsX(); i++) {
      dndeta->SetBinError(i,1./3*dndeta->GetBinError(i));
    }
  }
  Double_t lim = 0.5;
  TF1* f = new TF1("ff", "pol0", -lim, lim);
  TFitResultPtr r = dndeta->Fit(f, "QN0RS", "", -lim, lim);
  Printf("dNch/deta %20s: %6.1f +/- %6.1f (%5.2f - %5.3f+/-%5.3f)",
	 dir->GetName(),
	 r->Parameter(0), r->ParError(0), r->Chi2()/r->Ndf(),
	 scale, scaleE);


  // Multiply histograms by fiducial cut
  realMeas->Multiply(fiducial);
  simMeas ->Multiply(fiducial);
  realSig ->Multiply(fiducial);
  simSig  ->Multiply(fiducial); 
  realBg  ->Multiply(fiducial);
  simBg   ->Multiply(fiducial);
  realSig ->Multiply(fiducial);
  simSig  ->Multiply(fiducial);
  TH1* realAvgMeas = Avg(realMeas,realIPz, mode, "realAvgMeas",result);
  TH1* simAvgMeas  = Avg(simMeas, simIPz,  mode, "simAvgMeas", result);
  TH1* realAvgBg   = Avg(realBg,  realIPz, mode, "realAvgBg",  result);
  TH1* simAvgBg    = Avg(simBg,   simIPz,  mode, "simAvgBg",   result);
  TH1* realAvgSig  = Avg(realSig, realIPz, mode, "realAvgSig", result);
  TH1* simAvgSig   = Avg(simSig,  simIPz,  mode, "simAvgSig",  result);
  SetAttr(realAvgMeas, kRed+2, 21);
  SetAttr(realAvgMeas, kRed+2, 22);
  SetAttr(realAvgBg,   kGreen+2, 21);
  SetAttr(realAvgSig,  kGreen+2, 22);
  SetAttr(simAvgBg,    kBlue+2, 21);
  SetAttr(simAvgSig,   kBlue+2, 22);
  SetAttr(histK,       kGreen+2, 23);
  
  THStack* summary = new THStack("summary", dndeta->GetYaxis()->GetTitle());
  summary->Add(dndeta, "e2");
  summary->Add(truth,  "E");
  summary->Add(realAvgMeas);
  summary->Add(simAvgMeas);
  summary->Add(realAvgBg);
  summary->Add(simAvgBg);
  summary->Add(realAvgSig);
  summary->Add(simAvgSig);
  if (histK) {
    // histK->Scale(10);
    // histK->SetTitle(Form("%s\\times10",histK->GetTitle()));
    summary->Add(histK);
  }

  if (preScale) {
    // Scale vertex distributions
    realIPz->Scale(1/realNev);
    simIPz ->Scale(1/simNev);
    trueIPz->Scale(1/trueNev);
  }
  if (showOne){
    TCanvas* c = new TCanvas(Form("c%02d",b),
			     Form("Centrality bin %d",b),
			     1200, 1000);
    c->Divide(4,3); //,0,0);
    c->cd(1);  realMeas->DrawCopy("lego2 e");
    c->cd(2);  realBg  ->DrawCopy("lego2 e");
    c->cd(3);  realSig ->DrawCopy("lego2 e");
    c->cd(4);  realIPz ->DrawCopy();
    c->cd(4);  simIPz  ->DrawCopy("same");
    c->cd(4);  trueIPz ->DrawCopy("same");  
    c->cd(5);  simMeas ->DrawCopy("lego2 e");
    c->cd(6);  simBg   ->DrawCopy("lego2 e");
    c->cd(7);  simSig  ->DrawCopy("lego2 e");
    c->cd(8);  alpha   ->DrawCopy("lego2 e");
    c->cd(9);  fiducial->DrawCopy("lego2 e");
    c->cd(10); trueGen ->DrawCopy("lego2 e");
    c->cd(11); result  ->DrawCopy("lego2 e");
    c->cd(12); summary ->DrawClone("nostack");
    c->GetPad(4)->Modified();
    c->GetPad(12)->Modified();
    c->GetPad(4)->BuildLegend();
    c->GetPad(12)->BuildLegend();
    TLatex* ltx = new TLatex(.6, .3,
			     Form("%s|_{|\\eta|<%3.1f}"
				  "=%6.1f\\pm%6.1f (\\chi^2/\\nu=%5.2f)",
				  "dNch/deta", lim, r->Parameter(0),
				  r->ParError(0), r->Chi2()/r->Ndf()));
    ltx->SetTextAlign(23);
    ltx->SetNDC(true);
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.035);
    c->cd(12); ltx->Draw();
  }
  stack->Add(truth);
  stack->Add(dndeta, "e2");

  if (!dir) return;

  dir->cd();
  result ->Write();
  summary->Write();
  dndeta ->Write();
  if (histK) histK      ->Write();
  (new TParameter<double>("scale",  scale)) ->Write();
  (new TParameter<double>("scaleE", scaleE))->Write();

  TDirectory* det = dir->mkdir("details");
  det->cd();
  if (!preScale) {
    CorrectForIPz(realMeas, realIPz,                       full);
    CorrectForIPz(realBg,   (realComb ? simIPz : realIPz), full);
    CorrectForIPz(simMeas,  simIPz,                        full);
    CorrectForIPz(simBg,    simIPz,                        full);
    CorrectForIPz(trueGen,  trueIPz,                       full);
  }  
  realMeas ->Write();
  realBg   ->Write();
  realSig  ->Write();
  simMeas  ->Write();
  simBg    ->Write();
  simSig   ->Write();
  trueGen  ->Write();
  alpha    ->Write();
  fiducial ->Write();
  GetH1(realList, Form("b%d_TrData_WDist", b))->Write("realDataDelta");
  GetH1(realList, Form("b%d_TrInj_WDist",  b))->Write("realInjDelta");
  GetH1(simList,  Form("b%d_TrData_WDist", b))->Write("simDataDelta");
  GetH1(simList,  Form("b%d_TrInj_WDist",  b))->Write("simInjDelta");
  
  TDirectory* avs = dir->mkdir("averages");
  avs->cd();
  realAvgMeas->Write();
  realAvgSig ->Write();
  realAvgBg  ->Write();
  simAvgMeas ->Write();
  simAvgSig  ->Write();
  simAvgBg   ->Write();
  realIPz    ->Write();
  simIPz     ->Write();
  trueIPz    ->Write();
  truth      ->Write();
  
}

//____________________________________________________________________
/** 
 * Print the stats stored in hStat histogram 
 * 
 * @param list List to get histogram from 
 */
void
ShowStats(TList* list, const char* title, const char* path)
{
  Printf("=== %-15s - %s ===", title, path);
  TH1*   tmp   = GetH1(list, "hStat");
  if (!tmp) {
    Warning("", "No stats");
    return;
  }
  TH1*   stats = static_cast<TH1*>(tmp->Clone());
  stats->SetDirectory(0);
  TAxis* axis  = stats->GetXaxis();
  for (Int_t i = 1; i <= 4; i++) {
    Printf("%31s: %f", axis->GetBinLabel(i), stats->GetBinContent(i));
  }
  Double_t scale = stats->GetBinContent(4);
  stats->Scale(1/scale);
  for (Int_t i = 6; i <= 23; i++) {
    Printf("%31s: %f", axis->GetBinLabel(i), stats->GetBinContent(i));
  }
}  
  
/** 
 * Correct measured tracklet distributions using MC.  The flags select
 * how to do the correction.  It is a bit mask of options.
 *
 * - 0x1:  Use combinatorial background for real data 
 * - 0x2:  Use combinatorial background for simulated data 
 * - 0x4:  Pre-scale histograms to vertex bins 
 * - 0x8:  When pre-scaling, also propagate errors from vertex dist. 
 * - 0x10: Instead of real data, use the simulated data (closure test). 
 * 
 * The input is two files generated by running AliTrackletTaskMulti on
 * real and simulated data.  The files are expected to be in trdt.root
 * for real data and trmc.root for simulated data.
 * 
 * @param flags Flags for the processing 
 * @param nBins Number of bins to process 
 */
void SimpleCorrect(UShort_t flags,
		   const char* realFileName,
		   const char* simFileName,
		   Int_t       nBins=9,
		   const char* otherFileName="")
{
  TFile* realFile  = TFile::Open(realFileName,  "READ");
  TFile* simFile   = TFile::Open(simFileName,   "READ");
  TFile* otherFile = 0;
  if (otherFileName && otherFileName[0] != '\0') 
    otherFile = TFile::Open(otherFileName, "READ");
  TList* realList  = static_cast<TList*>(realFile->Get("clist"));
  TList* simList   = static_cast<TList*>(simFile ->Get("clist"));
  TObjArray* otherList = 0;
  if (otherFile)
    otherList =static_cast<TObjArray*>(otherFile->Get("TObjArray"));

  ShowStats(realList, "Real data",     realFile->GetName());
  ShowStats(simList,  "Simulated data",simFile ->GetName());

  TH1* realCent = GetH1(realList, "EvCentr");
  TH1* simCent  = GetH1(simList,  "EvCentr");
  if (!CheckConsistency(realCent, simCent)) {
    Warning("SimpleCorrection","Inconsistent centralities");
    return;
  }
  TH2* realIpz = GetH2(realList, "b0_TrData_ZvEtaCutT");
  TH2* simIpz  = GetH2(simList,  "b0_TrData_ZvEtaCutT");
  if (!CheckConsistency(realIpz, simIpz)) {
    Warning("SimpleCorrection", "Inconsistent IPz or eta axis");
    return;
  }
  
  
  Bool_t   realComb  = flags & 0x1;
  Bool_t   simComb   = flags & 0x2;
  Bool_t   preScale  = flags & 0x4;
  Bool_t   full      = flags & 0x8;
  Bool_t   showOne   = flags & 0x20;
  UShort_t scaleMode = (flags & 0x40 ? 0 : (flags & 0x80 ? 1 : 2));
  Printf("Real background:      %s\n"
	 "Simulated background: %s\n"
	 "Pre-scale:            %s\n"
	 "Full pre-scale:       %s\n"
	 "Show each bin:        %s\n"
	 "Scaling mode:         %s",
	 (realComb ? "comb" : "inj"),
	 (simComb  ? "comb" : "inj"),
	 (preScale ? "yes"  : "no"),
	 (full     ? "yes"  : "no"),
	 (showOne  ? "yes"  : "no"),
	 (scaleMode==0 ?"fixed":scaleMode==1?"c":"eta-c"));
	 
	 
  if (flags & 0x10) realList = simList;
  TString tit;
  tit.Form("Real BG: %s  Simulated BG: %s%s",
	   realComb ? "MC-labels" : "Injection",
	   simComb  ? "MC-labels" : "Injection",
	   preScale ? " (pre-scaled)" : "");

  TFile*   output = TFile::Open(Form("result_0x%x.root", flags&0x3),
				"RECREATE");
  realCent->Write("cent");
  THStack* res    = new THStack("result", tit);
  TH1*     cent   = GetH1(realList, "EvCentr");
  for (Int_t b = 0; b < nBins; b++) {
    Double_t c1 = cent->GetXaxis()->GetBinLowEdge(b+1);
    Double_t c2 = cent->GetXaxis()->GetBinUpEdge (b+1);
    TString  name;
    name.Form("cent%03dd%02d_%03dd%02d",
	      Int_t(c1), Int_t(c1*100)%100,
	      Int_t(c2), Int_t(c2*100)%100);
    TDirectory* dir = output->mkdir(name);
    CorrectOneBin(res,dir,b,realList,simList,realComb,simComb,preScale,full,
		  showOne, scaleMode);

    if (otherList) {
      TH1* od = GetH1(otherList, Form("bin%d_DataCorrSignal_Bla",b));
      TH1* ot = GetH1(otherList, Form("bin%d_MCTruth",b));
      TH1* td = static_cast<TH1*>(od->Clone("otherdNdeta"));
      TH1* tt = static_cast<TH1*>(ot->Clone("otherTruth"));
      td->SetDirectory(0);
      tt->SetDirectory(0);
      SetAttr(td, cc[b%11], 21, 1.2, 3004, 1, 1);
      SetAttr(tt, cc[b%11], 25, 1.4, 0, 7, 3);
      if (td->GetMinimum() < 1e-6) { td->SetMinimum(); td->SetMaximum(); }
      if (tt->GetMinimum() < 1e-6) { tt->SetMinimum(); tt->SetMaximum(); }
      res->Add(td, "e2");
      res->Add(tt);
      dir->cd();
      td->Write();
      tt->Write();
    }
  }
  // TCanvas* all = new TCanvas("all", "all");
  // res->Draw("nostack");
  res->SetMaximum(1.2*res->GetMaximum("nostack"));
  output->cd();
  res->Write();
  output->Write();

}
//
// EOF
//
