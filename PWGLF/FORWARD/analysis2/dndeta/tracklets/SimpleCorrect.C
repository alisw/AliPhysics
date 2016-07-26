class TCanvas; // Auto load 

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

//____________________________________________________________________
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

//____________________________________________________________________
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
	 const char* name) 
{
  if (!h) return 0;
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
    TArrayD iv(nIPz);
    TArrayD ie(nIPz);
    Int_t   j = 0;
    for (Int_t ipzBin = 1; ipzBin <= nIPz; ipzBin++) {
      Double_t bc = h->GetBinContent(etaBin, ipzBin);
      if (bc < 1e-9) continue; // Low value
      Double_t be = h->GetBinError  (etaBin, ipzBin);
      Double_t by = h->GetYaxis()->GetBinCenter(ipzBin);
      Int_t    ib = ipz ? ipz->FindBin(by) : 0;
      hv[j] = bc;
      he[j] = be;
      hr[j] = be/bc;
      // If we do not have the vertex distribution, then just count
      // number of observations. 
      iv[j] = !ipz ? 1 : ipz->GetBinContent(ib);
      ie[j] = !ipz ? 0 : ipz->GetBinError  (ib);
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
      Int_t    l   = idx[k]; // Sorted index
      Double_t hvv = hv[l];      
      Double_t hee = he[l];
      Double_t ivv = iv[l];
      Double_t iee = ie[l];
      Double_t hrr = hr[l];
      Double_t x = TMath::Sqrt(nsume+hee*hee)/(nsum+hvv);
#if 0
      Printf(" %10s - bin %3d,%3d %9.3f+/-%9.3f (%5.2f%%) [%9.3+/-%9.3f] "
	     "-> %9.3f (%9.3g)",
	     h->GetName(), etaBin, l, hvv, hee, 100*hrr, ivv, iee, x, rat);
#endif 
      if (x > rat) {
	continue; // Ignore - does not help
      }
      
      rat   =  x;
      nsum  += hvv;
      nsume += hee*hee;
      dsum  += ivv;
      dsume += iee*iee;
      n++;
    }
    if (k == 0 || n == 0) {
      ::Warning("Average", "Eta bin # %3d has no data",etaBin);
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

//____________________________________________________________________
void PrintH(TH2* h, Int_t prec=2)
{
  Printf("Content of %s - %s", h->GetName(), h->GetTitle());
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    printf("%3d: ", i);
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      printf("%.*f+/-%.*f ",
	     prec, h->GetBinContent(i,j),
	     prec, h->GetBinError(i,j));
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

//____________________________________________________________________
TH1* SetAttr(TH1* h,
	     Color_t  color,
	     Style_t  marker=20,
	     Double_t size=1.,
	     Style_t  fill=0,
	     Style_t  line=1,
	     Width_t  width=1)
{
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
TH2* GetBg(Bool_t comb, TSeqCollection* list, Int_t b)
{
  if (comb) {
    TH2* in  = GetH2(list,  Form("b%d_TrComb_ZvEtaCutT", b));
    TH2* ret = static_cast<TH2*>(in->Clone("bg"));
    ret->SetDirectory(0);
    ret->SetStats(0);
    return ret;
  }

  TH2* in  = GetH2(list,  Form("b%d_TrInj_ZvEtaCutT", b));
  TH2* ret = static_cast<TH2*>(in->Clone("bg"));
  in ->SetDirectory(0);
  ret->SetDirectory(0);
  ret->SetStats(0);
  TH1*  deltaRec = GetH1(list, Form("b%d_TrData_WDist", b));
  TH1*  deltaBg  = GetH1(list, Form("b%d_TrInj_WDist", b));
  Double_t top = deltaRec->GetXaxis()->GetXmax();
  Double_t eRec, eBg, eR;
  Double_t iRec  = Integrate(deltaRec, 5, top, eRec);
  Double_t iBg   = Integrate(deltaBg,  5, top, eBg);
  //Printf("Integral of reconstructed Delta: %14.1f +/- %14.2f", iRec, eRec);
  //Printf("Integral of injected Delta:      %14.1f +/- %14.2f", iBg,  eBg);
  Double_t r     = RatioE(iRec, eRec, iBg, eBg, eR);
  //Printf("Ratio of integrals:              %14.3f +/- %f14.3", r, eR);
  Double_t rr    = eR/r;
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= ret->GetNbinsY(); j++) {
      Double_t c  = ret->GetBinContent(i,j);
      Double_t e  = ret->GetBinError  (i,j);
      Double_t s  = (c > 0 ? e/c : 0);
      Double_t sc = r * c;
      Double_t se = sc*TMath::Sqrt(s*s+rr*rr);
      // Printf("Setting bin %2d,%2d to %f +/- %f", sc, se);
      ret->SetBinContent(i,j,sc);
      ret->SetBinError  (i,j,se);
    }
  }
  // PrintH(ret);
  return ret;
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
		   Bool_t          full=false)
{
  const Double_t k = (realList == simList ? 1 : 1.3);
  TH2* realMeas = GetH2(realList, Form("b%d_TrData_ZvEtaCutT", b));
  TH2* realIPzC = GetH2(realList, "zv");
  TH1* realIPz  = realIPzC->ProjectionX("realIPz", b+1,b+1, "e");
  TH2* simMeas  = GetH2(simList,  Form("b%d_TrData_ZvEtaCutT", b));
  TH2* simIPzC  = GetH2(simList,  "zv");
  TH1* simIPz   = simIPzC->ProjectionX("simIPz", b+1,b+1, "e");
  TH2* trueGen  = GetH2(simList,  Form("b%d_zvEtaPrimMC",      b));
  TH2* trueIPzC = GetH2(simList,  "zvMCNoPS");
  TH1* trueIPz  = simIPzC->ProjectionX("trueIPz", b+1,b+1, "e");
  TH2* simBg    = GetBg(simComb, simList, b);
  TH2* realBg   = 0;
  if (realComb) {
    // If we use MC combinatorics as real background, we simply copy
    // that and scale by constant.
    realBg = GetBg(simComb, simList, b);
    realBg->SetDirectory(0);
    realBg->Scale(k);
  }
  else
    realBg = GetBg(false, realList, b);
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
  else if {
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
  TH2* realSig = static_cast<TH2*>(realMeas->Clone("realSignal"));
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

  // Multiply histograms by fiducial cut
  realBg ->Multiply(fiducial);
  simBg  ->Multiply(fiducial);
  realSig->Multiply(fiducial);
  simSig ->Multiply(fiducial);
  alpha  ->Multiply(fiducial);

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
  
  TH1* realAvgBg  = Avg(realBg,  realIPz, mode, "realAvgBg");
  TH1* simAvgBg   = Avg(simBg,   simIPz,  mode, "simAvgBg");
  TH1* realAvgSig = Avg(realSig, realIPz, mode, "realAvgSig");
  TH1* simAvgSig  = Avg(simSig,  simIPz,  mode, "simAvgSig");
  SetAttr(realAvgBg,  kGreen+2, 21);
  SetAttr(realAvgSig, kGreen+2, 22);
  SetAttr(simAvgBg,   kBlue+2, 21);
  SetAttr(simAvgSig,  kBlue+2, 22);
  
  THStack* summary = new THStack("summary", dndeta->GetYaxis()->GetTitle());
  summary->Add(dndeta, "e2");
  summary->Add(truth,  "E");
  summary->Add(realAvgBg);
  summary->Add(simAvgBg);
  summary->Add(realAvgSig);
  summary->Add(simAvgSig);

  if (preScale) {
    // Scale vertex distributions
    realIPz->Scale(1/realNev);
    simIPz ->Scale(1/simNev);
    trueIPz->Scale(1/trueNev);
  }
#if 0
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
  c->GetPad(4)->BuildLegend();
  c->GetPad(12)->BuildLegend();
#endif 
  stack->Add(truth);
  stack->Add(dndeta, "e2");

  if (!dir) return;

  dir->cd();
  realMeas ->Write();
  realIPz  ->Write();
  realBg   ->Write();
  realSig  ->Write();
  simMeas  ->Write();
  simIPz   ->Write();
  simBg    ->Write();
  simSig   ->Write();
  trueGen  ->Write();
  trueIPz  ->Write();
  alpha    ->Write();
  fiducial ->Write();
  dndeta   ->Write();
  truth    ->Write();
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
  TFile* otherFile = TFile::Open(otherFileName, "READ");
  TList* realList  = static_cast<TList*>(realFile->Get("clist"));
  TList* simList   = static_cast<TList*>(simFile ->Get("clist"));
  TObjArray* otherList = 0;
  if (otherFile)
    otherList =static_cast<TObjArray*>(otherFile->Get("TObjArray"));

  Bool_t realComb = flags & 0x1;
  Bool_t simComb  = flags & 0x2;
  Bool_t preScale = flags & 0x4;
  Bool_t full     = flags & 0x8;
  if (flags & 0x10) realList = simList;
  TString tit;
  tit.Form("Real BG: %s  Simulated BG: %s%s",
	   realComb ? "MC-labels" : "Injection",
	   simComb  ? "MC-labels" : "Injection",
	   preScale ? " (pre-scaled)" : "");

  TFile*   output = TFile::Open("result.root", "RECREATE");
  THStack* res    = new THStack("result", tit);
  TH1*     cent   = GetH1(realList, "EvCentr");
  for (Int_t b = 0; b < nBins; b++) {
    Double_t c1 = cent->GetXaxis()->GetBinLowEdge(b+1);
    Double_t c2 = cent->GetXaxis()->GetBinUpEdge (b+1);
    TString  name;
    name.Form("cent%03d%02d_%03d%02d",
	      Int_t(c1), Int_t(c1*100)%100,
	      Int_t(c2), Int_t(c2*100)%100);
    TDirectory* dir = output->mkdir(name);
    CorrectOneBin(res,dir,b,realList,simList,realComb,simComb,preScale,full);

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
  TCanvas* all = new TCanvas("all", "all");
  res->Draw("nostack");
  res->SetMaximum(1.2*res->GetMaximum("nostack"));
  output->cd();
  res->Write();
  output->Write();

  new TBrowser;
}
//
// EOF
//
