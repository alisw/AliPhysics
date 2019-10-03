/*
  several function used for PbPb combined spectra
  Blast Wave is also implemented here
  further documentation will come
  
  author: Roberto Preghenella
  email : preghenella@bo.infn.it
*/


/*****************************************************************/
/* BOLTZMANN
/*****************************************************************/

Double_t
Boltzmann_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * mt * TMath::Exp(-mt / T);
}

TF1 *
Boltzmann(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
  
  TF1 *fBoltzmann = new TF1(name, Boltzmann_Func, 0., 10., 3);
  fBoltzmann->SetParameters(mass, T, norm);
  fBoltzmann->SetParNames("mass", "T", "norm");
  fBoltzmann->FixParameter(0, mass);
  return fBoltzmann;
}

/*****************************************************************/
/* LEVY-TSALLIS */
/*****************************************************************/

Double_t
LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 *
LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  return fLevyTsallis;
}

/*****************************************************************/
/* BOLTZMANN-GIBBS BLAST-WAVE */
/*****************************************************************/

static TF1 *fBGBlastWave_Integrand = NULL;
Double_t
BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{
  
  /* 
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
  */
  
  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];

  Double_t beta = beta_max * TMath::Power(r, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;
  Double_t rho = TMath::ATanH(beta);
  Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
  if (argI0 > 700.) argI0 = 700.;
  Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
  //  if (argI0 > 100 || argI0 < -100)
  //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
  return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
  
}

Double_t
BGBlastWave_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];
  
  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1.);
  return norm * pt * integral;
}

Double_t
BGBlastWave_Func_OneOverPt(const Double_t *x, const Double_t *p)
{
  /* 1/pt dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];
  
  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1.);

  return norm * integral;
}


TF1 *
BGBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 10.);
  return fBGBlastWave;
}

TF1 * BGBlastWave_OneOverPT(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func_OneOverPt, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 10.);
  return fBGBlastWave;
}

/*****************************************************************/
/* TSALLIS BLAST-WAVE */
/*****************************************************************/

static TF1 *fTsallisBlastWave_Integrand_r = NULL;
Double_t
TsallisBlastWave_Integrand_r(const Double_t *x, const Double_t *p)
{
  /* 
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
     p[5] -> q
     p[6] -> y (rapidity)
     p[7] -> phi (azimuthal angle)
  */
  
  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];
  Double_t q = p[5];
  Double_t y = p[6];
  Double_t phi = p[7];

  if (q <= 1.) return r;

  Double_t beta = beta_max * TMath::Power(r, n);
  Double_t rho = TMath::ATanH(beta);
  
  Double_t part1 = mt * TMath::CosH(y) * TMath::CosH(rho);
  Double_t part2 = pt * TMath::SinH(rho) * TMath::Cos(phi);
  Double_t part3 = part1 - part2;
  Double_t part4 = 1 + (q - 1.) * temp_1 * part3;
  Double_t expo = -1. / (q - 1.);
  //  printf("part1=%f, part2=%f, part3=%f, part4=%f, expo=%f\n", part1, part2, part3, part4, expo);
  Double_t part5 = TMath::Power(part4, expo);

  return r * part5;
}

static TF1 *fTsallisBlastWave_Integrand_phi = NULL;
Double_t
TsallisBlastWave_Integrand_phi(const Double_t *x, const Double_t *p)
{
  /* 
     x[0] -> phi (azimuthal angle)
  */
  
  Double_t phi = x[0];
  fTsallisBlastWave_Integrand_r->SetParameter(7, phi);
  Double_t integral = fTsallisBlastWave_Integrand_r->Integral(0., 1.);
  return integral;
}

static TF1 *fTsallisBlastWave_Integrand_y = NULL;
Double_t
TsallisBlastWave_Integrand_y(const Double_t *x, const Double_t *p)
{
  /* 
     x[0] -> y (rapidity)
  */

  Double_t y = x[0];
  fTsallisBlastWave_Integrand_r->SetParameter(6, y);
  Double_t integral = fTsallisBlastWave_Integrand_phi->Integral(-TMath::Pi(), TMath::Pi());
  return TMath::CosH(y) * integral;
}

Double_t
TsallisBlastWave_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t q = p[4];
  Double_t norm = p[5];

  if (!fTsallisBlastWave_Integrand_r)
    fTsallisBlastWave_Integrand_r = new TF1("fTsallisBlastWave_Integrand_r", TsallisBlastWave_Integrand_r, 0., 1., 8);
  if (!fTsallisBlastWave_Integrand_phi)
    fTsallisBlastWave_Integrand_phi = new TF1("fTsallisBlastWave_Integrand_phi", TsallisBlastWave_Integrand_phi, -TMath::Pi(), TMath::Pi(), 0);
  if (!fTsallisBlastWave_Integrand_y)
    fTsallisBlastWave_Integrand_y = new TF1("fTsallisBlastWave_Integrand_y", TsallisBlastWave_Integrand_y, -0.5, 0.5, 0);

  fTsallisBlastWave_Integrand_r->SetParameters(mt, pt, beta_max, temp, n, q, 0., 0.);
  Double_t integral = fTsallisBlastWave_Integrand_y->Integral(-0.5, 0.5);
  return norm * pt * integral;
}

TF1 *
TsallisBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t q = 2., Double_t norm = 1.e6)
{
  
  TF1 *fTsallisBlastWave = new TF1(name, TsallisBlastWave_Func, 0., 10., 6);
  fTsallisBlastWave->SetParameters(mass, beta_max, temp, n, q, norm);
  fTsallisBlastWave->SetParNames("mass", "beta_max", "T", "n", "q", "norm");
  fTsallisBlastWave->FixParameter(0, mass);
  fTsallisBlastWave->SetParLimits(1, 0.01, 0.99);
  fTsallisBlastWave->SetParLimits(2, 0.01, 1.);
  fTsallisBlastWave->SetParLimits(3, 0.1, 10.);
  fTsallisBlastWave->SetParLimits(4, 1., 10.);
  return fTsallisBlastWave;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/


TF1 *
BGBlastWave_SingleFit(TH1 *h, Double_t mass, Option_t *opt = "")
{

  TF1 *f = BGBlastWave(Form("fBGBW_%s", h->GetName()), mass);
  h->Fit(f);
  h->Fit(f);
  h->Fit(f, opt);
  return f;
  
}

Int_t nBW;
TF1 *fBGBW[1000];
TGraphErrors *gBW[1000];

TObjArray *
BGBlastWave_GlobalFit(TObjArray *data, Double_t *mass, Double_t profile = .9, Bool_t fixProfile = kFALSE)
{

  /* get data */
  nBW = data->GetEntries();
  for (Int_t idata = 0; idata < nBW; idata++) {
    gBW[idata] = (TGraphErrors *)data->At(idata);
    gBW[idata]->SetName(Form("gBW%d", idata));
  }

  /* init BG blast-wave functions */
  for (Int_t idata = 0; idata < nBW; idata++) {
    printf("init BG-BlastWave function #%d: mass = %f\n", idata, mass[idata]);
    fBGBW[idata] = BGBlastWave(Form("fBGBW%d", idata), mass[idata]);
  }

  /* display data */
  TCanvas *cBW = new TCanvas("cBW");
  cBW->Divide(nBW, 1);
  for (Int_t idata = 0; idata < nBW; idata++) {
    cBW->cd(idata + 1);
    gBW[idata]->Draw("ap*");
  }
  cBW->Update();

  /* init minuit: nBW normalizations + 3 (beta, T, n) BG-BlastWave params */
  const Int_t nbwpars = 3;
  const Int_t nfitpars = nBW + nbwpars;
  TMinuit *minuit = new TMinuit(nfitpars);
  minuit->SetFCN(BGBlastWave_FCN);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  for (Int_t idata = 0; idata < nBW; idata++)
    minuit->mnparm(idata, Form("norm%d", idata), 1.e6, 1., 0., 0., ierflg);
  // minuit->mnparm(nBW + 0, "<beta>", 0.65, 0.01, 0., 1., ierflg);
  // minuit->mnparm(nBW + 1, "T", 0.1, 0.01, 0., 1., ierflg);
  minuit->mnparm(nBW + 0, "<beta>", 0.55, 0.01, 0., 1., ierflg);
  minuit->mnparm(nBW + 1, "T", 0.13, 0.01, 0., 1., ierflg);
  minuit->mnparm(nBW + 2, "n", profile, 0.1, 0., 10., ierflg);
  if (fixProfile) minuit->FixParameter(nBW + 2);

  /* set strategy */
  arglist[0] = 1;
  minuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);

  /* start MIGRAD minimization */
  arglist[0] = 500000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  /* set strategy */
  arglist[0] = 2;
  minuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);

  /* start MIGRAD minimization */
  arglist[0] = 500000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  /* start IMPROVE minimization */
  arglist[0] = 500000;
  minuit->mnexcm("IMPROVE", arglist, 1, ierflg);

  /* start MINOS */
  arglist[0] = 500000;
  arglist[1] = nBW + 1;
  arglist[2] = nBW + 2;
  arglist[3] = nBW + 3;
  minuit->mnexcm("MINOS", arglist, 4, ierflg);

  /* print results */
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  minuit->mnprin(4, amin);

  /* get parameters */
  Double_t beta, betae, betaeplus, betaeminus, betagcc, temp, tempe, tempeplus, tempeminus, tempgcc, prof, profe, profeplus, profeminus, profgcc;
  minuit->GetParameter(nBW + 0, beta, betae);
  minuit->mnerrs(nBW + 0, betaeplus, betaeminus, betae, betagcc);
  minuit->GetParameter(nBW + 1, temp, tempe);
  minuit->mnerrs(nBW + 1, tempeplus, tempeminus, tempe, tempgcc);
  minuit->GetParameter(nBW + 2, prof, profe);
  minuit->mnerrs(nBW + 2, profeplus, profeminus, profe, profgcc);
  Double_t beta_max = 0.5 * (2. + prof) * beta;
  Double_t norm[1000], norme[1000];
  for (Int_t idata = 0; idata < nBW; idata++)
    minuit->GetParameter(idata, norm[idata], norme[idata]);

  /* printout */
  printf("*********************************\n");
  printf("beta_max = %f\n", beta_max);
  printf("<beta>   = %f +- %f (e+ = %f, e- = %f)\n", beta, betae, betaeplus, betaeminus);
  printf("T        = %f +- %f (e+ = %f, e- = %f)\n", temp, tempe, tempeplus, tempeminus);
  printf("n        = %f +- %f (e+ = %f, e- = %f)\n", prof, profe, profeplus, profeminus);

  /* 1-sigma contour */
  minuit->SetErrorDef(1);
  TGraph *gCont1 = NULL;
  //  gCont1 = (TGraph *) minuit->Contour(50, nBW + 0, nBW + 1);
  if (gCont1) gCont1->SetName("gCont1");

  /* 2-sigma contour */
  minuit->SetErrorDef(4);
  TGraph *gCont2 = NULL;
  //  gCont2 = (TGraph *) minuit->Contour(50, nBW + 0, nBW + 1);
  if (gCont2) gCont2->SetName("gCont2");

  /* display fits */
  for (Int_t idata = 0; idata < nBW; idata++) {
    cBW->cd(idata + 1);
    fBGBW[idata]->SetParameter(4, norm[idata]);
    fBGBW[idata]->SetParameter(1, beta_max);
    fBGBW[idata]->SetParameter(2, temp);
    fBGBW[idata]->SetParameter(3, prof);
    fBGBW[idata]->Draw("same");
  }
  cBW->Update();

  /* histo params */
  TH1D *hBW = new TH1D("hBW", "", 3, 0., 3.);
  hBW->SetBinContent(1, beta);
  hBW->SetBinError(1, betae);
  hBW->SetBinContent(2, temp);
  hBW->SetBinError(2, tempe);
  hBW->SetBinContent(3, prof);
  hBW->SetBinError(3, profe);

  /* BW graph */
  TGraphAsymmErrors *gBetaT = new TGraphAsymmErrors();
  gBetaT->SetName("gBetaT");
  gBetaT->SetPoint(0, beta, temp);
  gBetaT->SetPointEXlow(0, TMath::Abs(betaeminus));
  gBetaT->SetPointEXhigh(0, TMath::Abs(betaeplus));
  gBetaT->SetPointEYlow(0, TMath::Abs(tempeminus));
  gBetaT->SetPointEYhigh(0, TMath::Abs(tempeplus));

  /* prepare output array */
  TObjArray *outoa = new TObjArray();
  for (Int_t idata = 0; idata < nBW; idata++) {
    outoa->Add(gBW[idata]);
    outoa->Add(fBGBW[idata]);
  }
  outoa->Add(cBW);
  outoa->Add(hBW);
  outoa->Add(gBetaT);
  if (gCont1) outoa->Add(gCont1);
  if (gCont2) outoa->Add(gCont2);

  return outoa;

}

void 
BGBlastWave_FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  /* beta -> beta_max */
  Double_t beta = par[nBW+0];
  Double_t T = par[nBW+1];
  Double_t n = par[nBW+2];
  Double_t beta_max = 0.5 * (2. + n) * beta;
#if 0
  /* check beta_max */
  if (beta_max >= 1. || beta_max <= 0.) {
    f = kMaxInt;
    return;
  }
  /* check T */
  if (T <= 0.) {
    f = kMaxInt;
    return;
  }
#endif

  Double_t pt, pte, val, vale, func, pull, chi = 0;
  /* loop over all the data */
  for (Int_t iBW = 0; iBW < nBW; iBW++) {
    /* set BGBW parameters */
    fBGBW[iBW]->SetParameter(4, par[iBW]);
    fBGBW[iBW]->SetParameter(1, beta_max);
    fBGBW[iBW]->SetParameter(2, T);
    fBGBW[iBW]->SetParameter(3, n);
    /* loop over all the points */
    for (Int_t ipt = 0; ipt < gBW[iBW]->GetN(); ipt++) {
      pt = gBW[iBW]->GetX()[ipt];
      pte = gBW[iBW]->GetEX()[ipt];
      val = gBW[iBW]->GetY()[ipt];
      vale = gBW[iBW]->GetEY()[ipt];
      func = fBGBW[iBW]->Eval(pt);
      //      func = fBGBW[iBW]->Integral(pt - pte, pt + pte);
      pull = (val - func) / vale;
      chi += pull * pull;
    }
  }

  f = chi;
}

/*****************************************************************/

IntegratedProduction(TH1 *h, Double_t mass, Option_t *opt = "")
{

  Double_t yield, yielderr, yielderrcorr, mean, meanerr, meanerrcorr, partyield[3], partyielderr[3], partyielderrcorr[3];
  TF1 *f = BGBlastWave_SingleFit(h, mass, opt);
  GetYieldAndMean(h, f, yield, yielderr, yielderrcorr, mean, meanerr, meanerrcorr, 0., 10., partyield, partyielderr, partyielderrcorr);

  //  Double_t fint = f->Integral(0.,10.);
  //  Double_t finte = f->IntegralError(0.,10.);
  //  Double_t fmean = f->Mean(0., 10.);

  printf("dN/dy        = %f +- %f (%f)\n", yield, yielderr, yielderrcorr);
  printf("<pt>         = %f +- %f (%f)\n", mean, meanerr, meanerrcorr);
  printf("dN/dy (data) = %f +- %f (%f)\n", partyield[0], partyielderr[0], partyielderrcorr[0]);
  printf("dN/dy (low)  = %f +- %f (%f)\n", partyield[1], partyielderr[1], partyielderrcorr[1]);
  printf("dN/dy (high) = %f +- %f (%f)\n", partyield[2], partyielderr[2], partyielderrcorr[2]);
  //  printf("----\n");
  //  printf("dN/dy (func) = %f +- %f\n", fint, finte);
  //  printf("<pT> (func)  = %f +- %f\n", fmean, 0.);
  
  //  TH1 *hr = (TH1 *)h->Clone("hr");
  //  hr->Divide(f);
  //  new TCanvas;
  //  hr->Draw();

  //  TProfile *p = new TProfile("p", "", 100, 0., 10.);
  //  gROOT->LoadMacro("HistoUtils.C");
  //  HistoUtils_Function2Profile(f, p);
  //  p->Draw();
}

GetYieldAndMean(TH1 *h, TF1 *f, Double_t &yield, Double_t &yielderr, Double_t &yielderrcorr, Double_t &mean, Double_t &meanerr, Double_t &meanerrcorr, Double_t min, Double_t max, Double_t *partyield, Double_t *partyielderr, Double_t *partyielderrcorr)
{
  
  /* find lowest edge in histo */
  Int_t binlo;
  Double_t lo;
  for (Int_t ibin = 1; ibin < h->GetNbinsX() + 1; ibin++) {
    if (h->GetBinContent(ibin) != 0.) {
      binlo = ibin;
      lo = h->GetBinLowEdge(ibin);
      break;
    }
  }
  
  /* find highest edge in histo */
  Int_t binhi;
  Double_t hi;
  for (Int_t ibin = h->GetNbinsX(); ibin > 0; ibin--) {
    if (h->GetBinContent(ibin) != 0.) {
      binhi = ibin + 1;
      hi = h->GetBinLowEdge(ibin + 1);
      break;
    }
  }
  
  /* integrate the data */
  Double_t cont, err, width, cent, integral_data = 0., integralerr_data = 0., integralerrcorr_data = 0., meanintegral_data = 0., meanintegralerr_data = 0., meanintegralerrcorr_data = 0.;
  for (Int_t ibin = binlo; ibin < binhi; ibin++) {
    cent = h->GetBinCenter(ibin);
    width = h->GetBinWidth(ibin);
    cont = h->GetBinContent(ibin);
    err = h->GetBinError(ibin);
    /* check we didn't get an empty bin in between */
    if (cont != 0. && err != 0.) {
      /* all right, use data */
      integral_data += cont * width;
      integralerr_data += err * err * width * width;
      integralerrcorr_data += err * width;
      meanintegral_data += cont * width * cent;
      meanintegralerr_data += err * err * width * width * cent * cent;
      meanintegralerrcorr_data += err * width * cent;
    }
    else {
      /* missing data-point, complain and use function */
      printf("WARNING: missing data-point at %f\n", cent);
      printf("         using function as a patch\n");
      integral_data += f->Integral(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1));
      integralerr_data += f->IntegralError(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1), 0, 0, 1.e-6);
      meanintegral_data += f->Mean(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1)) * f->Integral(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1));
      meanintegralerr_data += f->Mean(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1)) * f->IntegralError(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1), 0, 0, 1.e-6);
    }
  }
  integralerr_data = TMath::Sqrt(integralerr_data);
  meanintegralerr_data = TMath::Sqrt(meanintegralerr_data);
  
  /* integrate below the data */
  Double_t integral_lo = min < lo ? f->Integral(min, lo) : 0.;
  Double_t integralerr_lo = min < lo ? f->IntegralError(min, lo, 0, 0, 1.e-6) : 0.;
  Double_t meanintegral_lo = min < lo ? f->Mean(min, lo) * integral_lo : 0.;
  Double_t meanintegralerr_lo = min < lo ? f->Mean(min, lo) * integralerr_lo : 0.;
  
  /* integrate above the data */
  Double_t integral_hi = max > hi ? f->Integral(hi, max) : 0.;
  Double_t integralerr_hi = max > hi ? f->IntegralError(hi, max, 0, 0, 1.e-6) : 0.;
  Double_t meanintegral_hi = max > hi ? f->Mean(hi, max) * integral_hi : 0.;
  Double_t meanintegralerr_hi = max > hi ? f->Mean(hi, max) * integralerr_hi : 0.;

  /* compute integrated yield */
  yield = integral_data + integral_lo + integral_hi;
  yielderr = TMath::Sqrt(integralerr_data * integralerr_data + 
			 integralerr_lo * integralerr_lo + 
			 integralerr_hi * integralerr_hi);
  yielderrcorr = TMath::Sqrt(integralerrcorr_data * integralerrcorr_data + 
			     integralerr_lo * integralerr_lo + 
			     integralerr_hi * integralerr_hi);
  
  /* compute integrated mean */
  mean = (meanintegral_data + meanintegral_lo + meanintegral_hi) / yield;
  meanerr = TMath::Sqrt(meanintegralerr_data * meanintegralerr_data + 
			meanintegralerr_lo * meanintegralerr_lo + 
			meanintegralerr_hi * meanintegralerr_hi) / yield;
  meanerrcorr = TMath::Sqrt(meanintegralerrcorr_data * meanintegralerrcorr_data + 
			    meanintegralerr_lo * meanintegralerr_lo + 
			    meanintegralerr_hi * meanintegralerr_hi) / yield;

  /* set partial yields */
  partyield[0] = integral_data;
  partyielderr[0] = integralerr_data;
  partyielderrcorr[0] = integralerrcorr_data;
  partyield[1] = integral_lo;
  partyielderr[1] = integralerr_lo;
  partyielderrcorr[1] = integralerr_lo;
  partyield[2] = integral_hi;
  partyielderr[2] = integralerr_hi;
  partyielderrcorr[2] = integralerr_hi;

}

/*****************************************************************/

Double_t 
y2eta(Double_t pt, Double_t mass, Double_t y){
  Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
  return TMath::ASinH(mt / pt * TMath::SinH(y));
}
Double_t 
eta2y(Double_t pt, Double_t mass, Double_t eta){
  Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

TH1 *
Convert_dNdy_1over2pipt_dNdeta(TH1 *hin, Double_t mass, Double_t eta = 0.8)
{

  TH1 *hout = hin->Clone("hout");
  hout->Reset();
  Double_t pt, mt, conv, val, vale;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    pt = hin->GetBinCenter(ibin + 1);
    conv = eta2y(pt, mass, eta) / eta;
    val = hin->GetBinContent(ibin + 1);
    vale = hin->GetBinError(ibin + 1);
    val /= (2. * TMath::Pi() * pt);
    vale /= (2. * TMath::Pi() * pt);
    val *= conv;
    vale *= conv;
    hout->SetBinContent(ibin + 1, val);
    hout->SetBinError(ibin + 1, vale);
  }
  return hout;
}

TH1 *
Convert_dNdy_1over2pipt_dNdy(TH1 *hin)
{

  TH1 *hout = hin->Clone("hout");
  hout->Reset();
  Double_t pt, mt, conv, val, vale;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    pt = hin->GetBinCenter(ibin + 1);
    val = hin->GetBinContent(ibin + 1);
    vale = hin->GetBinError(ibin + 1);
    val /= (2. * TMath::Pi() * pt);
    vale /= (2. * TMath::Pi() * pt);
    hout->SetBinContent(ibin + 1, val);
    hout->SetBinError(ibin + 1, vale);
  }
  return hout;
}

TH1 *
Convert_dNdy_dNdeta(TH1 *hin, Double_t mass, Double_t eta = 0.8)
{

  TH1 *hout = hin->Clone("hout");
  hout->Reset();
  Double_t pt, mt, conv, val, vale;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    pt = hin->GetBinCenter(ibin + 1);
    conv = eta2y(pt, mass, eta) / eta;
    val = hin->GetBinContent(ibin + 1);
    vale = hin->GetBinError(ibin + 1);
    val *= conv;
    vale *= conv;
    hout->SetBinContent(ibin + 1, val);
    hout->SetBinError(ibin + 1, vale);
  }
  return hout;
}

TGraph *
Convert_dNdy_dNdeta(TGraph *hin, Double_t mass, Double_t eta = 0.8)
{

  TGraph *hout = hin->Clone("hout");
  //  hout->Reset();
  Double_t pt, mt, conv, val, valelo, valehi;
  for (Int_t ibin = 0; ibin < hin->GetN(); ibin++) {
    pt = hin->GetX()[ibin];
    conv = eta2y(pt, mass, eta) / eta;
    val = hin->GetY()[ibin];
    valelo = hin->GetEYlow()[ibin];
    valehi = hin->GetEYhigh()[ibin];
    val *= conv;
    valelo *= conv;
    valehi *= conv;
    hout->GetX()[ibin] = pt;
    hout->GetY()[ibin] = val;
    hout->GetEYlow()[ibin] = valelo;
    hout->GetEYhigh()[ibin] = valehi;
  }
  return hout;
}

TH1 *
SummedId_1over2pipt_dNdeta(const Char_t *filename, Int_t icent)
{

  const Char_t *chargeName[2] = {
    "plus", "minus"
  };

  TFile *filein = TFile::Open(filename);
  TH1 *hy[AliPID::kSPECIES][2];
  TH1 *heta[AliPID::kSPECIES][2];
  for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      hy[ipart][icharge] = (TH1 *)filein->Get(Form("cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
      if (!hy[ipart][icharge]) {
	printf("cannot find cent%d_%s_%s\n", icent, AliPID::ParticleName(ipart), chargeName[icharge]);
	return NULL;
      }
      heta[ipart][icharge] = Convert_dNdy_1over2pipt_dNdeta(hy[ipart][icharge], AliPID::ParticleMass(ipart));
    }

  /* sum */
  TH1D *hsum = heta[2][0]->Clone("hsum");
  hsum->Reset();
  for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++)
      hsum->Add(heta[ipart][icharge]);

  return hsum;
}

TH1 *
SummedId_dNdeta(const Char_t *filename, Int_t icent)
{

  const Char_t *chargeName[2] = {
    "plus", "minus"
  };

  TFile *filein = TFile::Open(filename);
  TH1 *hy[AliPID::kSPECIES][2];
  TH1 *heta[AliPID::kSPECIES][2];
  for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      hy[ipart][icharge] = (TH1 *)filein->Get(Form("cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
      if (!hy[ipart][icharge]) {
	printf("cannot find cent%d_%s_%s\n", icent, AliPID::ParticleName(ipart), chargeName[icharge]);
	return NULL;
      }
      heta[ipart][icharge] = Convert_dNdy_dNdeta(hy[ipart][icharge], AliPID::ParticleMass(ipart));
    }

  /* sum */
  TH1D *hsum = heta[2][0]->Clone("hsum");
  hsum->Reset();
  for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++)
      hsum->Add(heta[ipart][icharge]);

  return hsum;
}

/*****************************************************************/

