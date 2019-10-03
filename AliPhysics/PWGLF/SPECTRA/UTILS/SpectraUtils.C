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

  Double_t part1 = (n - 1.) * (n - 2.);
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
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}

/*****************************************************************/
/* TRUE-TSALLIS */
/* See: L. Marques, E. Andrade-II, and A. Deppman  */
/* Phys. Rev. D 87, 114022 – Published 27 June 2013 */
/* for explanation why this is a more correct formulation */
/* Added by Lee Barnby - lbarnby@cern.ch
/*****************************************************************/
Double_t
TrueTsallis_Func(const Double_t *x, Double_t *p)
{
    /* dN/dpt */
    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t q = p[1];
    Double_t T = p[2];
    Double_t dNdy = p[3];
    
    Double_t part1 = pt*mt;
    Double_t part2 = part1/T;
    Double_t part3 = (2.-q)*(3.-2.*q);
    Double_t part4 = (2.-q)*mass*mass + 2.*mass*T +2.*T*T;
    Double_t part5 = part3/part4;
    
    Double_t part6 = 1. + (q-1.)*mass/T;
    Double_t part7 = TMath::Power(part6, 1./(q-1.));
    
    Double_t part8 = 1. + (q-1.)*mt/T;
    Double_t part9 = TMath::Power(part8, -q/(q-1.));
    return part2 * dNdy * part5 * part7 * part9;
}

TF1*
TrueTsallis(const Char_t *name, Double_t mass, Double_t q = 1.1., Double_t T = 0.1, Double_t dNdy = 0.2)
{
    TF1 *fTrueTsallis = new TF1(name, TrueTsallis_Func, 0., 10., 4);
    fTrueTsallis->SetParameters(mass, q, T, dNdy);
    fTrueTsallis->SetParNames("mass", "q", "T", "dN/dy");
    fTrueTsallis->FixParameter(0, mass);
    fTrueTsallis->SetParLimits(1, 1., 2.);
    fTrueTsallis->SetParLimits(2, 1.e-3, 10.e3);
    fTrueTsallis->SetParLimits(3, 1.e-8, 1.e8);
    return fTrueTsallis;
    
}

/*****************************************************************/
/* BOLTZMANN-GIBBS BLAST-WAVE */
/*****************************************************************/

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;
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
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., (Double_t *)0, 1.e-6);
  return norm * pt * integral;
}

Double_t
BGBlastWaveRatio_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max_num = p[1];
  Double_t temp_num = p[2];
  Double_t n_num = p[3];
  Double_t norm_num = p[4];
  Double_t beta_max_den = p[5];
  Double_t temp_den = p[6];
  Double_t n_den = p[7];
  Double_t norm_den = p[8];
  
  if (!fBGBlastWave_Integrand_num)
    fBGBlastWave_Integrand_num = new TF1("fBGBlastWave_Integrand_num", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_num->SetParameters(mt, pt, beta_max_num, temp_num, n_num);
  Double_t integral_num = fBGBlastWave_Integrand_num->Integral(0., 1.);

  if (!fBGBlastWave_Integrand_den)
    fBGBlastWave_Integrand_den = new TF1("fBGBlastWave_Integrand_den", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_den->SetParameters(mt, pt, beta_max_den, temp_den, n_den);
  Double_t integral_den = fBGBlastWave_Integrand_den->Integral(0., 1.);

  return (norm_num / norm_den) * (integral_num / integral_den);
}

Double_t
BGBlastWaveParticleRatio_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass_num = p[0];
  Double_t mass_den = p[1];
  Double_t mt_num = TMath::Sqrt(pt * pt + mass_num * mass_num);
  Double_t mt_den = TMath::Sqrt(pt * pt + mass_den * mass_den);
  Double_t beta_max = p[2];
  Double_t temp = p[3];
  Double_t n = p[4];
  Double_t norm_num = p[5];
  Double_t norm_den = p[6];
  
  if (!fBGBlastWave_Integrand_num)
    fBGBlastWave_Integrand_num = new TF1("fBGBlastWave_Integrand_num", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_num->SetParameters(mt_num, pt, beta_max, temp, n);
  Double_t integral_num = fBGBlastWave_Integrand_num->Integral(0., 1.);
  
  if (!fBGBlastWave_Integrand_den)
    fBGBlastWave_Integrand_den = new TF1("fBGBlastWave_Integrand_den", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_den->SetParameters(mt_den, pt, beta_max, temp, n);
  Double_t integral_den = fBGBlastWave_Integrand_den->Integral(0., 1.);

  return (norm_num / norm_den) * (integral_num / integral_den);
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
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., (Double_t *)0, 1.e-3);

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
  fBGBlastWave->SetParLimits(3, 0.01, 50.);
  return fBGBlastWave;
}

TF1 *
BGBlastWaveRatio(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWaveRatio_Func, 0., 10., 9);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max_num", "T_num", "n_num", "norm_num", "beta_max_den", "T_den", "n_den", "norm_den");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 10.);
  fBGBlastWave->SetParLimits(5, 0.01, 0.99);
  fBGBlastWave->SetParLimits(6, 0.01, 1.);
  fBGBlastWave->SetParLimits(7, 0.01, 10.);
  return fBGBlastWave;
}

TF1 *
BGBlastWaveParticleRatio(const Char_t *name, Double_t mass_num, Double_t mass_den, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm_num = 1.e6, Double_t norm_den = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWaveParticleRatio_Func, 0., 10., 7);
  fBGBlastWave->SetParameters(mass_num, mass_den, beta_max, temp, n, norm_num, norm_den);
  fBGBlastWave->SetParNames("mass_num", "mass_den", "beta_max", "T", "n", "norm_num", "norm_den");
  fBGBlastWave->FixParameter(0, mass_num);
  fBGBlastWave->FixParameter(1, mass_den);
  fBGBlastWave->SetParLimits(2, 0.01, 0.99);
  fBGBlastWave->SetParLimits(3, 0.01, 1.);
  fBGBlastWave->SetParLimits(4, 0.01, 10.);
  return fBGBlastWave;
}

TF1 *BGBlastWave_OneOverPT(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func_OneOverPt, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 50.);
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
TF1 *fBGBWratio[1000];
TGraphErrors *gBW[1000];

TObjArray *
BGBlastWave_GlobalFit(TObjArray *data, Double_t *mass, Double_t profile = 0.5, Bool_t computeCont = kFALSE, Bool_t fixProfile = kFALSE)
{

  /* get data */
  Int_t ndf = 0;
  nBW = data->GetEntries();
  for (Int_t idata = 0; idata < nBW; idata++) {
    gBW[idata] = (TGraphErrors *)data->At(idata);
    gBW[idata]->SetName(Form("gBW%d", idata));
    ndf += gBW[idata]->GetN();
  }

  /* init BG blast-wave functions */
  for (Int_t idata = 0; idata < nBW; idata++) {
    printf("init BG-BlastWave function #%d: mass = %f\n", idata, mass[idata]);
    fBGBW[idata] = BGBlastWave(Form("fBGBW%d", idata), mass[idata]);
  }

  if (computeCont)
    printf("-> compute contours requested\n");

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
  for (Int_t idata = 0; idata < nBW; idata++) {
    minuit->mnparm(idata, Form("norm%d", idata), 1.e6, 1., 0., 0., ierflg);
    ndf--;
  }
  //  minuit->mnparm(nBW + 0, "<beta>", 0.55, 0.01, 0., 1., ierflg);
  //  minuit->mnparm(nBW + 1, "T", 0.14, 0.01, 0., 1., ierflg);
  //  minuit->mnparm(nBW + 2, "n", profile, 0.1, 0., 10., ierflg);
  
  minuit->mnparm(nBW + 0, "<beta>", 0.7, 0.01, 0.2, 0.7, ierflg);
  minuit->mnparm(nBW + 1, "T", 0.07, 0.001, 0.07, 0.2, ierflg);
  minuit->mnparm(nBW + 2, "n", profile, 0.1, 0.6, 5., ierflg);

  ndf -= 3;
  if (fixProfile) {
    minuit->FixParameter(nBW + 2);
    ndf++;
  }

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
  printf("[x] *********************************\n");
  printf("[x] beta_max = %f\n", beta_max);
  printf("[x] <beta>   = %f +- %f (e+ = %f, e- = %f)\n", beta, betae, betaeplus, betaeminus);
  printf("[x] T        = %f +- %f (e+ = %f, e- = %f)\n", temp, tempe, tempeplus, tempeminus);
  printf("[x] n        = %f +- %f (e+ = %f, e- = %f)\n", prof, profe, profeplus, profeminus);
  printf("[x] chi2     = %f\n", amin);
  printf("[x] ndf      = %f\n", ndf);

  /* 1-sigma contour */
  minuit->SetErrorDef(1);
  TGraph *gCont1 = NULL;
  if (computeCont) gCont1 = (TGraph *) minuit->Contour(50, nBW + 0, nBW + 1);
  if (gCont1) gCont1->SetName("gCont1");

  /* 2-sigma contour */
  minuit->SetErrorDef(4);
  TGraph *gCont2 = NULL;
  //  if (computeCont) gCont2 = (TGraph *) minuit->Contour(50, nBW + 0, nBW + 1);
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
  TH1D *hBW = new TH1D("hBW", "", 4, 0., 4.);
  hBW->SetBinContent(1, beta);
  hBW->SetBinError(1, betae);
  hBW->SetBinContent(2, temp);
  hBW->SetBinError(2, tempe);
  hBW->SetBinContent(3, prof);
  hBW->SetBinError(3, profe);
  hBW->SetBinContent(4, amin/ndf);

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
#if 1
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

TObjArray *
BGBlastWave_GlobalFitRatio(TObjArray *data, Double_t *mass, Double_t profile = .7, Bool_t fixProfile = kFALSE)
{

  /* get data */
  Int_t ndf = 0;
  nBW = data->GetEntries();
  for (Int_t idata = 0; idata < nBW; idata++) {
    gBW[idata] = (TGraphErrors *)data->At(idata);
    gBW[idata]->SetName(Form("gBW%d", idata));
    ndf += gBW[idata]->GetN();
  }

  /* init BG blast-wave functions */
  for (Int_t idata = 0; idata < nBW; idata++) {
    printf("init BG-BlastWaveRatio function #%d: mass = %f\n", idata, mass[idata]);
    fBGBWratio[idata] = BGBlastWaveRatio(Form("fBGBWratio%d", idata), mass[idata]);
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
  const Int_t nfitpars = 2 * (nBW + nbwpars);
  TMinuit *minuit = new TMinuit(nfitpars);
  minuit->SetFCN(BGBlastWave_FCNRatio);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  for (Int_t idata = 0; idata < nBW; idata++) {
    minuit->mnparm(idata, Form("norm%d_num", idata), 1.e6, 1., 0., 0., ierflg);
    ndf--;
  }
  for (Int_t idata = nBW; idata < 2 * nBW; idata++) {
    minuit->mnparm(idata, Form("norm%d_den", idata), 1.e6, 1., 0., 0., ierflg);
    minuit->FixParameter(idata);
    ndf--;
  }
  minuit->mnparm(2 * nBW + 0, "<beta>_num", 0.65, 0.01, 0., 1., ierflg);
  minuit->mnparm(2 * nBW + 1, "T_num", 0.1, 0.01, 0., 1., ierflg);
  minuit->mnparm(2 * nBW + 2, "n_num", profile, 0.1, 0., 10., ierflg);
  minuit->mnparm(2 * nBW + 3, "<beta>_den", 0.65, 0.01, 0., 1., ierflg);
  minuit->mnparm(2 * nBW + 4, "T_den", 0.1, 0.01, 0., 1., ierflg);
  minuit->mnparm(2 * nBW + 5, "n_den", profile, 0.1, 0., 10., ierflg);
  ndf -= 3;

  if (fixProfile) {
    minuit->FixParameter(nBW + 2);
    minuit->FixParameter(nBW + 5);
    ndf++;
  }

  /* set strategy */
  arglist[0] = 1;
  minuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);

  printf("-->> STARTING MIGRAD with %d parameters\n", nfitpars);

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
  arglist[1] = 2 * nBW + 1;
  arglist[2] = 2 * nBW + 2;
  arglist[3] = 2 * nBW + 3;
  arglist[4] = 2 * nBW + 4;
  arglist[5] = 2 * nBW + 5;
  arglist[6] = 2 * nBW + 6;
  minuit->mnexcm("MINOS", arglist, 7, ierflg);

  /* print results */
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  minuit->mnprin(4, amin);

  /* get parameters */
  Double_t beta_num, betae_num, betaeplus_num, betaeminus_num, betagcc_num, temp_num, tempe_num, tempeplus_num, tempeminus_num, tempgcc_num, prof_num, profe_num, profeplus_num, profeminus_num, profgcc_num;
  Double_t beta_den, betae_den, betaeplus_den, betaeminus_den, betagcc_den, temp_den, tempe_den, tempeplus_den, tempeminus_den, tempgcc_den, prof_den, profe_den, profeplus_den, profeminus_den, profgcc_den;
  minuit->GetParameter(2*nBW + 0, beta_num, betae_num);
  minuit->mnerrs(nBW + 0, betaeplus_num, betaeminus_num, betae_num, betagcc_num);
  minuit->GetParameter(2*nBW + 1, temp_num, tempe_num);
  minuit->mnerrs(nBW + 1, tempeplus_num, tempeminus_num, tempe_num, tempgcc_num);
  minuit->GetParameter(2*nBW + 2, prof_num, profe_num);
  minuit->mnerrs(nBW + 2, profeplus_num, profeminus_num, profe_num, profgcc_num);
  minuit->GetParameter(2*nBW + 3, beta_den, betae_den);
  minuit->mnerrs(nBW + 3, betaeplus_den, betaeminus_den, betae_den, betagcc_den);
  minuit->GetParameter(2*nBW + 4, temp_den, tempe_den);
  minuit->mnerrs(nBW + 4, tempeplus_den, tempeminus_den, tempe_den, tempgcc_den);
  minuit->GetParameter(2*nBW + 5, prof_den, profe_den);
  minuit->mnerrs(nBW + 5, profeplus_den, profeminus_den, profe_den, profgcc_den);
  Double_t beta_max_num, beta_max_den;
  beta_max_num = 0.5 * (2. + prof_num) * beta_num;
  beta_max_den = 0.5 * (2. + prof_den) * beta_den;
  Double_t norm_num[1000], norme_num[1000];
  Double_t norm_den[1000], norme_den[1000];
  for (Int_t idata = 0; idata < nBW; idata++)
    minuit->GetParameter(idata, norm_num[idata], norme_num[idata]);
  for (Int_t idata = 0; idata < nBW; idata++) {
    minuit->GetParameter(nBW + idata, norm_den[idata], norme_den[idata]);
  }

  /* printout */
  printf("[x] *********************************\n");
  printf("[x] beta_max = %f\n", beta_max_num);
  printf("[x] <beta>   = %f +- %f (e+ = %f, e- = %f)\n", beta_num, betae_num, betaeplus_num, betaeminus_num);
  printf("[x] T        = %f +- %f (e+ = %f, e- = %f)\n", temp_num, tempe_num, tempeplus_num, tempeminus_num);
  printf("[x] n        = %f +- %f (e+ = %f, e- = %f)\n", prof_num, profe_num, profeplus_num, profeminus_num);
  printf("[x] *********************************\n");
  printf("[x] beta_max = %f\n", beta_max_den);
  printf("[x] <beta>   = %f +- %f (e+ = %f, e- = %f)\n", beta_den, betae_den, betaeplus_den, betaeminus_den);
  printf("[x] T        = %f +- %f (e+ = %f, e- = %f)\n", temp_den, tempe_den, tempeplus_den, tempeminus_den);
  printf("[x] n        = %f +- %f (e+ = %f, e- = %f)\n", prof_den, profe_den, profeplus_den, profeminus_den);
  printf("[x] *********************************\n");
  printf("[x] chi2     = %f\n", amin);
  printf("[x] ndf      = %f\n", ndf);

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
    fBGBWratio[idata]->SetParameter(4, norm_num[idata]);
    fBGBWratio[idata]->SetParameter(1, beta_max_num);
    fBGBWratio[idata]->SetParameter(2, temp_num);
    fBGBWratio[idata]->SetParameter(3, prof_num);
    fBGBWratio[idata]->SetParameter(8, norm_den[idata]);
    fBGBWratio[idata]->SetParameter(5, beta_max_den);
    fBGBWratio[idata]->SetParameter(6, temp_den);
    fBGBWratio[idata]->SetParameter(7, prof_den);
    fBGBWratio[idata]->Draw("same");
  }
  cBW->Update();

  return;

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
    outoa->Add(fBGBWratio[idata]);
  }
  outoa->Add(cBW);
  outoa->Add(hBW);
  outoa->Add(gBetaT);
  if (gCont1) outoa->Add(gCont1);
  if (gCont2) outoa->Add(gCont2);

  return outoa;

}

void 
BGBlastWave_FCNRatio(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  /* beta -> beta_max */
  Double_t beta_num = par[2*nBW+0];
  Double_t T_num = par[2*nBW+1];
  Double_t n_num = par[2*nBW+2];
  Double_t beta_max_num = 0.5 * (2. + n_num) * beta_num;
  Double_t beta_den = par[2*nBW+3];
  Double_t T_den = par[2*nBW+4];
  Double_t n_den = par[2*nBW+5];
  Double_t beta_max_den = 0.5 * (2. + n_den) * beta_den;
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
    fBGBWratio[iBW]->SetParameter(4, par[iBW]);
    fBGBWratio[iBW]->SetParameter(1, beta_max_num);
    fBGBWratio[iBW]->SetParameter(2, T_num);
    fBGBWratio[iBW]->SetParameter(3, n_num);
    fBGBWratio[iBW]->SetParameter(8, par[nBW+iBW]);
    fBGBWratio[iBW]->SetParameter(5, beta_max_den);
    fBGBWratio[iBW]->SetParameter(6, T_den);
    fBGBWratio[iBW]->SetParameter(7, n_den);
    /* loop over all the points */
    for (Int_t ipt = 0; ipt < gBW[iBW]->GetN(); ipt++) {
      pt = gBW[iBW]->GetX()[ipt];
      pte = gBW[iBW]->GetEX()[ipt];
      val = gBW[iBW]->GetY()[ipt];
      vale = gBW[iBW]->GetEY()[ipt];
      func = fBGBWratio[iBW]->Eval(pt);
      //      func = fBGBW[iBW]->Integral(pt - pte, pt + pte);
      pull = (val - func) / vale;
      chi += pull * pull;
    }
  }

  f = chi;
}

/*****************************************************************/

Int_t ratioPartNum[9] = {2, 3, 4, 3, 3, 3, 4, 4, 4};
Int_t ratioChargeNum[9] = {1, 1, 1, 0, 1, 2, 0, 1, 2};
Int_t ratioPartDen[9] = {2, 3, 4, 2, 2, 2, 2, 2, 2};
Int_t ratioChargeDen[9] = {0, 0, 0, 0, 1, 2, 0, 1, 2};

Char_t *ratioPartNumName[9] = {"pion", "kaon", "proton", "kaon", "kaon", "kaon", "proton", "proton", "proton"};
Char_t *ratioChargeNumName[9] = {"minus", "minus", "minus", "plus", "minus", "both", "plus", "minus", "both"};
Char_t *ratioPartDenName[9] = {"pion", "kaon", "proton", "pion", "pion", "pion", "pion", "pion", "pion"};
Char_t *ratioChargeDenName[9] = {"plus", "plus", "plus", "plus", "minus", "both", "plus", "minus", "both"};

IntegratedRatioError(const Char_t *spectrafilename, const Char_t *ratiosfilename)
{

  TFile *fspectra = TFile::Open(spectrafilename);
  TFile *fratios = TFile::Open(ratiosfilename);
  
  Int_t icent = 0;
  for (Int_t irat = 3; irat < 9; irat++) {

    Int_t ipnum = ratioPartNum[irat];
    Int_t icnum = ratioChargeNum[irat];

    Int_t ipden = ratioPartDen[irat];
    Int_t icden = ratioChargeDen[irat];

    if (icnum == 2 || icden == 2) continue;

    printf("%s\n", Form("sys_cent%d_%s_%s_%s_%s", icent, ratioPartNumName[irat], ratioChargeNumName[irat], ratioPartDenName[irat], ratioChargeDenName[irat]));
    TH1 *hrat_sys = (TH1 *)fratios->Get(Form("sys_cent%d_%s_%s_%s_%s", icent, ratioPartNumName[irat], ratioChargeNumName[irat], ratioPartDenName[irat], ratioChargeDenName[irat]));

    printf("%d %d %d %d\n", ipnum, icnum, ipden, icden);
    TH1 *hnum = (TH1 *)fspectra->Get(Form("cent%d_%s_%s", icent, ratioPartNumName[irat], ratioChargeNumName[irat]));
    TH1 *hnum_stat = (TH1 *)fspectra->Get(Form("stat_cent%d_%s_%s", icent, ratioPartNumName[irat], ratioChargeNumName[irat]));
    printf("%s = %p\n", Form("cent%d_%s_%s", icent, ratioPartNumName[irat], ratioChargeNumName[irat]), hnum);
  
    TH1 *hden = (TH1 *)fspectra->Get(Form("cent%d_%s_%s", icent, ratioPartDenName[irat], ratioChargeDenName[irat]));
    TH1 *hden_stat = (TH1 *)fspectra->Get(Form("stat_cent%d_%s_%s", icent, ratioPartDenName[irat], ratioChargeDenName[irat]));
    printf("%s = %p\n", Form("cent%d_%s_%s", icent, ratioPartDenName[irat], ratioChargeDenName[irat]), hden);
    
    TH1 *hnum_plus = (TH1 *)hnum->Clone("hnum_plus");
    TH1 *hnum_minus = (TH1 *)hnum->Clone("hnum_minus");
    TH1 *hden_plus = (TH1 *)hden->Clone("hden_plus");
    TH1 *hden_minus = (TH1 *)hden->Clone("hden_minus");
    TH1 *hnum_stat_plus = (TH1 *)hnum_stat->Clone("hnum_stat_plus");
    TH1 *hnum_stat_minus = (TH1 *)hnum_stat->Clone("hnum_stat_minus");
    TH1 *hden_stat_plus = (TH1 *)hden_stat->Clone("hden_stat_plus");
    TH1 *hden_stat_minus = (TH1 *)hden_stat->Clone("hden_stat_minus");
    for (Int_t ibin = 0; ibin < hnum_plus->GetNbinsX(); ibin++) {
      Double_t val = hrat_sys->GetBinContent(ibin + 1);
      Double_t vale = hrat_sys->GetBinError(ibin + 1);
      if (vale <= 0.) continue;
      Double_t syse = vale / val;

      val = hnum_plus->GetBinContent(ibin + 1);
      vale = hnum_plus->GetBinError(ibin + 1);
      val *= (1. + syse);
      vale *= (1. + syse);
      hnum_plus->SetBinContent(ibin + 1, val);
      hnum_plus->SetBinError(ibin + 1, vale);

      val = hnum_minus->GetBinContent(ibin + 1);
      vale = hnum_minus->GetBinError(ibin + 1);
      val *= (1. - syse);
      vale *= (1. - syse);
      hnum_minus->SetBinContent(ibin + 1, val);
      hnum_minus->SetBinError(ibin + 1, vale);

      val = hden_plus->GetBinContent(ibin + 1);
      vale = hden_plus->GetBinError(ibin + 1);
      val *= (1. + syse);
      vale *= (1. + syse);
      hden_plus->SetBinContent(ibin + 1, val);
      hden_plus->SetBinError(ibin + 1, vale);

      val = hden_minus->GetBinContent(ibin + 1);
      vale = hden_minus->GetBinError(ibin + 1);
      val *= (1. - syse);
      vale *= (1. - syse);
      hden_minus->SetBinContent(ibin + 1, val);
      hden_minus->SetBinError(ibin + 1, vale);



      val = hnum_stat_plus->GetBinContent(ibin + 1);
      vale = hnum_stat_plus->GetBinError(ibin + 1);
      val *= (1. + syse);
      vale *= (1. + syse);
      hnum_stat_plus->SetBinContent(ibin + 1, val);
      hnum_stat_plus->SetBinError(ibin + 1, vale);

      val = hnum_stat_minus->GetBinContent(ibin + 1);
      vale = hnum_stat_minus->GetBinError(ibin + 1);
      val *= (1. - syse);
      vale *= (1. - syse);
      hnum_stat_minus->SetBinContent(ibin + 1, val);
      hnum_stat_minus->SetBinError(ibin + 1, vale);

      val = hden_stat_plus->GetBinContent(ibin + 1);
      vale = hden_stat_plus->GetBinError(ibin + 1);
      val *= (1. + syse);
      vale *= (1. + syse);
      hden_stat_plus->SetBinContent(ibin + 1, val);
      hden_stat_plus->SetBinError(ibin + 1, vale);

      val = hden_stat_minus->GetBinContent(ibin + 1);
      vale = hden_stat_minus->GetBinError(ibin + 1);
      val *= (1. - syse);
      vale *= (1. - syse);
      hden_stat_minus->SetBinContent(ibin + 1, val);
      hden_stat_minus->SetBinError(ibin + 1, vale);
    }

    TCanvas *c = new TCanvas;
    hnum_stat->SetLineColor(8);
    hnum_stat->SetMarkerColor(8);
    hnum_stat->SetMarkerSize(1);
    hnum_stat->SetMarkerStyle(20);
    hnum_stat->Draw();
    
    hnum_stat_plus->SetLineColor(2);
    hnum_stat_plus->SetMarkerColor(2);
    hnum_stat_plus->SetMarkerSize(1);
    hnum_stat_plus->SetMarkerStyle(20);
    hnum_stat_plus->Draw("same");
    
    hnum_stat_minus->SetLineColor(4);
    hnum_stat_minus->SetMarkerColor(4);
    hnum_stat_minus->SetMarkerSize(1);
    hnum_stat_minus->SetMarkerStyle(20);
    hnum_stat_minus->Draw("same");
    
    c->Update();

    c = new TCanvas();
    TH1 *hrat = (TH1 *)hnum_stat->Clone("hrat");
    TH1 *hrat_plus = (TH1 *)hnum_stat_plus->Clone("hrat_plus");
    TH1 *hrat_minus = (TH1 *)hnum_stat_minus->Clone("hrat_minus");

    hrat->Divide(hden_stat);
    hrat_plus->Divide(hden_stat);
    hrat_minus->Divide(hden_stat);

    hrat->Draw();
    hrat_plus->Draw("same");
    hrat_minus->Draw("same");


    c->Update();
    
    getchar();
    
    
    /* integrate as they are */
    printf("***** NUM:\n");
    IntegratedProduction(hnum, AliPID::ParticleMass(ratioPartNum[irat]), 0, "0q");
    printf("***** NUM PLUS:\n");
    IntegratedProduction(hnum_plus, AliPID::ParticleMass(ratioPartNum[irat]), 0, "0q");
    printf("***** NUM MINUS:\n");
    IntegratedProduction(hnum_minus, AliPID::ParticleMass(ratioPartNum[irat]), 0, "0q");
    
    printf("***** DEN:\n");
    IntegratedProduction(hden, AliPID::ParticleMass(ratioPartDen[irat]), 0, "0q");
    printf("***** DEN PLUS:\n");
    IntegratedProduction(hden_plus, AliPID::ParticleMass(ratioPartDen[irat]), 0, "0q");
    printf("***** DEN MINUS:\n");
    IntegratedProduction(hden_minus, AliPID::ParticleMass(ratioPartDen[irat]), 0, "0q");

  }
  
}

enum EIntegratedData_t {
  kInt_YieldData, kInt_YieldDataErr, kInt_YieldDataLin,
  kInt_MeanData, kInt_MeanDataErr, kInt_MeanDataLin,
  kInt_YieldLow, kInt_YieldLowErr,
  kInt_MeanLow, kInt_MeanLowErr,
  kInt_YieldHigh, kInt_YieldHighErr,
  kInt_MeanHigh, kInt_MeanHighErr,
  kInt_NData
};

IntegratedProduction_measurement(const Char_t *filename, Option_t *opt = "q0R")
{

  for (Int_t ipart = 2; ipart < 5; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      IntegratedProduction(filename, ipart, icharge, 0, 0., 10., opt);
      IntegratedProduction(filename, ipart, icharge, 1, 0., 10., opt);
    }

}

IntegratedProduction_systematics(const Char_t *filename, Option_t *opt = "q0R")
{

  Double_t min[5] = {0., 0., 0., 0., 0.};
  Double_t max[5] = {0., 0., 0.5, 1.0, 2.0};
  for (Int_t ipart = 2; ipart < 5; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++)
      for (Int_t ifunc = 1; ifunc < 6; ifunc++)
	IntegratedProduction(filename, ipart, icharge, ifunc, min[ipart], max[ipart], opt);

}

IntegratedProduction_check(const Char_t *filename, Option_t *opt = "q0R")
{

  for (Int_t ipart = 2; ipart < 5; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++)
      for (Int_t ifunc = 1; ifunc < 6; ifunc++)
	IntegratedProduction(filename, ipart, icharge, ifunc, 0., 10., opt);

}

IntegratedProduction(const Char_t *filename, Int_t ipart, Int_t icharge, Int_t ifunc = 0, Float_t min = 0., Float_t max = 10., Option_t *opt = "q0R")
{

  Char_t *centName[10] = {
    "0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-100%", "minimum bias"
  };
  
  
  Char_t *chargeName[3] = {"plus", "minus", "sum"};
  
  Char_t *partChargeName[5][2] = {
    "", "",
    "", "",
    "#pi^{+}", "#pi^{-}",
    "K^{+}", "K^{-}",
    "p", "#bar{p}"
  };
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  /* PWGTools */
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  AliPWGFunc pwgfunc;
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);

  Double_t mass = AliPID::ParticleMass(ipart);

  TF1 *f = NULL;
  switch (ifunc) {
  case 0:
    f = BGBlastWave("BlastWave", mass);
    f->SetLineColor(2);
    break;
  case 1:
    f = LevyTsallis("Levy-Tsallis", mass, 5., mass);
    f->SetLineColor(8);
    break;
  case 2:
    f = Boltzmann("Boltzmann", mass);
    f->SetLineColor(kPink+1);
    break;
  case 3:
    f = pwgfunc.GetMTExp(mass, 0.1, 1., "m_{T}-exponential");
    f->SetLineColor(kViolet+1);
    break;
  case 4:
    f = pwgfunc.GetPTExp(0.1, 1., "p_{T}-exponential");
    f->SetLineColor(kAzure+1);
    break;
  case 5:
    f = pwgfunc.GetBoseEinstein(mass, 0.1, 1., "Bose-Einstein");
    f->SetLineColor(kYellow+1);
    break;
  default:
  };
  f->SetTitle(f->GetName());
  f->SetLineWidth(2);
  f->SetFillColor(0);
  f->SetMarkerColor(f->GetLineColor());

  /* define low-pt fit range according to the number of free parameters */
  Int_t nfreeparams = f->GetNumberFreeParameters();

  TCanvas *cc = new TCanvas("cCanvasFit", "", 1200, 800);
  cc->Divide(4, 2);

  //  if (max < 0)
  //  TFile *fout = TFile::Open(Form("%s.IntegratedProduction.LowPtFit_%s_%s_%s.root", filename, f->GetName(), AliPID::ParticleName(ipart), chargeName[icharge]), "RECREATE");
  //  else if (min < 0)
  //  TFile *fout = TFile::Open(Form("%s.IntegratedProduction.HighPtFit_%s_%s_%s.root", filename, f->GetName(), AliPID::ParticleName(ipart), chargeName[icharge]), "RECREATE");
  //  else
    TFile *fout = TFile::Open(Form("%s.IntegratedProduction_%s_%s_%s.root", filename, f->GetName(), AliPID::ParticleName(ipart), chargeName[icharge]), "RECREATE");
  
  Double_t integrated_data[kInt_NData];

  TFile *filein = TFile::Open(filename);
  TH1F *hyield_data_stat = new TH1F(Form("hPartYield_data_stat_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_data_sys = new TH1F(Form("hPartYield_data_sys_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_data = new TH1F(Form("hPartYield_data_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_low_stat = new TH1F(Form("hPartYield_low_stat_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_low_sys = new TH1F(Form("hPartYield_low_sys_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_low = new TH1F(Form("hPartYield_low_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_high_stat = new TH1F(Form("hPartYield_high_stat_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_high_sys = new TH1F(Form("hPartYield_high_sys_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_high = new TH1F(Form("hPartYield_high_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_data_stat = new TH1F(Form("hPartMean_data_stat_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_data_sys = new TH1F(Form("hPartMean_data_sys_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_data = new TH1F(Form("hPartMean_data_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_low_stat = new TH1F(Form("hPartMean_low_stat_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_low_sys = new TH1F(Form("hPartMean_low_sys_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_low = new TH1F(Form("hPartMean_low_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_high_stat = new TH1F(Form("hPartMean_high_stat_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_high_sys = new TH1F(Form("hPartMean_high_sys_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_high = new TH1F(Form("hPartMean_high_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_stat = new TH1F(Form("hYield_stat_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield_sys = new TH1F(Form("hYield_sys_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hyield = new TH1F(Form("hYield_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_stat = new TH1F(Form("hMean_stat_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean_sys = new TH1F(Form("hMean_sys_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  TH1F *hmean = new TH1F(Form("hMean_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 7, 0, 7);
  for (Int_t icent = 0; icent < 7; icent++) {
    cc->cd(icent + 1)->SetLogy();
    printf("*************************\n");
    printf("cent=%d part=%d charge=%d\n", icent, ipart, icharge);
    TH1 *h = (TH1 *)filein->Get(Form("cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    //    h = Convert_dNdy_dNdeta(h, AliPID::ParticleMass(ipart), 0.5);
    h->SetTitle(Form("%s (%s);p_{T} (GeV/c);dN/dp_{T}", partChargeName[ipart][icharge], centName[icent]));
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1);
    h->SetLineWidth(1);
    h->SetFillColor(0);
    h->SetFillStyle(0);
    h->SetLineColor(1);
    h->SetMarkerColor(4);
    TH1 *hstat = (TH1 *)filein->Get(Form("stat_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    //    hstat = Convert_dNdy_dNdeta(hstat, AliPID::ParticleMass(ipart), 0.5);
    hstat->SetMarkerStyle(20);
    hstat->SetMarkerSize(1);
    hstat->SetLineWidth(1);
    hstat->SetLineColor(1);
    hstat->SetMarkerColor(4);
    TH1 *hsys = (TH1 *)filein->Get(Form("sys_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    //    hsys = Convert_dNdy_dNdeta(hsys, AliPID::ParticleMass(ipart), 0.5);
    hsys->SetMarkerStyle(20);
    hsys->SetMarkerSize(1);
    hsys->SetLineWidth(1);
    hsys->SetLineColor(1);
    hsys->SetMarkerColor(4);

    if (max < 0.) {
      /* define fit range according to number of free parameter */
      Int_t gotpoints = 0;
      for (Int_t ipt = 0; ipt < h->GetNbinsX(); ipt++) {
	if (h->GetBinError(ipt + 1) <= 0.) continue;
	gotpoints++;
	if (gotpoints < nfreeparams*2) continue;
	max = h->GetXaxis()->GetBinUpEdge(ipt + 1);
	break;
      }
      printf("fitting up to %f\n", max);
    }    
    if (min < 0.) {
      /* define fit range according to number of free parameter */
      Int_t gotpoints = 0;
      for (Int_t ipt = h->GetNbinsX(); ipt > 0; ipt--) {
	if (h->GetBinError(ipt) <= 0.) continue;
	gotpoints++;
	if (gotpoints < nfreeparams+3) continue;
	min = h->GetXaxis()->GetBinLowEdge(ipt);
	break;
      }
      printf("fitting from %f\n", min);
    }    

    printf("-> processing stat-only data\n");

    IntegratedProduction(hstat, f, opt, min, max, integrated_data);
    fout->cd();
    hstat->Write(Form("IntegratedProduction_stat_cent%d_%s_%s.root", icent, AliPID::ParticleName(ipart), chargeName[icharge]));

    hyield_data_stat->SetBinContent(icent + 1, integrated_data[kInt_YieldData]);
    hyield_data_stat->SetBinError(icent + 1, integrated_data[kInt_YieldDataErr]);
    hyield_low_stat->SetBinContent(icent + 1, integrated_data[kInt_YieldLow]);
    hyield_low_stat->SetBinError(icent + 1, integrated_data[kInt_YieldLowErr]);

    hyield_high_stat->SetBinContent(icent + 1, integrated_data[kInt_YieldHigh]);
    hyield_high_stat->SetBinError(icent + 1, integrated_data[kInt_YieldHighErr]);

    hmean_data_stat->SetBinContent(icent + 1, integrated_data[kInt_MeanData]);
    hmean_data_stat->SetBinError(icent + 1, integrated_data[kInt_MeanDataErr]);

    hmean_low_stat->SetBinContent(icent + 1, integrated_data[kInt_MeanLow]);
    hmean_low_stat->SetBinError(icent + 1, integrated_data[kInt_MeanLowErr]);
    
    printf("-> processing sys-only data\n");

    IntegratedProduction(hsys, f, opt, min, max, integrated_data);
    fout->cd();
    hsys->Write(Form("IntegratedProduction_sys_cent%d_%s_%s.root", icent, AliPID::ParticleName(ipart), chargeName[icharge]));

    hyield_data_sys->SetBinContent(icent + 1, integrated_data[kInt_YieldData]);
    hyield_data_sys->SetBinError(icent + 1, integrated_data[kInt_YieldDataLin]);

    hyield_low_sys->SetBinContent(icent + 1, integrated_data[kInt_YieldLow]);
    hyield_low_sys->SetBinError(icent + 1, integrated_data[kInt_YieldLowErr]);

    hyield_high_sys->SetBinContent(icent + 1, integrated_data[kInt_YieldHigh]);
    hyield_high_sys->SetBinError(icent + 1, integrated_data[kInt_YieldHighErr]);

    hmean_data_sys->SetBinContent(icent + 1, integrated_data[kInt_MeanData]);
    hmean_data_sys->SetBinError(icent + 1, integrated_data[kInt_MeanDataLin]);

    hmean_low_sys->SetBinContent(icent + 1, integrated_data[kInt_MeanLow]);
    hmean_low_sys->SetBinError(icent + 1, integrated_data[kInt_MeanLowErr]);
    
    printf("-> processing stat+sys data\n");

    IntegratedProduction(h, f, opt, min, max, integrated_data);
    fout->cd();
    h->Write(Form("IntegratedProduction_cent%d_%s_%s.root", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    f->SetRange(min, max);
    f->Write(Form("%s_cent%d_%s_%s.root", f->GetName(), icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    h->DrawCopy();
    f->DrawCopy("same");
    TLegend *l = gPad->BuildLegend(0.15, 0.12, 0.5, 0.3);
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetFillStyle(0);

    hyield_data->SetBinContent(icent + 1, integrated_data[kInt_YieldData]);
    hyield_data->SetBinError(icent + 1, integrated_data[kInt_YieldDataLin]);

    hyield_low->SetBinContent(icent + 1, integrated_data[kInt_YieldLow]);
    hyield_low->SetBinError(icent + 1, integrated_data[kInt_YieldLowErr]);

    hyield_high->SetBinContent(icent + 1, integrated_data[kInt_YieldHigh]);
    hyield_high->SetBinError(icent + 1, integrated_data[kInt_YieldHighErr]);

    hmean_data->SetBinContent(icent + 1, integrated_data[kInt_MeanData]);
    hmean_data->SetBinError(icent + 1, integrated_data[kInt_MeanDataLin]);

    hmean_low->SetBinContent(icent + 1, integrated_data[kInt_MeanLow]);
    hmean_low->SetBinError(icent + 1, integrated_data[kInt_MeanLowErr]);
    
    hmean_high->SetBinContent(icent + 1, integrated_data[kInt_MeanHigh]);
    hmean_high->SetBinError(icent + 1, integrated_data[kInt_MeanHighErr]);
    
    
  }

  hyield_stat->Add(hyield_data_stat);
  hyield_stat->Add(hyield_low_stat);
  hyield_stat->Add(hyield_high_stat);

  hyield_sys->Add(hyield_data_sys);
  hyield_sys->Add(hyield_low_sys);
  hyield_sys->Add(hyield_high_sys);

  hyield->Add(hyield_data);
  hyield->Add(hyield_low);
  hyield->Add(hyield_high);
  
  
  hmean_stat->Add(hmean_data_stat);
  hmean_stat->Add(hmean_low_stat);
  hmean_stat->Add(hmean_high_stat);
  //      hmean_stat->Divide(hyield_stat);
  hmean_sys->Add(hmean_data_sys);
  hmean_sys->Add(hmean_low_sys);
  hmean_sys->Add(hmean_high_sys);
  //      hmean_sys->Divide(hyield_sys);
  hmean->Add(hmean_data);
  hmean->Add(hmean_low);
  hmean->Add(hmean_high);
  
  fout->cd();

  hyield->Write();
  hmean->Write();

  hyield_stat->Write();
  hmean_stat->Write();

  hyield_sys->Write();
  hmean_sys->Write();

  hyield_data->Write();
  hmean_data->Write();

  hyield_data_stat->Write();
  hmean_data_stat->Write();

  hyield_data_sys->Write();
  hmean_data_sys->Write();

  hyield_low->Write();
  hmean_low->Write();

  hyield_high->Write();
  hmean_high->Write();

  hyield_low_stat->Write();
  hmean_low_stat->Write();

  hyield_high_stat->Write();
  hmean_high_stat->Write();

  hmean_low_sys->Write();
  hyield_low_sys->Write();
  
  hmean_high_sys->Write();
  hyield_high_sys->Write();

  cc->Write();
  
  fout->Close();
}

IntegratedProduction_pp(const Char_t *filename, Char_t *what = "", Int_t ifunc = 1, Option_t *opt = "0qI")
{
  TFile *filein = TFile::Open(filename);
  for (Int_t ipart = 0; ipart < 6; ipart++) {
    printf("*************************\n");
    printf("part=%d\n", ipart);
    TH1 *h = (TH1 *)filein->Get(Form("hComb_ITSsa%d_TPC%d_TOF%d_HMPID%d", ipart, ipart, ipart, ipart));
    IntegratedProduction(h, AliPID::ParticleMass(ipart % 3 + 2), ifunc, opt);
    if (gPad) gPad->SaveAs(Form("IntegratedProduction_pp_%d.root", ipart));
  }
}

IntegratedProduction(TH1 *h, TF1 *f = NULL, Option_t *opt = "0q", Float_t min = 0., Float_t max = 10., Double_t *integrated_data = NULL, Bool_t verbose = kFALSE)
{

  Double_t yield, yielderr, yielderrcorr, mean, meanerr, meanerrcorr, partyield[3], partyielderr[3], partyielderrcorr[3], partmean[3], partmeanerr[3], partmeanerrcorr[3];

  TVirtualFitter::SetMaxIterations(1000000);

  if (f) {
    Int_t fres = 1;
    Int_t fatt = 0;
    while (fres != 0 && fatt < 10) {
      fres = h->Fit(f, opt, "", min, max);
      fres = h->Fit(f, opt, "", min, max);
      printf("fit res = %d\n", fres);
      fatt++;
    }
  }

  // return;

  GetYieldAndMean(h, f, yield, yielderr, yielderrcorr, mean, meanerr, meanerrcorr, 0., 10., partyield, partyielderr, partyielderrcorr, partmean, partmeanerr, partmeanerrcorr);

  //  Double_t fint = f ? f->Integral(0.,10.) : 0.;
  //  Double_t finte = f ? f->IntegralError(0.,10.) : 0.;
  //  Double_t fmean = f ? f->Mean(0., 10.) : 0.;

  if (verbose) {
  printf("----\n");
  printf("dN/dy        = %f +- %f (%f)\n", yield, yielderr, yielderrcorr);
  printf("<pt>         = %f +- %f (%f)\n", mean, meanerr, meanerrcorr);
  printf("----\n");
  printf("dN/dy (data) = %f +- %f (%f)\n", partyield[0], partyielderr[0], partyielderrcorr[0]);
  printf("dN/dy (low)  = %f +- %f (%f)\n", partyield[1], partyielderr[1], partyielderrcorr[1]);
  printf("dN/dy (high) = %f +- %f (%f)\n", partyield[2], partyielderr[2], partyielderrcorr[2]);
  printf("<pt> (data)  = %f +- %f (%f)\n", partmean[0], partmeanerr[0], partmeanerrcorr[0]);
  printf("<pt> (low)   = %f +- %f (%f)\n", partmean[1], partmeanerr[1], partmeanerrcorr[1]);
  printf("<pt> (high)  = %f +- %f (%f)\n", partmean[2], partmeanerr[2], partmeanerrcorr[2]);
  printf("----\n");
  //  printf("dN/dy (func) = %f +- %f\n", fint, finte);
  //  printf("<pT> (func)  = %f +- %f\n", fmean, 0.);
  }

  if (!integrated_data) return;

  integrated_data[kInt_YieldData] = partyield[0];
  integrated_data[kInt_YieldDataErr] = partyielderr[0];
  integrated_data[kInt_YieldDataLin] = partyielderrcorr[0];
  integrated_data[kInt_MeanData] = partmean[0];
  integrated_data[kInt_MeanDataErr] = partmeanerr[0];
  integrated_data[kInt_MeanDataLin] = partmeanerrcorr[0];
  integrated_data[kInt_YieldLow] = partyield[1];
  integrated_data[kInt_YieldLowErr] = partyielderr[1];
  integrated_data[kInt_MeanLow] = partmean[1];
  integrated_data[kInt_MeanLowErr] = partmeanerr[1];
  integrated_data[kInt_YieldHigh] = partyield[2];
  integrated_data[kInt_YieldHighErr] = partyielderr[2];
  integrated_data[kInt_MeanHigh] = partmean[2];
  integrated_data[kInt_MeanHighErr] = partmeanerr[2];

  //  TH1 *hr = (TH1 *)h->Clone("hr");
  //  hr->Divide(f);
  //  new TCanvas;
  //  hr->Draw();

  //  TProfile *p = new TProfile("p", "", 100, 0., 10.);
  //  gROOT->LoadMacro("HistoUtils.C");
  //  HistoUtils_Function2Profile(f, p);
  //  p->Draw();
}

GetYieldAndMean(TH1 *h, TF1 *f, Double_t &yield, Double_t &yielderr, Double_t &yielderrcorr, Double_t &mean, Double_t &meanerr, Double_t &meanerrcorr, Double_t min, Double_t max, Double_t *partyield, Double_t *partyielderr, Double_t *partyielderrcorr, Double_t *partmean, Double_t *partmeanerr, Double_t *partmeanerrcorr)
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
  if (!f) {
    min = lo;
    max = hi;
    printf("WARNING: no function provided! Use only data with no extrapolation\n");
  }
  Double_t integral_lo = min < lo ? f->Integral(min, lo, (Double_t *)0, 1.e-6) : 0.;
  Double_t integralerr_lo = min < lo ? f->IntegralError(min, lo, 0, 0, 1.e-6) : 0.;
  Double_t meanintegral_lo = min < lo ? f->Mean(min, lo, (Double_t *)0, 1.e-6) * integral_lo : 0.;
  Double_t meanintegralerr_lo = min < lo ? f->Mean(min, lo, (Double_t *)0, 1.e-6) * integralerr_lo : 0.;
  
  /* integrate above the data */
  Double_t integral_hi = max > hi ? f->Integral(hi, max, (Double_t *)0, 1.e-6) : 0.;
  Double_t integralerr_hi = max > hi ? f->IntegralError(hi, max, 0, 0, 1.e-6) : 0.;
  Double_t meanintegral_hi = max > hi ? f->Mean(hi, max, (Double_t *)0, 1.e-6) * integral_hi : 0.;
  Double_t meanintegralerr_hi = max > hi ? f->Mean(hi, max, (Double_t *)0, 1.e-6) * integralerr_hi : 0.;

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

  /* set partial means */
  partmean[0] = meanintegral_data;
  partmeanerr[0] = meanintegralerr_data;
  partmeanerrcorr[0] = meanintegralerrcorr_data;
  partmean[1] = meanintegral_lo;
  partmeanerr[1] = meanintegralerr_lo;
  partmeanerrcorr[1] = meanintegralerr_lo;
  partmean[2] = meanintegral_hi;
  partmeanerr[2] = meanintegralerr_hi;
  partmeanerrcorr[2] = meanintegralerr_hi;
  
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
Convert_dNdy_1overpt_dNdy(TH1 *hin)
{

  TH1 *hout = hin->Clone("hout");
  hout->Reset();
  Double_t pt, mt, conv, val, vale;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    pt = hin->GetBinCenter(ibin + 1);
    val = hin->GetBinContent(ibin + 1);
    vale = hin->GetBinError(ibin + 1);
    val /= pt;
    vale /= pt;
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
SummedId_1over2pipt_dNdeta(const Char_t *filename, Int_t icent, Float_t etarange, Float_t scalePi = 1.)
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
      heta[ipart][icharge] = Convert_dNdy_1over2pipt_dNdeta(hy[ipart][icharge], AliPID::ParticleMass(ipart), etarange);
      if (ipart == 2) heta[ipart][icharge]->Scale(scalePi);
    }

  /* sum */
  TH1D *hsum = heta[2][0]->Clone("hsum");
  hsum->Reset();
  for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      hsum->Add(heta[ipart][icharge]);
    }
  for (Int_t ipt = 0; ipt < hsum->GetNbinsX(); ipt++) {
    Double_t err = 0.;
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < 2; icharge++) {
	err += heta[ipart][icharge]->GetBinError(ipt + 1);
      }
    hsum->SetBinError(ipt + 1, err);
  }
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

TH1 *
ReturnExtremeHighHisto(TH1 *hin, TH1 *herr = NULL)
{
  TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehigh", hin->GetName()));
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    if (herr) err = herr->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val + err);
  }
  return hout;
}

TH1 *
ReturnExtremeLowHisto(TH1 *hin, TH1 *herr = NULL)
{
  TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremelow", hin->GetName()));
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    if (herr) err = herr->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val - err);
  }
  return hout;
}

TH1 *
ReturnExtremeSoftHisto(TH1 *hin, TH1 *herr = NULL)
{
  return ReturnExtremeHisto(hin, herr, -1.);
}

TH1 *
ReturnExtremeHardHisto(TH1 *hin, TH1 *herr = NULL)
{
  return ReturnExtremeHisto(hin, herr, 1.);
}

TH1 *
ReturnExtremeHisto(TH1 *hin, TH1 *herr = NULL, Float_t sign = 1.)
{
  Double_t ptlow, pthigh;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    ptlow = hin->GetBinLowEdge(ibin + 1);
    break;
  }
  for (Int_t ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    pthigh = hin->GetBinLowEdge(ibin + 2);
    break;
  }

  Double_t mean = hin->GetMean();
  Double_t maxdiff = 0.;
  TH1 *hmax = NULL;
  for (Int_t inode = 0; inode < hin->GetNbinsX(); inode++) {

    Double_t ptnode = hin->GetBinCenter(inode + 1);
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehard", hin->GetName()));
    
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      Double_t val = hin->GetBinContent(ibin + 1);
      Double_t err = hin->GetBinError(ibin + 1);
      if (herr) err = herr->GetBinError(ibin + 1);
      Double_t cen = hin->GetBinCenter(ibin + 1);
      //    err *= -1. - (cen - ptlow) * (-1. - 1.) / (pthigh - ptlow);
      if (cen < ptnode)
        err *= -1. + (cen - ptlow) / (ptnode - ptlow);
      else
        err *= (cen - ptnode) / (pthigh - ptnode);

      hout->SetBinContent(ibin + 1, val + sign * err);
    }

    Double_t diff = TMath::Abs(mean - hout->GetMean());
    if (diff > maxdiff) {
      //      printf("found max at %f\n", ptnode);
      if (hmax) delete hmax;
      hmax = (TH1 *)hout->Clone("hmax");
      maxdiff = diff;
    }
    delete hout;
  }
  return hmax;
}

TGraphErrors *
ReturnExtremeSoftGraph(TGraphErrors *hin, TGraphErrors *herr = NULL)
{
  TGraphErrors *hout = (TGraphErrors *)hin->Clone(Form("%s_extremesoft", hin->GetName()));
  Double_t ptlow, pthigh;
  Double_t ptlow = hin->GetX()[0];
  Double_t pthigh = hin->GetX()[hin->GetN() - 1];

  for (Int_t ibin = 0; ibin < hin->GetN(); ibin++) {
    Double_t val = hin->GetY()[ibin];
    Double_t err = hin->GetEY()[ibin];
    if (herr) err = herr->GetEY()[ibin];
    Double_t cen = hin->GetX()[ibin];
    err *= 1. + (cen - ptlow) * (-1. - 1.) / (pthigh - ptlow);
    hout->SetPoint(ibin, cen, val + err);
  }
  return hout;
}


TGraphErrors *
ReturnExtremeHardGraph(TGraphErrors *hin, TGraphErrors *herr = NULL)
{
  TGraphErrors *hout = (TGraphErrors *)hin->Clone(Form("%s_extremehard", hin->GetName()));
  Double_t ptlow, pthigh;
  Double_t ptlow = hin->GetX()[0];
  Double_t pthigh = hin->GetX()[hin->GetN() - 1];

  for (Int_t ibin = 0; ibin < hin->GetN(); ibin++) {
    Double_t val = hin->GetY()[ibin];
    Double_t err = hin->GetEY()[ibin];
    if (herr) err = herr->GetEY()[ibin];
    Double_t cen = hin->GetX()[ibin];
    err *= -1. - (cen - ptlow) * (-1. - 1.) / (pthigh - ptlow);
    hout->SetPoint(ibin, cen, val + err);
  }
  return hout;
}

