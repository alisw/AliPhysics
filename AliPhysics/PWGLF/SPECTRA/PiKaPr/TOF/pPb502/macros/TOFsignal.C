#if 0
TF1 *
TOFsignal(Float_t norm, Float_t mean, Float_t sigma, Float_t tail)
{

  TF1 *f = new TF1("fTOFsignal", "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * TMath::Exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))");

  f->SetParameter(0, norm);
  f->SetParameter(1, mean);
  f->SetParameter(2, sigma);
  f->SetParameter(3, tail);
  f->SetRange(mean - 10. * sigma, mean + 10. * sigma);
  return f;
}
#endif

Double_t
TOFsignal(Double_t *x, Double_t *par)
{
  Double_t norm = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t tail = par[3];
  
  if (x[0] <= (tail + mean))
    return norm * TMath::Gaus(x[0], mean, sigma);
  else
    return norm * TMath::Gaus(tail + mean, mean, sigma) * TMath::Exp(-tail * (x[0] - tail - mean) / (sigma * sigma));
}

Double_t
TOFsignal_double(Double_t *x, Double_t *par)
{
  return TOFsignal(x, par) + TOFsignal(x, &par[4]);
}

Double_t
TOFsignal_triple(Double_t *x, Double_t *par)
{
  return TOFsignal(x, par) + TOFsignal(x, &par[4]) + TOFsignal(x, &par[8]);
}

static TF1 *fIntegrand = NULL;
Double_t 
TOFsignal_Integrand(const Double_t *x, const Double_t *par)
{
  Double_t f = TOFsignal(x, par);
  Double_t g = /*TMath::Abs(par[5] - x[0]) < 1.e-1 ? 1. : 0.;*/TMath::Gaus(x[0], par[5], par[4], kTRUE);
  return f * g;
}

Double_t
TOFsignal_convolution(const Double_t *x, const Double_t *par)
{
  if (!fIntegrand)
    fIntegrand = new TF1("fIntegrand", TOFsignal_Integrand, -10., 10., 6);
  fIntegrand->SetParameters(par[0], par[1], par[2], par[3], par[4], x[0]);
  Double_t integral = fIntegrand->Integral(-5. * par[4], 5. * par[4]);
  return integral;
}
