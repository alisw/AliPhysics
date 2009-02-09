////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDirectYlm - Correlation function that is binned in Ylms    //
// directly. Provides a way to store the numerator and denominator            //
// in Ylms directly and correctly calculate the correlation                   //
// function from them.                                                        //
//                                                                            //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCorrFctnDirectYlm.h"
#include <TMath.h>
#include <iostream>

using namespace std;

AliFemtoCorrFctnDirectYlm::AliFemtoCorrFctnDirectYlm(const char *name, int maxl, int ibin=30, double vmin=0.0, double vmax=0.3):
  fnumsreal(0),
  fnumsimag(0),
  fdensreal(0),
  fdensimag(0),
  fbinctn(0),
  fbinctd(0),
  fcovnum(0),
  fcovden(0),
  fcovmnum(0),
  fcovmden(0),
  fMaxL(0),
  fMaxJM(0),
  fels(0),
  fems(0),
  felsi(0),
  femsi(0),
  fYlmBuffer(0),
  factorials(0),
  fSout(0.0),
  fSside(0.0),
  fSlong(0.0)
{
  // Main constructor
  fMaxL = maxl;
  fMaxJM = (maxl+1)*(maxl+1);

  cout <<  "Size is " << sizeof(double) << " " << sizeof(complex<double>) << endl;

  // Fill in factorials table
  factorials = (double *) malloc(sizeof(double) * (4 * (maxl + 1)));
  int fac = 1;
  factorials[0] = 1;
  for (int iter=1; iter<4*(maxl+1); iter++)
    {
      fac *= iter;
      factorials[iter] = fac;
    }

  // Fill in els and ems table
  int el = 0;
  int em = 0;
  int il = 0;
  fels = (double *) malloc(sizeof(double) * (fMaxJM));
  fems = (double *) malloc(sizeof(double) * (fMaxJM));
  felsi = (int *) malloc(sizeof(int) * (fMaxJM));
  femsi = (int *) malloc(sizeof(int) * (fMaxJM));
  do {
    fels[il] = el;
    fems[il] = em;
    felsi[il] = (int) el;
    femsi[il] = (int) em;

    cout << "il el em " << il << " " << felsi[il] << " " << femsi[il] << endl;
    em++;
    il++;
    if (em > el) {
      el++;
      em = -el;
    }
  }
  while (el <= maxl);
  
  for (il=0; il<fMaxJM; il++)
    cout << "il el em " << il << " " << felsi[il] << " " << femsi[il] << endl;

  // Create numerator and denominator historgrams
  //  int sthp = sizeof(TH1D *);
  //  fnumsreal = (TH1D **) malloc(sthp * fMaxJM);
//   fnumsreal = new TH1D * [fMaxJM];
//   fnumsimag = new TH1D * [fMaxJM];
//   fdensreal = new TH1D * [fMaxJM];
//   fdensimag = new TH1D * [fMaxJM];
  fnumsreal = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fnumsimag = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fdensreal = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fdensimag = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  
  char bufname[200];
  for (int ihist=0; ihist<fMaxJM; ihist++) {
    sprintf(bufname, "NumReYlm%i%i%s", felsi[ihist], femsi[ihist]<0 ? felsi[ihist]-femsi[ihist] : femsi[ihist], name);
    fnumsreal[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);
    sprintf(bufname, "NumImYlm%i%i%s", felsi[ihist], femsi[ihist]<0 ? felsi[ihist]-femsi[ihist] : femsi[ihist], name);
    fnumsimag[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);
    sprintf(bufname, "DenReYlm%i%i%s", felsi[ihist], femsi[ihist]<0 ? felsi[ihist]-femsi[ihist] : femsi[ihist], name);
    fdensreal[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);
    sprintf(bufname, "DenImYlm%i%i%s", felsi[ihist], femsi[ihist]<0 ? felsi[ihist]-femsi[ihist] : femsi[ihist], name);
    fdensimag[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);

    fnumsreal[ihist]->Sumw2();
    fnumsimag[ihist]->Sumw2();
    fdensreal[ihist]->Sumw2();
    fdensimag[ihist]->Sumw2();
  }

  sprintf(bufname, "BinCountNum%s", name);
  fbinctn = new TH1D(bufname, bufname, ibin, vmin, vmax);

  sprintf(bufname, "BinCountDen%s", name);
  fbinctd = new TH1D(bufname, bufname, ibin, vmin, vmax);

  fYlmBuffer = (complex<double> *) malloc(sizeof(complex<double>) * fMaxJM);
  
  // Covariance matrices
  fcovmnum = (double *) malloc(sizeof(double) * fMaxJM * fMaxJM * 4 * ibin);
  fcovmden = (double *) malloc(sizeof(double) * fMaxJM * fMaxJM * 4 * ibin);

  fcovnum = 0;
  fcovden = 0;

  AliFemtoYlm::InitializeYlms();
}


AliFemtoCorrFctnDirectYlm::AliFemtoCorrFctnDirectYlm():
  fnumsreal(0),
  fnumsimag(0),
  fdensreal(0),
  fdensimag(0),
  fbinctn(0),
  fbinctd(0),
  fcovnum(0),
  fcovden(0),
  fcovmnum(0),
  fcovmden(0),
  fMaxL(0),
  fMaxJM(0),
  fels(0),
  fems(0),
  felsi(0),
  femsi(0),
  fYlmBuffer(0),
  factorials(0),
  fSout(0.0),
  fSside(0.0),
  fSlong(0.0)
{
  // Default constructor
  AliFemtoCorrFctnDirectYlm("AliFemtoCorrFctnDirectYlm",2);
}

AliFemtoCorrFctnDirectYlm::AliFemtoCorrFctnDirectYlm(const AliFemtoCorrFctnDirectYlm& aCorrFctn):
  AliFemtoCorrFctn(),
  fnumsreal(0),
  fnumsimag(0),
  fdensreal(0),
  fdensimag(0),
  fbinctn(0),
  fbinctd(0),
  fcovnum(0),
  fcovden(0),
  fcovmnum(0),
  fcovmden(0),
  fMaxL(0),
  fMaxJM(0),
  fels(0),
  fems(0),
  felsi(0),
  femsi(0),
  fYlmBuffer(0),
  factorials(0),
  fSout(0.0),
  fSside(0.0),
  fSlong(0.0)
{
  // Copy constructor
  int ibin = aCorrFctn.fbinctn->GetNbinsX();

  fMaxL = aCorrFctn.fMaxL;
  fMaxJM = (fMaxL+1)*(fMaxL+1);

  // Fill in factorials table
  factorials = (double *) malloc(sizeof(double) * (4 * (fMaxL + 1)));
  for (int iter=1; iter<4*(fMaxL+1); iter++)
    {
      factorials[iter] = aCorrFctn.factorials[iter];
    }

  // Fill in els and ems table
  int el = 0;
  int em = 0;
  int il = 0;
  fels = (double *) malloc(sizeof(double) * (fMaxJM));
  fems = (double *) malloc(sizeof(double) * (fMaxJM));
  felsi = (int *) malloc(sizeof(int) * (fMaxJM));
  femsi = (int *) malloc(sizeof(int) * (fMaxJM));
  do {
    fels[il] = el;
    fems[il] = em;
    felsi[il] = (int) el;
    femsi[il] = (int) em;

    em++;
    il++;
    if (em > el) {
      el++;
      em = -el;
    }
  }
  while (el <= fMaxL);
  
  fnumsreal = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fnumsimag = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fdensreal = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fdensimag = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  
  for (int ihist=0; ihist<fMaxJM; ihist++) {
    if (aCorrFctn.fnumsreal[ihist])
      fnumsreal[ihist] = new TH1D(*aCorrFctn.fnumsreal[ihist]);
    else
      fnumsreal[ihist] = 0;
    if (aCorrFctn.fnumsimag[ihist])
      fnumsimag[ihist] = new TH1D(*aCorrFctn.fnumsimag[ihist]);
    else
      fnumsimag[ihist] = 0;
    if (aCorrFctn.fdensreal[ihist])
      fdensreal[ihist] = new TH1D(*aCorrFctn.fdensreal[ihist]);
    else
      fdensreal[ihist] = 0;
    if (aCorrFctn.fdensimag[ihist])
      fdensimag[ihist] = new TH1D(*aCorrFctn.fdensimag[ihist]);
    else
      fdensimag[ihist] = 0;
  }

  if (aCorrFctn.fbinctn) 
    fbinctn = new TH1D(*aCorrFctn.fbinctn);
  else
    fbinctn = 0;
  if (aCorrFctn.fbinctd) 
    fbinctd = new TH1D(*aCorrFctn.fbinctd);
  else
    fbinctd = 0;

  fYlmBuffer = (complex<double> *) malloc(sizeof(complex<double>) * fMaxJM);
  
  // Covariance matrices
  fcovmnum = (double *) malloc(sizeof(double) * fMaxJM * fMaxJM * 4 * ibin);
  fcovmden = (double *) malloc(sizeof(double) * fMaxJM * fMaxJM * 4 * ibin);

  for (int iter=0; iter<fMaxJM * fMaxJM * 4 * ibin; iter++) {
    fcovmnum[iter] = aCorrFctn.fcovmnum[iter];
    fcovmden[iter] = aCorrFctn.fcovmden[iter];
  }

  if (aCorrFctn.fcovnum)
    fcovnum = new TH3D(*aCorrFctn.fcovnum);
  else
    fcovnum = 0;
  if (aCorrFctn.fcovden)
    fcovden = new TH3D(*aCorrFctn.fcovden);
  else
    fcovden = 0;

  fSout = aCorrFctn.fSout;
  fSside = aCorrFctn.fSside;
  fSlong = aCorrFctn.fSlong;
}

AliFemtoCorrFctnDirectYlm& AliFemtoCorrFctnDirectYlm::operator=(const AliFemtoCorrFctnDirectYlm& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;
  
  int ibin = aCorrFctn.fbinctn->GetNbinsX();

  fMaxL = aCorrFctn.fMaxL;
  fMaxJM = (fMaxL+1)*(fMaxL+1);

  // Fill in factorials table
  factorials = (double *) malloc(sizeof(double) * (4 * (fMaxL + 1)));
  for (int iter=1; iter<4*(fMaxL+1); iter++)
    {
      factorials[iter] = aCorrFctn.factorials[iter];
    }

  // Fill in els and ems table
  int el = 0;
  int em = 0;
  int il = 0;
  fels = (double *) malloc(sizeof(double) * (fMaxJM));
  fems = (double *) malloc(sizeof(double) * (fMaxJM));
  felsi = (int *) malloc(sizeof(int) * (fMaxJM));
  femsi = (int *) malloc(sizeof(int) * (fMaxJM));
  do {
    fels[il] = el;
    fems[il] = em;
    felsi[il] = (int) el;
    femsi[il] = (int) em;

    em++;
    il++;
    if (em > el) {
      el++;
      em = -el;
    }
  }
  while (el <= fMaxL);
  
  fnumsreal = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fnumsimag = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fdensreal = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  fdensimag = (TH1D **) malloc(sizeof(TH1D *) * fMaxJM);
  
  for (int ihist=0; ihist<fMaxJM; ihist++) {
    if (aCorrFctn.fnumsreal[ihist])
      fnumsreal[ihist] = new TH1D(*aCorrFctn.fnumsreal[ihist]);
    else
      fnumsreal[ihist] = 0;
    if (aCorrFctn.fnumsimag[ihist])
      fnumsimag[ihist] = new TH1D(*aCorrFctn.fnumsimag[ihist]);
    else
      fnumsimag[ihist] = 0;
    if (aCorrFctn.fdensreal[ihist])
      fdensreal[ihist] = new TH1D(*aCorrFctn.fdensreal[ihist]);
    else
      fdensreal[ihist] = 0;
    if (aCorrFctn.fdensimag[ihist])
      fdensimag[ihist] = new TH1D(*aCorrFctn.fdensimag[ihist]);
    else
      fdensimag[ihist] = 0;
  }

  if (aCorrFctn.fbinctn) 
    fbinctn = new TH1D(*aCorrFctn.fbinctn);
  else
    fbinctn = 0;
  if (aCorrFctn.fbinctd) 
    fbinctd = new TH1D(*aCorrFctn.fbinctd);
  else
    fbinctd = 0;

  fYlmBuffer = (complex<double> *) malloc(sizeof(complex<double>) * fMaxJM);
  
  // Covariance matrices
  fcovmnum = (double *) malloc(sizeof(double) * fMaxJM * fMaxJM * 4 * ibin);
  fcovmden = (double *) malloc(sizeof(double) * fMaxJM * fMaxJM * 4 * ibin);

  for (int iter=0; iter<fMaxJM * fMaxJM * 4 * ibin; iter++) {
    fcovmnum[iter] = aCorrFctn.fcovmnum[iter];
    fcovmden[iter] = aCorrFctn.fcovmden[iter];
  }

  if (aCorrFctn.fcovnum)
    fcovnum = new TH3D(*aCorrFctn.fcovnum);
  else
    fcovnum = 0;
  if (aCorrFctn.fcovden)
    fcovden = new TH3D(*aCorrFctn.fcovden);
  else
    fcovden = 0;

  fSout = aCorrFctn.fSout;
  fSside = aCorrFctn.fSside;
  fSlong = aCorrFctn.fSlong;

  return *this;
}

AliFemtoCorrFctnDirectYlm::~AliFemtoCorrFctnDirectYlm()
{
  // Destructor
  for (int ihist=0; ihist<fMaxJM; ihist++) {
    delete fnumsreal[ihist];
    delete fnumsimag[ihist];
    delete fdensreal[ihist];
    delete fdensimag[ihist];
  }

  delete fbinctn;
  delete fbinctd;

  //  delete fnumsreal;
  //  delete fnumsimag;
  //  delete fdensreal;
  //  delete fdensimag;

  free( fnumsreal);
  free( fnumsimag);
  free( fdensreal);
  free( fdensimag);

  free(factorials);
  free(fels);
  free(fems);
  free(felsi);
  free(femsi);
  free(fYlmBuffer);

  free(fcovmnum);
  free(fcovmden);

  if (fcovnum) delete fcovnum;
  if (fcovden) delete fcovden;
}

double AliFemtoCorrFctnDirectYlm::ClebschGordan(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm)
{
  // Calculate Clebsh-Gordan coefficient
  int mint, maxt;
  double cgc = 0.0;
  int titer;
  double coef;

  maxt = lrint(aJot1 + aJot2 - aJot);
  mint = 0;
  if (lrint(aJot1 - aEm1) < maxt) maxt = lrint(aJot1 - aEm1);
  if (lrint(aJot2 + aEm2) < maxt) maxt = lrint(aJot2 + aEm2);
  if (lrint(-(aJot-aJot2+aEm1)) > mint) mint = lrint(-(aJot-aJot2+aEm1));
  if (lrint(-(aJot-aJot1-aEm2)) > mint) mint = lrint(-(aJot-aJot1-aEm2));

  for (titer = mint; titer<=maxt; titer ++)
    {
      coef = TMath::Power(-1, titer);
      coef *= TMath::Sqrt((2*aJot+1)*
			  factorials[lrint(aJot1+aEm1)] *
			  factorials[lrint(aJot1-aEm1)] *
			  factorials[lrint(aJot2+aEm2)] *
			  factorials[lrint(aJot2-aEm2)] *
			  factorials[lrint(aJot+aEm)] *
			  factorials[lrint(aJot-aEm)]);
      coef /= (factorials[titer] *
	       factorials[lrint(aJot1+aJot2-aJot-titer)] *
	       factorials[lrint(aJot1-aEm1-titer)] *
	       factorials[lrint(aJot2+aEm2-titer)] *
	       factorials[lrint(aJot-aJot2+aEm1+titer)] *
	       factorials[lrint(aJot-aJot1-aEm2+titer)]);
      
      cgc += coef;
    }

  cgc *= DeltaJ(aJot1, aJot2, aJot);

  return cgc;
}

double AliFemtoCorrFctnDirectYlm::DeltaJ(double aJot1, double aJot2, double aJot)
{
  // Calculate J for the Clebsh-Gordan coefficient
  if ((aJot1+aJot2-aJot) < 0) {
    //    cout << "J1+J2-J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((aJot1-aJot2+aJot) < 0) {
    //    cout << "J1-J2+J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((-aJot1+aJot2+aJot) < 0) {
    //    cout << "-J1+J2+J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((aJot1+aJot2+aJot+1) < 0) {
    //    cout << "J1+J2+J3+1 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  double res = TMath::Sqrt(1.0 * 
			   factorials[lrint(aJot1+aJot2-aJot)] * 
			   factorials[lrint(aJot1-aJot2+aJot)] * 
			   factorials[lrint(-aJot1+aJot2+aJot)] / 
			   factorials[lrint(aJot1+aJot2+aJot+1)]);
  
  return res;
}

double AliFemtoCorrFctnDirectYlm::WignerSymbol(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm)
{
  // Get Wigner symbol
  if (lrint(aEm1+aEm2+aEm) != 0.0) 
    return 0.0;
  double cge = ClebschGordan(aJot1, aEm1, aJot2, aEm2, aJot, -aEm);
  if (lrint(abs(aJot1 - aJot2 - aEm)) % 2) 
    cge *= -1.0;
  cge /= sqrt(2*aJot + 1);

  if (cge == -0.0) cge = 0.0;

  return cge;
}


void AliFemtoCorrFctnDirectYlm::GetMtilde(complex<double> *aMat, double *aMTilde)
{
  // Create the Mtilde for a given q bin
  double lzero, mzero;
  double lprim, mprim;
  double lbis, mbis;
 
  int lzeroi, mzeroi;
  int lprimi, mprimi;
  int lbisi, mbisi;

  complex<double> mcomp;

  for (int izero = 0; izero<GetMaxJM(); izero++) {
    GetElEmForIndex(izero, &lzero, &mzero);
    GetElEmForIndex(izero, &lzeroi, &mzeroi);
    for (int ibis = 0; ibis<GetMaxJM(); ibis++) {
      GetElEmForIndex(ibis, &lbis, &mbis);
      GetElEmForIndex(ibis, &lbisi, &mbisi);
      complex<double> val = complex<double>(0.0, 0.0);
      for (int iprim = 0; iprim<GetMaxJM(); iprim++) {

	GetElEmForIndex(iprim, &lprim, &mprim);
	GetElEmForIndex(iprim, &lprimi, &mprimi);

	if (abs(mzeroi) % 2) mcomp = complex<double>(-1.0, 0.0); // (-1)^m
	else mcomp = complex<double>(1.0, 0.0);
	
	mcomp *= sqrt((2*lzero+1)*(2*lprim+1)*(2*lbis+1));   // P1
	mcomp *= WignerSymbol(lzero, 0, lprim, 0, lbis, 0); // W1
	mcomp *= WignerSymbol(lzero, -mzero, lprim, mprim, lbis, mbis); // W2
	mcomp *= aMat[iprim];
	val += mcomp;
      }
      aMTilde[(izero*2)*(2*GetMaxJM()) + (ibis*2)]     =  real(val);
      aMTilde[(izero*2+1)*(2*GetMaxJM()) + (ibis*2)]   =  imag(val);
      if (imag(val) != 0.0)
	aMTilde[(izero*2)*(2*GetMaxJM()) + (ibis*2+1)]   = -imag(val);
      else 
	aMTilde[(izero*2)*(2*GetMaxJM()) + (ibis*2+1)]   = 0.0;
      aMTilde[(izero*2+1)*(2*GetMaxJM()) + (ibis*2+1)] =  real(val);
      
    }
  }
}

int  AliFemtoCorrFctnDirectYlm::GetMaxJM() const
{ return fMaxJM; }

void AliFemtoCorrFctnDirectYlm::GetElEmForIndex(int aIndex, double *aEl, double *aEm) const
{
  // Get l,m for a given index
  *aEl = fels[aIndex];
  *aEm = fems[aIndex];
}

void AliFemtoCorrFctnDirectYlm::GetElEmForIndex(int aIndex, int *aEl, int *aEm) const
{
  // Get l,m for a given index
  *aEl = felsi[aIndex];
  *aEm = femsi[aIndex];
}

int AliFemtoCorrFctnDirectYlm::GetBin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag)
{
  return (qbin*GetMaxJM()*GetMaxJM()*4 +
	  (ilmprim*2 + primimag) * GetMaxJM()*2 +
	  ilmzero*2 + zeroimag);
}

void AliFemtoCorrFctnDirectYlm::AddRealPair(double qout, double qside, double qlong, double weight)
{
  // Fill numerator
  double kv = sqrt(qout*qout + qside*qside + qlong*qlong);
  int nqbin = fbinctn->GetXaxis()->FindFixBin(kv) - 1;
  
  // Use saved ylm values for same qout, qside, qlong
  if ((qout != fSout) || (qside != fSside) || (qlong != fSlong)) {
    AliFemtoYlm::YlmUpToL(fMaxL, qout, qside, qlong, fYlmBuffer);
    fSout = qout; fSside = qside; fSlong = qlong;
  }
  for (int ilm=0; ilm<GetMaxJM(); ilm++) {
    //    fYlmBuffer[ilm] = AliFemtoYlm::Ylm(elsi[ilm], emsi[ilm], qout, qside, qlong);

    fnumsreal[ilm]->Fill(kv, real(fYlmBuffer[ilm])*weight);
    fnumsimag[ilm]->Fill(kv, -imag(fYlmBuffer[ilm])*weight);

    fbinctn->Fill(kv, 1.0);
  }

  // Fill in the error matrix
  //  int tabshift = nqbin*GetMaxJM()*GetMaxJM()*4;
  if (nqbin < fbinctn->GetNbinsX())
    for (int ilmzero=0; ilmzero<GetMaxJM(); ilmzero++)
      for (int ilmprim=0; ilmprim<GetMaxJM(); ilmprim++) {
	fcovmnum[GetBin(nqbin, ilmzero, 0, ilmprim, 0)] += real(fYlmBuffer[ilmzero])*real(fYlmBuffer[ilmprim])*weight*weight;
	fcovmnum[GetBin(nqbin, ilmzero, 0, ilmprim, 1)] += real(fYlmBuffer[ilmzero])*-imag(fYlmBuffer[ilmprim])*weight*weight;
	fcovmnum[GetBin(nqbin, ilmzero, 1, ilmprim, 0)] += -imag(fYlmBuffer[ilmzero])*real(fYlmBuffer[ilmprim])*weight*weight;
	fcovmnum[GetBin(nqbin, ilmzero, 1, ilmprim, 1)] += -imag(fYlmBuffer[ilmzero])*-imag(fYlmBuffer[ilmprim])*weight*weight;
	
      }
  
}

void AliFemtoCorrFctnDirectYlm::AddMixedPair(double qout, double qside, double qlong, double weight)
{
  // Fill denominator
  double kv = sqrt(qout*qout + qside*qside + qlong*qlong);
  
  // Use saved ylm values for same qout, qside, qlong
  if ((qout != fSout) || (qside != fSside) || (qlong != fSlong)) {
    AliFemtoYlm::YlmUpToL(fMaxL, qout, qside, qlong, fYlmBuffer);
    fSout = qout; fSside = qside; fSlong = qlong;
  }
  for (int ilm=0; ilm<GetMaxJM(); ilm++) {
    //    fYlmBuffer[ilm] = AliFemtoYlm::Ylm(elsi[ilm], emsi[ilm], qout, qside, qlong);

    fdensreal[ilm]->Fill(kv, real(fYlmBuffer[ilm])*weight);
    fdensimag[ilm]->Fill(kv, -imag(fYlmBuffer[ilm])*weight);

    fbinctd->Fill(kv, 1.0);
  }

  // Fill in the error matrix
  int nqbin = fbinctn->GetXaxis()->FindFixBin(kv) - 1;
  //  int tabshift = nqbin*GetMaxJM()*GetMaxJM()*4;
  if (nqbin < fbinctn->GetNbinsX())
    for (int ilmzero=0; ilmzero<GetMaxJM(); ilmzero++)
      for (int ilmprim=0; ilmprim<GetMaxJM(); ilmprim++) {
	fcovmden[GetBin(nqbin, ilmzero, 0, ilmprim, 0)] += real(fYlmBuffer[ilmzero])*real(fYlmBuffer[ilmprim]);
	fcovmden[GetBin(nqbin, ilmzero, 0, ilmprim, 1)] += real(fYlmBuffer[ilmzero])*-imag(fYlmBuffer[ilmprim]);
	fcovmden[GetBin(nqbin, ilmzero, 1, ilmprim, 0)] += -imag(fYlmBuffer[ilmzero])*real(fYlmBuffer[ilmprim]);
	fcovmden[GetBin(nqbin, ilmzero, 1, ilmprim, 1)] += -imag(fYlmBuffer[ilmzero])*-imag(fYlmBuffer[ilmprim]);
	
    }
}

void AliFemtoCorrFctnDirectYlm::AddRealPair(double *qvec, double weight) {
  AddRealPair(qvec[0], qvec[1], qvec[2], weight);
}

void AliFemtoCorrFctnDirectYlm::AddMixedPair(double *qvec, double weight) {
  AddMixedPair(qvec[0], qvec[1], qvec[2], weight);
}

void AliFemtoCorrFctnDirectYlm::Finish()
{
  PackCovariances();
}

void AliFemtoCorrFctnDirectYlm::Write()
{
  // Write out output histograms
  for (int ilm=0; ilm<fMaxJM; ilm++) {
    fnumsreal[ilm]->Write();
    fdensreal[ilm]->Write();
    fnumsimag[ilm]->Write();
    fdensimag[ilm]->Write();
  }
  if (fcovnum) fcovnum->Write();
  if (fcovden) fcovden->Write();
}

TList* AliFemtoCorrFctnDirectYlm::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  for (int ilm=0; ilm<fMaxJM; ilm++) {
    tOutputList->Add(fnumsreal[ilm]);
    tOutputList->Add(fdensreal[ilm]);
    tOutputList->Add(fnumsimag[ilm]);
    tOutputList->Add(fdensimag[ilm]);
  }
  if (fcovnum) tOutputList->Add(fcovnum);
  if (fcovden) tOutputList->Add(fcovden);

  return tOutputList;
}


void AliFemtoCorrFctnDirectYlm::ReadFromFile(TFile *infile, const char *name, int maxl)
{
  // Raad in the numerator and denominator from file
  if (maxl != fMaxL) {
    cout << "Cannot read function for L " << maxl << " into a container with L "<< fMaxL << endl;
    return;
  }
  cout << "Reading in numerators and denominators" << endl;

  char bufname[200];
  for (int ihist=0; ihist<fMaxJM; ihist++) {
    sprintf(bufname, "NumReYlm%i%i%s", felsi[ihist], femsi[ihist]<0 ? felsi[ihist]-femsi[ihist] : femsi[ihist], name);
    if (fnumsreal[ihist]) delete fnumsreal[ihist];
    fnumsreal[ihist] = new TH1D(*((TH1D *) infile->Get(bufname)));

    sprintf(bufname, "NumImYlm%i%i%s", felsi[ihist], femsi[ihist]<0 ? felsi[ihist]-femsi[ihist] : femsi[ihist], name);
    if (fnumsimag[ihist]) delete fnumsimag[ihist];
    fnumsimag[ihist] = new TH1D(*((TH1D *) infile->Get(bufname)));

    sprintf(bufname, "DenReYlm%i%i%s", felsi[ihist], femsi[ihist]<0 ? felsi[ihist]-femsi[ihist] : femsi[ihist], name);
    if (fdensreal[ihist]) delete fdensreal[ihist];
    fdensreal[ihist] = new TH1D(*((TH1D *) infile->Get(bufname)));

    sprintf(bufname, "DenImYlm%i%i%s", felsi[ihist], femsi[ihist]<0 ? felsi[ihist]-femsi[ihist] : femsi[ihist], name);
    if (fdensimag[ihist]) delete fdensimag[ihist];
    fdensimag[ihist] = new TH1D(*((TH1D *) infile->Get(bufname)));
  }

  if (fcovnum) delete fcovnum;
  sprintf(bufname, "covNum%s", name);
  fcovnum = new TH3D (*((TH3D *) infile->Get(bufname)));

  if (fcovden) delete fcovden;
  sprintf(bufname, "CovDen%s", name);
  fcovden = new TH3D (*((TH3D *) infile->Get(bufname)));

  if ((fcovnum) && (fcovden)) {
    cout << "Unpacking covariance matrices from file " << endl;
    UnpackCovariances();
  }
  else {

    cout << "Creating fake covariance matrices" << endl;
  
    for (int ibin=1; ibin<=fnumsreal[0]->GetNbinsX(); ibin++) {
      double nent = fnumsreal[0]->GetEntries();
      double nentd = fdensreal[0]->GetEntries();
      for (int ilmx=0; ilmx<GetMaxJM(); ilmx++) {
	for (int ilmy=0; ilmy<GetMaxJM(); ilmy++) {
	  double t1t2rr = fnumsreal[ilmx]->GetBinContent(ibin)*fnumsreal[ilmy]->GetBinContent(ibin)/nent/nent;
	  double t1t2ri = fnumsreal[ilmx]->GetBinContent(ibin)*fnumsimag[ilmy]->GetBinContent(ibin)/nent/nent;
	  double t1t2ir = fnumsimag[ilmx]->GetBinContent(ibin)*fnumsreal[ilmy]->GetBinContent(ibin)/nent/nent;
	  double t1t2ii = fnumsimag[ilmx]->GetBinContent(ibin)*fnumsimag[ilmy]->GetBinContent(ibin)/nent/nent;
	  if (ilmx == ilmy) {
	    fcovmnum[GetBin(ibin-1, ilmx, 0, ilmy, 0)] = nent*(TMath::Power(fnumsreal[ilmx]->GetBinError(ibin)/nent,2)*(nent-1) + t1t2rr);
	    fcovmnum[GetBin(ibin-1, ilmx, 0, ilmy, 1)] = nent*t1t2ri;
	    fcovmnum[GetBin(ibin-1, ilmx, 1, ilmy, 0)] = nent*t1t2ir;
	    fcovmnum[GetBin(ibin-1, ilmx, 1, ilmy, 1)] = nent*(TMath::Power(fnumsimag[ilmx]->GetBinError(ibin)/nent,2)*(nent-1) + t1t2rr);
	  }
	  else {
	    fcovmnum[GetBin(ibin-1, ilmx, 0, ilmy, 0)] = nent*t1t2rr;
	    fcovmnum[GetBin(ibin-1, ilmx, 0, ilmy, 1)] = nent*t1t2ri;
	    fcovmnum[GetBin(ibin-1, ilmx, 1, ilmy, 0)] = nent*t1t2ir;
	    fcovmnum[GetBin(ibin-1, ilmx, 1, ilmy, 1)] = nent*t1t2ii;
	  }
	  t1t2rr = fdensreal[ilmx]->GetBinContent(ibin)*fdensreal[ilmy]->GetBinContent(ibin)/nentd/nentd;
	  t1t2ri = fdensreal[ilmx]->GetBinContent(ibin)*fdensimag[ilmy]->GetBinContent(ibin)/nentd/nentd;
	  t1t2ir = fdensimag[ilmx]->GetBinContent(ibin)*fdensreal[ilmy]->GetBinContent(ibin)/nentd/nentd;
	  t1t2ii = fdensimag[ilmx]->GetBinContent(ibin)*fdensimag[ilmy]->GetBinContent(ibin)/nentd/nentd;
	  
	  fcovmden[GetBin(ibin-1, ilmx, 0, ilmy, 0)] = nentd*t1t2rr;
	  fcovmden[GetBin(ibin-1, ilmx, 0, ilmy, 1)] = nentd*t1t2ri;
	  fcovmden[GetBin(ibin-1, ilmx, 1, ilmy, 0)] = nentd*t1t2ir;
	  fcovmden[GetBin(ibin-1, ilmx, 1, ilmy, 1)] = nentd*t1t2ii;
	}
      }
    }
  }

  // Recalculating the correlation functions
  Finish();
}

int AliFemtoCorrFctnDirectYlm::PackYlmVector(const double *invec, double *outvec)
{
  // Pack a vector in l,m into an array using
  // only independent components
  int ioutcount = 0;
  int em, el;
  for (int ilm=0; ilm<GetMaxJM(); ilm++) {
    GetElEmForIndex(ilm, &el, &em);
    outvec[ioutcount++] = invec[ilm*2];
    if (em == 0)
      continue;
    outvec[ioutcount++] = invec[ilm*2 + 1];
  }
  
  return ioutcount;
}

int AliFemtoCorrFctnDirectYlm::PackYlmMatrix(const double *inmat, double *outmat)
{
  // Pack a matrix in l,m x l,m into an array using
  // only independent components
  int ioutcountz = 0;
  int ioutcountp = 0;
  int emz, elz;
  int emp, elp;
  int finalsize = 0;

  for (int ilm=0; ilm<GetMaxJM(); ilm++) {
    GetElEmForIndex(ilm, &elz, &emz);
    finalsize++;
    if (emz == 0) continue;
    finalsize++;
  }

  for (int ilmz=0; ilmz<GetMaxJM(); ilmz++) {
    GetElEmForIndex(ilmz, &elz, &emz);
    ioutcountp = 0;
    for (int ilmp=0; ilmp<GetMaxJM(); ilmp++) {
      GetElEmForIndex(ilmp, &elp, &emp);
      outmat[ioutcountz*finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 0, ilmp, 0)];
      ioutcountp++;
      if (emp == 0) continue;
      outmat[ioutcountz*finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 0, ilmp, 1)];
      ioutcountp++;
    }
    ioutcountz++;

    if (emz == 0) continue;
    ioutcountp = 0;
    for (int ilmp=0; ilmp<GetMaxJM(); ilmp++) {
      GetElEmForIndex(ilmp, &elp, &emp);
      outmat[ioutcountz*finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 1, ilmp, 0)];
      ioutcountp++;
      if (emp == 0) continue;
      outmat[ioutcountz*finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 1, ilmp, 1)];
      ioutcountp++;
    }
    ioutcountz++;    
  }	
  
  return ioutcountz;  
}

void AliFemtoCorrFctnDirectYlm::PackCovariances()
{
  // Migrate the covariance matrix into a 3D histogram for storage
  char bufname[200];
  sprintf(bufname, "CovNum%s", fnumsreal[0]->GetName()+10);

  if (fcovnum) delete fcovnum;
  fcovnum = new TH3D(bufname,bufname, 
		    fnumsreal[0]->GetNbinsX(), fnumsreal[0]->GetXaxis()->GetXmin(), fnumsreal[0]->GetXaxis()->GetXmax(),
		    GetMaxJM()*2, -0.5, GetMaxJM()*2 - 0.5,
		    GetMaxJM()*2, -0.5, GetMaxJM()*2 - 0.5);
  
  for (int ibin=1; ibin<=fcovnum->GetNbinsX(); ibin++)
    for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++)
      for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++)
	fcovnum->SetBinContent(ibin, ilmz+1, ilmp+1, fcovmnum[GetBin(ibin-1, ilmz/2, ilmz%2, ilmp/2, ilmp%2)]);

  if (fcovden) delete fcovden;
  sprintf(bufname, "CovDen%s", fnumsreal[0]->GetName()+10);
  fcovden  = new TH3D(bufname,bufname, 
		     fdensreal[0]->GetNbinsX(), fdensreal[0]->GetXaxis()->GetXmin(), fdensreal[0]->GetXaxis()->GetXmax(),
		     GetMaxJM()*2, -0.5, GetMaxJM()*2 - 0.5,
		     GetMaxJM()*2, -0.5, GetMaxJM()*2 - 0.5);
		     
  for (int ibin=1; ibin<=fcovden->GetNbinsX(); ibin++)
    for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++)
      for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++)
	fcovden->SetBinContent(ibin, ilmz+1, ilmp+1, fcovmden[GetBin(ibin-1, ilmz/2, ilmz%2, ilmp/2, ilmp%2)]);

}

void AliFemtoCorrFctnDirectYlm::UnpackCovariances()
{
  // Extract the covariance matrices from storage
  if (fcovnum) {
    for (int ibin=1; ibin<=fcovnum->GetNbinsX(); ibin++)
      for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++)
	for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++)
	  fcovmnum[GetBin(ibin-1, ilmz/2, ilmz%2, ilmp/2, ilmp%2)] = fcovnum->GetBinContent(ibin, ilmz+1, ilmp+1);
    
  }
  if (fcovden) {
    for (int ibin=1; ibin<=fcovden->GetNbinsX(); ibin++)
      for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++)
	for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++)
	  fcovmden[GetBin(ibin-1, ilmz/2, ilmz%2, ilmp/2, ilmp%2)] = fcovden->GetBinContent(ibin, ilmz+1, ilmp+1);
  }
}

int AliFemtoCorrFctnDirectYlm::GetIndexForLM(int el, int em) const
{
  // Get array index for a given l,m
  for (int iter=0; iter<fMaxJM; iter++)
    if ((el == felsi[iter]) && (em == femsi[iter]))
      return iter;
  return -1;
}

TH1D *AliFemtoCorrFctnDirectYlm::GetNumRealHist(int el, int em)
{
  // Get numerator hist for a given l,m
  if (GetIndexForLM(el, em)>=0)
    return fnumsreal[GetIndexForLM(el, em)];
  else 
    return 0;
}
TH1D *AliFemtoCorrFctnDirectYlm::GetNumImagHist(int el, int em)
{
  // Get numerator hist for a given l,m
  if (GetIndexForLM(el, em)>=0)
    return fnumsimag[GetIndexForLM(el, em)];
  else 
    return 0;
}

TH1D *AliFemtoCorrFctnDirectYlm::GetDenRealHist(int el, int em)
{
  // Get denominator hist for a given l,m
  if (GetIndexForLM(el, em)>=0)
    return fdensreal[GetIndexForLM(el, em)];
  else 
    return 0;
}
TH1D *AliFemtoCorrFctnDirectYlm::GetDenImagHist(int el, int em)
{
  // Get denominator hist for a given l,m
  if (GetIndexForLM(el, em)>=0)
    return fdensimag[GetIndexForLM(el, em)];
  else 
    return 0;
}

AliFemtoString AliFemtoCorrFctnDirectYlm::Report()
{
  return "AliFemtoCorrFctnDirectYlm::Finish";
}

void AliFemtoCorrFctnDirectYlm::AddRealPair(AliFemtoPair* aPair)
{
  AddRealPair(aPair->QOutPf(), aPair->QSidePf(), aPair->QLongPf(), 1.0);
}
void AliFemtoCorrFctnDirectYlm::AddMixedPair(AliFemtoPair* aPair)
{
  AddMixedPair(aPair->QOutPf(), aPair->QSidePf(), aPair->QLongPf(), 1.0);
}

