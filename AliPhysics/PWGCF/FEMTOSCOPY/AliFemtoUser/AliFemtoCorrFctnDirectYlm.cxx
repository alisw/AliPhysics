///
/// \file AliFemtoUser/AliFemtoCorrFctnDirectYml.cxx
///


#include "AliFemtoCorrFctnDirectYlm.h"
#include <TMath.h>
#include <iostream>
#include <algorithm>

using namespace std;


AliFemtoCorrFctnDirectYlm::AliFemtoCorrFctnDirectYlm(const char *name,
                                                     int maxl,
                                                     int ibin,
                                                     double vmin,
                                                     double vmax,
                                                     int aUseLCMS)
  : AliFemtoCorrFctn()
  , fMaxL(maxl)
  , fMaxJM((maxl + 1) * (maxl + 1))
  , fnumsreal(fMaxJM, nullptr)
  , fnumsimag(fMaxJM, nullptr)
  , fdensreal(fMaxJM, nullptr)
  , fdensimag(fMaxJM, nullptr)
  , fbinctn(nullptr)
  , fbinctd(nullptr)
  , fcovnum(nullptr)
  , fcovden(nullptr)
  , fcovmnum(fMaxJM * fMaxJM * 4 * ibin)
  , fcovmden(fMaxJM * fMaxJM * 4 * ibin)
  , fels(fMaxJM)
  , fems(fMaxJM)
  , felsi(fMaxJM)
  , femsi(fMaxJM)
  , fYlmBuffer(fMaxJM)
  , factorials(4 * (maxl + 1))
  , fSout(0.0)
  , fSside(0.0)
  , fSlong(0.0)
  , fUseLCMS(aUseLCMS)
{
  // *DEB*  cout <<  "Size is " << sizeof(double) << " " << sizeof(complex<double>) << endl;

  // Fill in factorials table
  factorials[0] = 1;
  for (size_t iter = 1; iter < factorials.size(); iter++) {
    factorials[iter] = iter * factorials[iter - 1];
  }

  // Fill in els and ems table
  int el = 0;
  int em = 0;
  int il = 0;
  do {
    fels[il] = el;
    fems[il] = em;
    felsi[il] = (int)el;
    femsi[il] = (int)em;

    // *DEB*    cout << "il el em " << il << " " << felsi[il] << " " << femsi[il] << endl;
    em++;
    il++;
    if (em > el) {
      el++;
      em = -el;
    }
  } while (el <= maxl);

  // *DEB*  for (il=0; il<fMaxJM; il++)
  // *DEB*    cout << "il el em " << il << " " << felsi[il] << " " << femsi[il] << endl;

  // Create numerator and denominator historgrams
  for (int ihist = 0; ihist < fMaxJM; ihist++) {
    const TString suffix = TString::Format("Ylm%i%i%s",
                                           felsi[ihist],
                                           femsi[ihist] < 0 ? felsi[ihist] - femsi[ihist] : femsi[ihist],
                                           name);
    fnumsreal[ihist] = new TH1D("NumRe" + suffix, "NumRe" + suffix, ibin, vmin, vmax);
    fnumsimag[ihist] = new TH1D("NumIm" + suffix, "NumIm" + suffix, ibin, vmin, vmax);
    fdensreal[ihist] = new TH1D("DenRe" + suffix, "DenRe" + suffix, ibin, vmin, vmax);
    fdensimag[ihist] = new TH1D("DenIm" + suffix, "DenIm" + suffix, ibin, vmin, vmax);

    fnumsreal[ihist]->Sumw2();
    fnumsimag[ihist]->Sumw2();
    fdensreal[ihist]->Sumw2();
    fdensimag[ihist]->Sumw2();
  }

  fbinctn = new TH1D(TString("BinCountNum") + name, "Bin Occupation (Numerator)", ibin, vmin, vmax);
  fbinctd = new TH1D(TString("BinCountDen") + name, "Bin Occupation (Denominator)", ibin, vmin, vmax);

  AliFemtoYlm::InitializeYlms();
}

AliFemtoCorrFctnDirectYlm::AliFemtoCorrFctnDirectYlm()
  : AliFemtoCorrFctnDirectYlm("AliFemtoCorrFctnDirectYlm", 2)
{
}

AliFemtoCorrFctnDirectYlm::AliFemtoCorrFctnDirectYlm(const AliFemtoCorrFctnDirectYlm& aCorrFctn)
  : AliFemtoCorrFctn(aCorrFctn)
  , fMaxL(aCorrFctn.fMaxL)
  , fMaxJM((fMaxL + 1) * (fMaxL + 1))
  , fnumsreal(fMaxJM, nullptr)
  , fnumsimag(fMaxJM, nullptr)
  , fdensreal(fMaxJM, nullptr)
  , fdensimag(fMaxJM, nullptr)
  , fbinctn(new TH1D(*aCorrFctn.fbinctn))
  , fbinctd(new TH1D(*aCorrFctn.fbinctd))
  , fcovnum(aCorrFctn.fcovnum ? new TH3D(*aCorrFctn.fcovnum) : nullptr)
  , fcovden(aCorrFctn.fcovden ? new TH3D(*aCorrFctn.fcovden) : nullptr)
  , fcovmnum(aCorrFctn.fcovmnum)
  , fcovmden(aCorrFctn.fcovmden)
  , fels(aCorrFctn.fels)
  , fems(aCorrFctn.fems)
  , felsi(aCorrFctn.felsi)
  , femsi(aCorrFctn.femsi)
  , fYlmBuffer(aCorrFctn.fYlmBuffer)
  , factorials(aCorrFctn.factorials)
  , fSout(aCorrFctn.fSout)
  , fSside(aCorrFctn.fSside)
  , fSlong(aCorrFctn.fSlong)
  , fUseLCMS(aCorrFctn.fUseLCMS)
{
  // Copy constructor
  for (int ihist = 0; ihist < fMaxJM; ihist++) {
    fnumsreal[ihist] = new TH1D(*aCorrFctn.fnumsreal[ihist]);
    fnumsimag[ihist] = new TH1D(*aCorrFctn.fnumsimag[ihist]);
    fdensreal[ihist] = new TH1D(*aCorrFctn.fdensreal[ihist]);
    fdensimag[ihist] = new TH1D(*aCorrFctn.fdensimag[ihist]);
  }
}

AliFemtoCorrFctnDirectYlm&
AliFemtoCorrFctnDirectYlm::operator=(const AliFemtoCorrFctnDirectYlm& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  AliFemtoCorrFctn::operator=(aCorrFctn);

  fMaxL = aCorrFctn.fMaxL;
  fMaxJM = (fMaxL + 1) * (fMaxL + 1);

  // Copy factorials table if necessary
  if (factorials.size() != aCorrFctn.factorials.size()) {
    factorials = aCorrFctn.factorials;
  }

  // Copy in els and ems table
  if (fels.size() != aCorrFctn.fels.size()) {
    fels = aCorrFctn.fels;
    fems = aCorrFctn.fems;
    felsi = aCorrFctn.felsi;
    femsi = aCorrFctn.femsi;
  }


  if (fnumsreal.size() == aCorrFctn.fnumsreal.size()) {
    // simple assignments if same number of histograms

    for (int ihist = 0; ihist < fMaxJM; ihist++) {
      *fnumsreal[ihist] = *aCorrFctn.fnumsreal[ihist];
      *fnumsimag[ihist] = *aCorrFctn.fnumsimag[ihist];
      *fdensreal[ihist] = *aCorrFctn.fdensreal[ihist];
      *fdensimag[ihist] = *aCorrFctn.fdensimag[ihist];
    }

  } else {
    // delete and resize to match
    for (size_t i = 0; i < fnumsreal.size(); ++i) {
      delete fnumsreal[i];
      delete fnumsimag[i];
      delete fdensreal[i];
      delete fdensimag[i];
    }

    fnumsreal.resize(fMaxJM, nullptr);
    fnumsimag.resize(fMaxJM, nullptr);
    fdensreal.resize(fMaxJM, nullptr);
    fdensimag.resize(fMaxJM, nullptr);

    for (int ihist = 0; ihist < fMaxJM; ihist++) {
      fnumsreal[ihist] = new TH1D(*aCorrFctn.fnumsreal[ihist]);
      fnumsimag[ihist] = new TH1D(*aCorrFctn.fnumsimag[ihist]);
      fdensreal[ihist] = new TH1D(*aCorrFctn.fdensreal[ihist]);
      fdensimag[ihist] = new TH1D(*aCorrFctn.fdensimag[ihist]);
    }
  }

  fbinctn = new TH1D(*aCorrFctn.fbinctn);
  fbinctd = new TH1D(*aCorrFctn.fbinctd);

  fYlmBuffer = aCorrFctn.fYlmBuffer;

  // Covariance matrices
  fcovmnum = aCorrFctn.fcovmnum;
  fcovmden = aCorrFctn.fcovmden;

  if (aCorrFctn.fcovnum && fcovnum) {
    *fcovnum = *aCorrFctn.fcovnum;
  } else if (aCorrFctn.fcovnum) {
    fcovnum = new TH3D(*aCorrFctn.fcovnum);
  } else {
    delete fcovnum;
    fcovnum = nullptr;
  }

  if (aCorrFctn.fcovden && fcovden) {
    *fcovden = *aCorrFctn.fcovden;
  } else if (aCorrFctn.fcovden) {
    fcovden = new TH3D(*aCorrFctn.fcovden);
  } else {
    delete fcovden;
    fcovden = nullptr;
  }

  fSout = aCorrFctn.fSout;
  fSside = aCorrFctn.fSside;
  fSlong = aCorrFctn.fSlong;

  fUseLCMS = aCorrFctn.fUseLCMS;

  return *this;
}

AliFemtoCorrFctnDirectYlm::~AliFemtoCorrFctnDirectYlm()
{
  // Destructor
  for (int ihist = 0; ihist < fMaxJM; ihist++) {
    delete fnumsreal[ihist];
    delete fnumsimag[ihist];
    delete fdensreal[ihist];
    delete fdensimag[ihist];
  }

  delete fbinctn;
  delete fbinctd;

  delete fcovnum;
  delete fcovden;
}

double
AliFemtoCorrFctnDirectYlm::ClebschGordan(double aJot1,
                                         double aEm1,
                                         double aJot2,
                                         double aEm2,
                                         double aJot,
                                         double aEm)
{
  // Calculate Clebsh-Gordan coefficient
  double cgc = 0.0;

  const int maxt = std::min({ lrint(aJot1 + aJot2 - aJot),
                              lrint(aJot2 + aEm2),
                              lrint(aJot1 - aEm1) }),

            mint = std::max({ 0l,
                              lrint(-(aJot - aJot2 + aEm1)),
                              lrint(-(aJot - aJot1 - aEm2)) });

  for (int titer = mint; titer <= maxt; titer++) {
    double coef = titer % 2 ? -1.0 : 1.0;  // (-1)^titer
    coef *= TMath::Sqrt((2 * aJot + 1)
                        * factorials[lrint(aJot1 + aEm1)]
                        * factorials[lrint(aJot1 - aEm1)]
                        * factorials[lrint(aJot2 + aEm2)]
                        * factorials[lrint(aJot2 - aEm2)]
                        * factorials[lrint(aJot + aEm)]
                        * factorials[lrint(aJot - aEm)]);
    coef /= (factorials[titer]
             * factorials[lrint(aJot1 + aJot2 - aJot - titer)]
             * factorials[lrint(aJot1 - aEm1 - titer)]
             * factorials[lrint(aJot2 + aEm2 - titer)]
             * factorials[lrint(aJot - aJot2 + aEm1 + titer)]
             * factorials[lrint(aJot - aJot1 - aEm2 + titer)]);

    cgc += coef;
  }

  cgc *= DeltaJ(aJot1, aJot2, aJot);

  return cgc;
}

double
AliFemtoCorrFctnDirectYlm::DeltaJ(double aJot1, double aJot2, double aJot)
{
  // Calculate J for the Clebsh-Gordan coefficient
  if ((aJot1 + aJot2 - aJot) < 0) {
    //    cout << "J1+J2-J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((aJot1 - aJot2 + aJot) < 0) {
    //    cout << "J1-J2+J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((-aJot1 + aJot2 + aJot) < 0) {
    //    cout << "-J1+J2+J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((aJot1 + aJot2 + aJot + 1) < 0) {
    //    cout << "J1+J2+J3+1 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  double res = TMath::Sqrt(1.0 * factorials[lrint(aJot1 + aJot2 - aJot)]
                               * factorials[lrint(aJot1 - aJot2 + aJot)]
                               * factorials[lrint(-aJot1 + aJot2 + aJot)]
                               / factorials[lrint(aJot1 + aJot2 + aJot + 1)]);

  return res;
}

double
AliFemtoCorrFctnDirectYlm::WignerSymbol(double aJot1,
                                        double aEm1,
                                        double aJot2,
                                        double aEm2,
                                        double aJot,
                                        double aEm)
{
  // Get Wigner symbol
  if (lrint(aEm1 + aEm2 + aEm) != 0.0) {
    return 0.0;
  }
  double cge = ClebschGordan(aJot1, aEm1, aJot2, aEm2, aJot, -aEm);
  if (lrint(abs(aJot1 - aJot2 - aEm)) % 2) {
    cge *= -1.0;
  }
  cge /= sqrt(2 * aJot + 1);

  if (cge == -0.0) {
    cge = 0.0;
  }

  return cge;
}

void
AliFemtoCorrFctnDirectYlm::GetMtilde(complex<double> *aMat, double *aMTilde)
{
  // Create the Mtilde for a given q bin
  double lzero, mzero;
  double lprim, mprim;
  double lbis, mbis;

  int lzeroi, mzeroi;
  int lprimi, mprimi;
  int lbisi, mbisi;

  for (int izero = 0; izero < GetMaxJM(); izero++) {
    GetElEmForIndex(izero, &lzero, &mzero);
    GetElEmForIndex(izero, &lzeroi, &mzeroi);

    const double neg_one_to_m = abs(mzeroi) % 2 ? -1.0 : 1.0;

    for (int ibis = 0; ibis < GetMaxJM(); ibis++) {
      GetElEmForIndex(ibis, &lbis, &mbis);
      GetElEmForIndex(ibis, &lbisi, &mbisi);

      std::complex<double> val(0.0, 0.0);
      for (int iprim = 0; iprim < GetMaxJM(); iprim++) {

        GetElEmForIndex(iprim, &lprim, &mprim);
        GetElEmForIndex(iprim, &lprimi, &mprimi);

        std::complex<double> mcomp(neg_one_to_m, 0.0);  // (-1)^m

        mcomp *= sqrt((2 * lzero + 1) * (2 * lprim + 1) * (2 * lbis + 1)); // P1
        mcomp *= WignerSymbol(lzero, 0, lprim, 0, lbis, 0);                // W1
        mcomp *= WignerSymbol(lzero, -mzero, lprim, mprim, lbis, mbis);    // W2
        mcomp *= aMat[iprim];
        val += mcomp;
      }
      aMTilde[(izero * 2) * (2 * GetMaxJM()) + (ibis * 2)] = real(val);
      aMTilde[(izero * 2 + 1) * (2 * GetMaxJM()) + (ibis * 2)] = imag(val);
      aMTilde[(izero * 2) * (2 * GetMaxJM()) + (ibis * 2 + 1)] = imag(val) ? -val.imag() : 0.0;
      aMTilde[(izero * 2 + 1) * (2 * GetMaxJM()) + (ibis * 2 + 1)] = real(val);
    }
  }
}

int
AliFemtoCorrFctnDirectYlm::GetMaxJM() const
{
  return fMaxJM;
}

void
AliFemtoCorrFctnDirectYlm::GetElEmForIndex(int aIndex, double *aEl, double *aEm) const
{
  // Get l,m for a given index
  *aEl = fels[aIndex];
  *aEm = fems[aIndex];
}

void
AliFemtoCorrFctnDirectYlm::GetElEmForIndex(int aIndex, int *aEl, int *aEm) const
{
  // Get l,m for a given index
  *aEl = felsi[aIndex];
  *aEm = femsi[aIndex];
}

int
AliFemtoCorrFctnDirectYlm::GetBin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag)
{
  return qbin * GetMaxJM() * GetMaxJM() * 4
         + (ilmprim * 2 + primimag) * GetMaxJM() * 2
         + ilmzero * 2
         + zeroimag;
}

void
AliFemtoCorrFctnDirectYlm::AddRealPair(double qout, double qside, double qlong, double weight)
{
  // Fill numerator
  double kv = sqrt(qout * qout + qside * qside + qlong * qlong);
  int nqbin = fbinctn->GetXaxis()->FindFixBin(kv) - 1;

  // Use saved ylm values for same qout, qside, qlong
  if ((qout != fSout) || (qside != fSside) || (qlong != fSlong)) {
    AliFemtoYlm::YlmUpToL(fMaxL, qout, qside, qlong, fYlmBuffer.data());
    fSout = qout;
    fSside = qside;
    fSlong = qlong;
  }
  for (int ilm = 0; ilm < GetMaxJM(); ilm++) {
    //    fYlmBuffer[ilm] = AliFemtoYlm::Ylm(elsi[ilm], emsi[ilm], qout, qside, qlong);

    fnumsreal[ilm]->Fill(kv, real(fYlmBuffer[ilm]) * weight);
    fnumsimag[ilm]->Fill(kv, -imag(fYlmBuffer[ilm]) * weight);

    fbinctn->Fill(kv, 1.0);
  }

  // Fill in the error matrix
  //  int tabshift = nqbin*GetMaxJM()*GetMaxJM()*4;
  if (nqbin < fbinctn->GetNbinsX())
    for (int ilmzero = 0; ilmzero < GetMaxJM(); ilmzero++)
      for (int ilmprim = 0; ilmprim < GetMaxJM(); ilmprim++) {
        fcovmnum[GetBin(nqbin, ilmzero, 0, ilmprim, 0)] +=
          real(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]) * weight * weight;
        fcovmnum[GetBin(nqbin, ilmzero, 0, ilmprim, 1)] +=
          real(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]) * weight * weight;
        fcovmnum[GetBin(nqbin, ilmzero, 1, ilmprim, 0)] +=
          -imag(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]) * weight * weight;
        fcovmnum[GetBin(nqbin, ilmzero, 1, ilmprim, 1)] +=
          -imag(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]) * weight * weight;
      }
}

void
AliFemtoCorrFctnDirectYlm::AddMixedPair(double qout, double qside, double qlong, double weight)
{
  // Fill denominator
  double kv = sqrt(qout * qout + qside * qside + qlong * qlong);

  // Use saved ylm values for same qout, qside, qlong
  if ((qout != fSout) || (qside != fSside) || (qlong != fSlong)) {
    AliFemtoYlm::YlmUpToL(fMaxL, qout, qside, qlong, fYlmBuffer.data());
    fSout = qout;
    fSside = qside;
    fSlong = qlong;
  }
  for (int ilm = 0; ilm < GetMaxJM(); ilm++) {
    //    fYlmBuffer[ilm] = AliFemtoYlm::Ylm(elsi[ilm], emsi[ilm], qout, qside, qlong);

    fdensreal[ilm]->Fill(kv, real(fYlmBuffer[ilm]) * weight);
    fdensimag[ilm]->Fill(kv, -imag(fYlmBuffer[ilm]) * weight);

    fbinctd->Fill(kv, 1.0);
  }

  // Fill in the error matrix
  int nqbin = fbinctn->GetXaxis()->FindFixBin(kv) - 1;
  //  int tabshift = nqbin*GetMaxJM()*GetMaxJM()*4;
  if (nqbin < fbinctn->GetNbinsX())
    for (int ilmzero = 0; ilmzero < GetMaxJM(); ilmzero++)
      for (int ilmprim = 0; ilmprim < GetMaxJM(); ilmprim++) {
        fcovmden[GetBin(nqbin, ilmzero, 0, ilmprim, 0)] +=
          real(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]);
        fcovmden[GetBin(nqbin, ilmzero, 0, ilmprim, 1)] +=
          real(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]);
        fcovmden[GetBin(nqbin, ilmzero, 1, ilmprim, 0)] +=
          -imag(fYlmBuffer[ilmzero]) * real(fYlmBuffer[ilmprim]);
        fcovmden[GetBin(nqbin, ilmzero, 1, ilmprim, 1)] +=
          -imag(fYlmBuffer[ilmzero]) * -imag(fYlmBuffer[ilmprim]);
      }
}

void
AliFemtoCorrFctnDirectYlm::AddRealPair(double *qvec, double weight)
{
  AddRealPair(qvec[0], qvec[1], qvec[2], weight);
}

void
AliFemtoCorrFctnDirectYlm::AddMixedPair(double *qvec, double weight)
{
  AddMixedPair(qvec[0], qvec[1], qvec[2], weight);
}

void
AliFemtoCorrFctnDirectYlm::Finish()
{
  PackCovariances();
}

void
AliFemtoCorrFctnDirectYlm::Write()
{
  // Write out output histograms
  if ((!fcovnum) || (!fcovden))
    PackCovariances();

  for (int ilm = 0; ilm < fMaxJM; ilm++) {
    fnumsreal[ilm]->Write();
    fdensreal[ilm]->Write();
    fnumsimag[ilm]->Write();
    fdensimag[ilm]->Write();
  }
  if (fcovnum)
    fcovnum->Write();
  if (fcovden)
    fcovden->Write();
}

TList*
AliFemtoCorrFctnDirectYlm::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  if ((!fcovnum) || (!fcovden))
    PackCovariances();

  TList *tOutputList = new TList();

  for (int ilm = 0; ilm < fMaxJM; ilm++) {
    tOutputList->Add(fnumsreal[ilm]);
    tOutputList->Add(fdensreal[ilm]);
    tOutputList->Add(fnumsimag[ilm]);
    tOutputList->Add(fdensimag[ilm]);
  }
  if (fcovnum)
    tOutputList->Add(fcovnum);
  if (fcovden)
    tOutputList->Add(fcovden);

  return tOutputList;
}

void
AliFemtoCorrFctnDirectYlm::ReadFromFile(TFile *infile, const char *name, int maxl)
{
  // Raad in the numerator and denominator from file
  if (maxl != fMaxL) {
    cout << "Cannot read function for L " << maxl << " into a container with L " << fMaxL << endl;
    return;
  }
  cout << "Reading in numerators and denominators" << endl;

  for (int ihist = 0; ihist < fMaxJM; ihist++) {

    TString suffix = TString::Format("Ylm%i%i%s",
                                     felsi[ihist],
                                     femsi[ihist] < 0 ? felsi[ihist] - femsi[ihist] : femsi[ihist],
                                     name);
    if (fnumsreal[ihist])
      delete fnumsreal[ihist];
    fnumsreal[ihist] = new TH1D(*((TH1D*)infile->Get("NumRe" + suffix)));

    if (fnumsimag[ihist])
      delete fnumsimag[ihist];
    fnumsimag[ihist] = new TH1D(*((TH1D*)infile->Get("NumIm" + suffix)));

    if (fdensreal[ihist])
      delete fdensreal[ihist];
    fdensreal[ihist] = new TH1D(*((TH1D*)infile->Get("DenRe" + suffix)));

    if (fdensimag[ihist])
      delete fdensimag[ihist];
    fdensimag[ihist] = new TH1D(*((TH1D*)infile->Get("DenIm" + suffix)));
  }

  if (fcovnum)
    delete fcovnum;
  fcovnum = new TH3D(*((TH3D*)infile->Get(Form("CovNum%s", name))));

  if (fcovden)
    delete fcovden;
  fcovden = new TH3D(*((TH3D*)infile->Get(Form("CovDen%s", name))));

  if ((fcovnum) && (fcovden)) {
    cout << "Unpacking covariance matrices from file " << endl;
    UnpackCovariances();
  } else {

    cout << "Creating fake covariance matrices" << endl;

    for (int ibin = 1; ibin <= fnumsreal[0]->GetNbinsX(); ibin++) {
      double nent = fnumsreal[0]->GetEntries();
      double nentd = fdensreal[0]->GetEntries();
      for (int ilmx = 0; ilmx < GetMaxJM(); ilmx++) {
        for (int ilmy = 0; ilmy < GetMaxJM(); ilmy++) {
          double t1t2rr =
            fnumsreal[ilmx]->GetBinContent(ibin) * fnumsreal[ilmy]->GetBinContent(ibin) / nent / nent;
          double t1t2ri =
            fnumsreal[ilmx]->GetBinContent(ibin) * fnumsimag[ilmy]->GetBinContent(ibin) / nent / nent;
          double t1t2ir =
            fnumsimag[ilmx]->GetBinContent(ibin) * fnumsreal[ilmy]->GetBinContent(ibin) / nent / nent;
          double t1t2ii =
            fnumsimag[ilmx]->GetBinContent(ibin) * fnumsimag[ilmy]->GetBinContent(ibin) / nent / nent;
          if (ilmx == ilmy) {
            fcovmnum[GetBin(ibin - 1, ilmx, 0, ilmy, 0)] =
              nent * (TMath::Power(fnumsreal[ilmx]->GetBinError(ibin) / nent, 2) * (nent - 1) + t1t2rr);
            fcovmnum[GetBin(ibin - 1, ilmx, 0, ilmy, 1)] = nent * t1t2ri;
            fcovmnum[GetBin(ibin - 1, ilmx, 1, ilmy, 0)] = nent * t1t2ir;
            fcovmnum[GetBin(ibin - 1, ilmx, 1, ilmy, 1)] =
              nent * (TMath::Power(fnumsimag[ilmx]->GetBinError(ibin) / nent, 2) * (nent - 1) + t1t2rr);
          } else {
            fcovmnum[GetBin(ibin - 1, ilmx, 0, ilmy, 0)] = nent * t1t2rr;
            fcovmnum[GetBin(ibin - 1, ilmx, 0, ilmy, 1)] = nent * t1t2ri;
            fcovmnum[GetBin(ibin - 1, ilmx, 1, ilmy, 0)] = nent * t1t2ir;
            fcovmnum[GetBin(ibin - 1, ilmx, 1, ilmy, 1)] = nent * t1t2ii;
          }
          t1t2rr =
            fdensreal[ilmx]->GetBinContent(ibin) * fdensreal[ilmy]->GetBinContent(ibin) / nentd / nentd;
          t1t2ri =
            fdensreal[ilmx]->GetBinContent(ibin) * fdensimag[ilmy]->GetBinContent(ibin) / nentd / nentd;
          t1t2ir =
            fdensimag[ilmx]->GetBinContent(ibin) * fdensreal[ilmy]->GetBinContent(ibin) / nentd / nentd;
          t1t2ii =
            fdensimag[ilmx]->GetBinContent(ibin) * fdensimag[ilmy]->GetBinContent(ibin) / nentd / nentd;

          fcovmden[GetBin(ibin - 1, ilmx, 0, ilmy, 0)] = nentd * t1t2rr;
          fcovmden[GetBin(ibin - 1, ilmx, 0, ilmy, 1)] = nentd * t1t2ri;
          fcovmden[GetBin(ibin - 1, ilmx, 1, ilmy, 0)] = nentd * t1t2ir;
          fcovmden[GetBin(ibin - 1, ilmx, 1, ilmy, 1)] = nentd * t1t2ii;
        }
      }
    }
  }

  // Recalculating the correlation functions
  Finish();
}

int
AliFemtoCorrFctnDirectYlm::PackYlmVector(const double *invec, double *outvec)
{
  // Pack a vector in l,m into an array using
  // only independent components
  int ioutcount = 0;
  int em, el;
  for (int ilm = 0; ilm < GetMaxJM(); ilm++) {
    GetElEmForIndex(ilm, &el, &em);
    outvec[ioutcount++] = invec[ilm * 2];
    if (em == 0)
      continue;
    outvec[ioutcount++] = invec[ilm * 2 + 1];
  }

  return ioutcount;
}

int
AliFemtoCorrFctnDirectYlm::PackYlmMatrix(const double *inmat, double *outmat)
{
  // Pack a matrix in l,m x l,m into an array using
  // only independent components
  int ioutcountz = 0;
  int ioutcountp = 0;
  int emz, elz;
  int emp, elp;
  int finalsize = 0;

  for (int ilm = 0; ilm < GetMaxJM(); ilm++) {
    GetElEmForIndex(ilm, &elz, &emz);
    finalsize++;
    if (emz == 0)
      continue;
    finalsize++;
  }

  for (int ilmz = 0; ilmz < GetMaxJM(); ilmz++) {
    GetElEmForIndex(ilmz, &elz, &emz);
    ioutcountp = 0;
    for (int ilmp = 0; ilmp < GetMaxJM(); ilmp++) {
      GetElEmForIndex(ilmp, &elp, &emp);
      outmat[ioutcountz * finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 0, ilmp, 0)];
      ioutcountp++;
      if (emp == 0)
        continue;
      outmat[ioutcountz * finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 0, ilmp, 1)];
      ioutcountp++;
    }
    ioutcountz++;

    if (emz == 0)
      continue;
    ioutcountp = 0;
    for (int ilmp = 0; ilmp < GetMaxJM(); ilmp++) {
      GetElEmForIndex(ilmp, &elp, &emp);
      outmat[ioutcountz * finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 1, ilmp, 0)];
      ioutcountp++;
      if (emp == 0)
        continue;
      outmat[ioutcountz * finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 1, ilmp, 1)];
      ioutcountp++;
    }
    ioutcountz++;
  }

  return ioutcountz;
}

void
AliFemtoCorrFctnDirectYlm::PackCovariances()
{
  // Migrate the covariance matrix into a 3D histogram for storage

  //  if (fcovnum) delete fcovnum;
  if (!fcovnum) {
    auto *nax = fnumsreal[0]->GetXaxis();
    const char *bufname = Form("CovNum%s", fnumsreal[0]->GetName() + 10);
    fcovnum = new TH3D(bufname, bufname,
                       fnumsreal[0]->GetNbinsX(), nax->GetXmin(), nax->GetXmax(),
                       GetMaxJM() * 2, -0.5, GetMaxJM() * 2 - 0.5,
                       GetMaxJM() * 2, -0.5, GetMaxJM() * 2 - 0.5);
  }

  for (int ibin = 1; ibin <= fcovnum->GetNbinsX(); ibin++)
    for (int ilmz = 0; ilmz < GetMaxJM() * 2; ilmz++)
      for (int ilmp = 0; ilmp < GetMaxJM() * 2; ilmp++) {
        auto bin = GetBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
        auto value = fcovmnum[bin];
        fcovnum->SetBinContent(ibin, ilmz + 1, ilmp + 1, value);
      }

  //  if (fcovden) delete fcovden;
  if (!fcovden) {
    auto *dax = fdensreal[0]->GetXaxis();
    const char *bufname = Form("CovDen%s", fnumsreal[0]->GetName() + 10);
    fcovden = new TH3D(bufname, bufname,
                       fdensreal[0]->GetNbinsX(), dax->GetXmin(), dax->GetXmax(),
                       GetMaxJM() * 2, -0.5, GetMaxJM() * 2 - 0.5,
                       GetMaxJM() * 2, -0.5, GetMaxJM() * 2 - 0.5);
  }

  for (int ibin = 1; ibin <= fcovden->GetNbinsX(); ibin++)
    for (int ilmz = 0; ilmz < GetMaxJM() * 2; ilmz++)
      for (int ilmp = 0; ilmp < GetMaxJM() * 2; ilmp++) {
        auto bin = GetBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2);
        auto value = fcovmden[bin];
        fcovden->SetBinContent(ibin, ilmz + 1, ilmp + 1, value);
      }
}

void
AliFemtoCorrFctnDirectYlm::UnpackCovariances()
{
  // Extract the covariance matrices from storage
  if (fcovnum) {
    for (int ibin = 1; ibin <= fcovnum->GetNbinsX(); ibin++)
      for (int ilmz = 0; ilmz < GetMaxJM() * 2; ilmz++)
        for (int ilmp = 0; ilmp < GetMaxJM() * 2; ilmp++)
          fcovmnum[GetBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2)] =
            fcovnum->GetBinContent(ibin, ilmz + 1, ilmp + 1);
  }
  if (fcovden) {
    for (int ibin = 1; ibin <= fcovden->GetNbinsX(); ibin++)
      for (int ilmz = 0; ilmz < GetMaxJM() * 2; ilmz++)
        for (int ilmp = 0; ilmp < GetMaxJM() * 2; ilmp++)
          fcovmden[GetBin(ibin - 1, ilmz / 2, ilmz % 2, ilmp / 2, ilmp % 2)] =
            fcovden->GetBinContent(ibin, ilmz + 1, ilmp + 1);
  }
}

int
AliFemtoCorrFctnDirectYlm::GetIndexForLM(int el, int em) const
{
  // Get array index for a given l,m
  for (int iter = 0; iter < fMaxJM; iter++) {
    if ((el == felsi[iter]) && (em == femsi[iter])) {
      return iter;
    }
  }
  return -1;
}

TH1D*
AliFemtoCorrFctnDirectYlm::GetNumRealHist(int el, int em)
{
  // Get numerator hist for a given l,m
  if (GetIndexForLM(el, em) >= 0) {
    return fnumsreal[GetIndexForLM(el, em)];
  }
  return nullptr;
}
TH1D*
AliFemtoCorrFctnDirectYlm::GetNumImagHist(int el, int em)
{
  // Get numerator hist for a given l,m
  if (GetIndexForLM(el, em) >= 0) {
    return fnumsimag[GetIndexForLM(el, em)];
  }
  return nullptr;
}

TH1D*
AliFemtoCorrFctnDirectYlm::GetDenRealHist(int el, int em)
{
  // Get denominator hist for a given l,m
  if (GetIndexForLM(el, em) >= 0) {
    return fdensreal[GetIndexForLM(el, em)];
  }
  return nullptr;
}
TH1D*
AliFemtoCorrFctnDirectYlm::GetDenImagHist(int el, int em)
{
  // Get denominator hist for a given l,m
  if (GetIndexForLM(el, em) >= 0) {
    return fdensimag[GetIndexForLM(el, em)];
  }
  return nullptr;
}

AliFemtoString
AliFemtoCorrFctnDirectYlm::Report()
{
  return "AliFemtoCorrFctnDirectYlm::Finish";
}

void
AliFemtoCorrFctnDirectYlm::AddRealPair(AliFemtoPair *aPair)
{
  // Fill in the numerator
  if (fPairCut && !fPairCut->Pass(aPair)) {
    return;
  }

  if (fUseLCMS)
    AddRealPair(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), 1.0);
  else
    AddRealPair(aPair->KOut(), aPair->KSide(), aPair->KLong(), 1.0);
}

void
AliFemtoCorrFctnDirectYlm::AddMixedPair(AliFemtoPair *aPair)
{
  // Fill in the denominator
  if (fPairCut && !fPairCut->Pass(aPair)) {
    return;
  }

  if (fUseLCMS)
    AddMixedPair(aPair->QOutCMS(), aPair->QSideCMS(), aPair->QLongCMS(), 1.0);
  else
    AddMixedPair(aPair->KOut(), aPair->KSide(), aPair->KLong(), 1.0);
}

void
AliFemtoCorrFctnDirectYlm::SetUseLCMS(int aUseLCMS)
{
  fUseLCMS = aUseLCMS;
}

int
AliFemtoCorrFctnDirectYlm::GetUseLCMS()
{
  return fUseLCMS;
}
