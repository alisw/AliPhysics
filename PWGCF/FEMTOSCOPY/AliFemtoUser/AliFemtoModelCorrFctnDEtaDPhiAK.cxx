///
/// \file AliFemtoModelCorrFctnDEtaDPhiAK.cxx
///

#include "AliFemtoModelCorrFctnDEtaDPhiAK.h"

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"

#include <TMath.h>

#include <tuple>
#include <cstdio>


AliFemtoModelCorrFctnDEtaDPhiAK
  ::AliFemtoModelCorrFctnDEtaDPhiAK(const Params &p)
  : AliFemtoModelCorrFctnDEtaDPhiAK(p.suffix,
                                    p.nbins_phi,
                                    p.nbins_eta)
{}

AliFemtoModelCorrFctnDEtaDPhiAK
  ::AliFemtoModelCorrFctnDEtaDPhiAK(const char* suffix,
                                    const UInt_t nbins_phi,
                                    const UInt_t nbins_eta)
  : AliFemtoModelCorrFctn()
  , fDPhiDEtaNumeratorTrue(nullptr)
  , fDPhiDEtaNumeratorFake(nullptr)
  , fDPhiDEtaDenominator(nullptr)
  , fDPhiDEtaColNumerator(nullptr)
  , fDPhiDEtaColDenominator(nullptr)
  , fDPhiNumeratorTrue(nullptr)
  , fDPhiNumeratorFake(nullptr)
  , fDPhiDenominator(nullptr)
  , fDCosNumeratorTrue(nullptr)
  , fDCosNumeratorFake(nullptr)
  , fDCosDenominator(nullptr)
  , fDPhiPtNumerator(nullptr)
  , fDPhiPtDenominator(nullptr)
  , fDCosPtNumerator(nullptr)
  , fDCosPtDenominator(nullptr)
{
  const double
    lo_phi = -M_PI_2,
    hi_phi = 3 * M_PI_2,
    lo_eta = -2.0,
    hi_eta = 2.0,
    lo_pt = 0.0,
    hi_pt = 3.0;

  const Int_t
    nbins_pt = 30;

  // Δɸ, Δη

  fDPhiDEtaNumeratorTrue = new TH2D(Form("NumDPhiDEtaTrue%s", suffix),
                                    "Numerator (same event); #Delta#phi; #Delta#eta",
                                    nbins_phi, lo_phi, hi_phi,
                                    nbins_eta, lo_eta, hi_eta);

  fDPhiDEtaNumeratorFake = new TH2D(Form("NumDPhiDEtaFake%s", suffix),
                                    "Numerator (mixed event); #Delta#phi; #Delta#eta",
                                    nbins_phi, lo_phi, hi_phi,
                                    nbins_eta, lo_eta, hi_eta);

  fDPhiDEtaDenominator = new TH2D(Form("DenDPhiDEta%s", suffix),
                                  "Denominator; #Delta#phi; #Delta#eta",
                                  nbins_phi, lo_phi, hi_phi,
                                  nbins_eta, lo_eta, hi_eta);

  fDPhiDEtaColNumerator = new TH2D(Form("NumDPhiDEtaCol%s", suffix),
                                   "Numerator; #Delta#phi; #Delta#eta Colinear",
                                   nbins_phi, lo_phi, hi_phi,
                                   nbins_eta, lo_eta, hi_eta);

  fDPhiDEtaColDenominator = new TH2D(Form("DenDPhiDEtaCol%s", suffix),
                                     "Denominator; #Delta#phi; #Delta#eta Colinear",
                                     nbins_phi, lo_phi, hi_phi,
                                     nbins_eta, lo_eta, hi_eta);

  // Δɸ

  fDPhiNumeratorTrue = new TH1D(Form("NumDPhiTrue%s", suffix),
                                "Numerator (same event); #Delta#phi",
                                2*nbins_phi, lo_phi, hi_phi);

  fDPhiNumeratorFake = new TH1D(Form("NumDPhiFake%s", suffix),
                                "Numerator (mixed event); #Delta#phi",
                                2*nbins_phi, lo_phi, hi_phi);

  fDPhiDenominator = new TH1D(Form("DenDPhi%s", suffix),
                              "Denominator; #Delta#phi",
                              2*nbins_phi, lo_phi, hi_phi);

  // cos(Δɸ)

  fDCosNumeratorTrue = new TH1D(Form("NumDCosTrue%s", suffix),
                                "Numerator Cosine Phi (same event); cos(#Delta#phi);",
                                2*nbins_phi, -1.0, 1.0);

  fDCosNumeratorFake = new TH1D(Form("NumDCosFake%s", suffix),
                                "Numerator Cosine Phi (mixed event); cos(#Delta#phi);",
                                2*nbins_phi, -1.0, 1.0);

  fDCosDenominator = new TH1D(Form("DenDCos%s", suffix),
                              "Denominator Cosine Phi; cos(#Delta#phi);",
                              2*nbins_phi, -1.0, 1.0);

  // Δɸ, pT

  fDPhiPtNumerator = new TH2D(Form("NumDPhiPt%s", suffix),
                              "Numerator; #DeltaPhi; p_{T}",
                              2*nbins_phi, lo_phi, hi_phi,
                              nbins_pt, lo_pt, hi_pt);

  fDPhiPtDenominator = new TH2D(Form("DenDPhiPt%s", suffix),
                                "Denominator; #DeltaPhi; p_{T}",
                                2*nbins_phi, lo_phi, hi_phi,
                                nbins_pt, lo_pt, hi_pt);

  // cos(Δɸ), pT

  fDCosPtNumerator = new TH2D(Form("NumDCosPt%s", suffix),
                              "Numerator; Cos(#DeltaPhi); p_{$}",
                              2*nbins_phi, -1.0, 1.0,
                              nbins_pt, lo_pt, hi_pt);

  fDCosPtDenominator = new TH2D(Form("DenDCosPt%s", suffix),
                                "Denominator; Cos(#DeltaPhi); p_{$}",
                                2*nbins_phi, -1.0, 1.0,
                                nbins_pt, lo_pt, hi_pt);

  // to enable error bar calculation (only enable for numerators)
  fDPhiDEtaNumeratorTrue->Sumw2();
  fDPhiDEtaNumeratorFake->Sumw2();
  fDPhiDEtaColNumerator->Sumw2();
  fDPhiNumeratorTrue->Sumw2();
  fDPhiNumeratorFake->Sumw2();
  fDCosNumeratorTrue->Sumw2();
  fDCosNumeratorFake->Sumw2();
  fDPhiPtNumerator->Sumw2();
  fDCosPtNumerator->Sumw2();
}

AliFemtoModelCorrFctnDEtaDPhiAK
  ::AliFemtoModelCorrFctnDEtaDPhiAK(const AliFemtoModelCorrFctnDEtaDPhiAK& orig)
  : AliFemtoModelCorrFctn(orig)
  , fDPhiDEtaNumeratorTrue(nullptr)
  , fDPhiDEtaNumeratorFake(nullptr)
  , fDPhiDEtaDenominator(nullptr)
  , fDPhiDEtaColNumerator(nullptr)
  , fDPhiDEtaColDenominator(nullptr)
  , fDPhiNumeratorTrue(nullptr)
  , fDPhiNumeratorFake(nullptr)
  , fDPhiDenominator(nullptr)
  , fDCosNumeratorTrue(nullptr)
  , fDCosNumeratorFake(nullptr)
  , fDCosDenominator(nullptr)
  , fDPhiPtNumerator(nullptr)
  , fDPhiPtDenominator(nullptr)
  , fDCosPtNumerator(nullptr)
  , fDCosPtDenominator(nullptr)
{
  fDPhiDEtaNumeratorTrue = new TH2D(*orig.fDPhiDEtaNumeratorTrue);
  fDPhiDEtaNumeratorFake = new TH2D(*orig.fDPhiDEtaNumeratorFake);

  fDPhiDEtaDenominator = new TH2D(*orig.fDPhiDEtaDenominator);
  fDPhiDEtaColNumerator = new TH2D(*orig.fDPhiDEtaColNumerator);
  fDPhiDEtaColDenominator = new TH2D(*orig.fDPhiDEtaColDenominator);

  fDPhiNumeratorTrue = new TH1D(*orig.fDPhiNumeratorTrue);

  fDPhiNumeratorFake = new TH1D(*orig.fDPhiNumeratorFake);
  fDPhiDenominator = new TH1D(*orig.fDPhiDenominator);

  fDCosNumeratorTrue = new TH1D(*orig.fDCosNumeratorTrue);
  fDCosNumeratorFake = new TH1D(*orig.fDCosNumeratorFake);

  fDCosDenominator = new TH1D(*orig.fDCosDenominator);

  fDPhiPtNumerator = new TH2D(*orig.fDPhiPtNumerator);
  fDPhiPtDenominator = new TH2D(*orig.fDPhiPtDenominator);

  fDCosPtNumerator = new TH2D(*orig.fDCosPtNumerator);
  fDCosPtDenominator = new TH2D(*orig.fDCosPtDenominator);
}

AliFemtoModelCorrFctnDEtaDPhiAK::~AliFemtoModelCorrFctnDEtaDPhiAK()
{
  delete fDPhiDEtaNumeratorTrue;
  delete fDPhiDEtaNumeratorFake;
  delete fDPhiDEtaDenominator;
  delete fDPhiDEtaColNumerator;
  delete fDPhiDEtaColDenominator;
  delete fDPhiNumeratorTrue;
  delete fDPhiNumeratorFake;
  delete fDPhiDenominator;
  delete fDCosNumeratorTrue;
  delete fDCosNumeratorFake;
  delete fDCosDenominator;
  delete fDPhiPtNumerator;
  delete fDPhiPtDenominator;
  delete fDCosPtNumerator;
  delete fDCosPtDenominator;
}

AliFemtoModelCorrFctnDEtaDPhiAK&
AliFemtoModelCorrFctnDEtaDPhiAK::operator=(const AliFemtoModelCorrFctnDEtaDPhiAK& rhs)
{
  if (this == &rhs) {
    return *this;
  }

  *fDPhiDEtaNumeratorTrue = *rhs.fDPhiDEtaNumeratorTrue;
  *fDPhiDEtaNumeratorFake = *rhs.fDPhiDEtaNumeratorFake;

  *fDPhiDEtaDenominator = *rhs.fDPhiDEtaDenominator;
  *fDPhiDEtaColNumerator = *rhs.fDPhiDEtaColNumerator;
  *fDPhiDEtaColDenominator = *rhs.fDPhiDEtaColDenominator;

  *fDPhiNumeratorTrue = *rhs.fDPhiNumeratorTrue;
  *fDPhiNumeratorFake = *rhs.fDPhiNumeratorFake;

  *fDPhiDenominator = *rhs.fDPhiDenominator;

  *fDCosNumeratorTrue = *rhs.fDCosNumeratorTrue;
  *fDCosNumeratorFake = *rhs.fDCosNumeratorFake;

  *fDCosDenominator = *rhs.fDCosDenominator;

  *fDPhiPtNumerator = *rhs.fDPhiPtNumerator;
  *fDPhiPtDenominator = *rhs.fDPhiPtDenominator;

  *fDCosPtNumerator = *rhs.fDCosPtNumerator;
  *fDCosPtDenominator = *rhs.fDCosPtDenominator;

  return *this;
}

void
AliFemtoModelCorrFctnDEtaDPhiAK::Finish()
{
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

AliFemtoString
AliFemtoModelCorrFctnDEtaDPhiAK::Report()
{
  // create report
  AliFemtoString report = "TPC Ncls Correlation Function Report:\n";
  report += Form("Number of entries in numerator true:\t%E\n", fDPhiDEtaNumeratorTrue->GetEntries());
  report += Form("Number of entries in numerator fake:\t%E\n", fDPhiDEtaNumeratorFake->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDPhiDEtaDenominator->GetEntries());

  return report;
}

inline
std::tuple<double, double, double, double, double>
calculate_quantities(const AliFemtoPair &pair)
{
  const auto
    &track1 = *pair.Track1()->Track(),
    &track2 = *pair.Track2()->Track();

  const auto
    p1 = track1.P(),
    p2 = track2.P();

  const double
    phi1 = p1.Phi(),
    phi2 = p2.Phi(),
    eta1 = p1.PseudoRapidity(),
    eta2 = p2.PseudoRapidity(),
    deta = eta1 - eta2,

    pt1 = p1.Perp(),
    pt2 = p2.Perp(),

    ptmin = std::min(pt1, pt2),
    cosphi = p1.Dot(p2) / std::sqrt(p1.Mag2() * p2.Mag2()),

    coleta = cosphi > 0.0 ? deta : -(eta1 + eta2);

  double dphi = phi1 - phi2;
  while (__builtin_expect(dphi < -M_PI_2, 0)) {
    dphi += 2 * M_PI;
  }
  while (__builtin_expect(dphi > 3 * M_PI_2, 0)) {
    dphi -= 2 * M_PI;
  }

  return std::make_tuple(dphi, deta, coleta, ptmin, cosphi);
}

void
AliFemtoModelCorrFctnDEtaDPhiAK::AddRealPair(const AliFemtoPair& pair)
{
  double dphi, deta, col_eta, pt, cosphi;
  std::tie(dphi, deta, col_eta, pt, cosphi) = calculate_quantities(pair);

  Double_t weight = fManager->GetWeight(const_cast<AliFemtoPair*>(&pair));
  fDPhiDEtaNumeratorTrue->Fill(dphi, deta, weight);
  fDPhiDEtaColNumerator->Fill(dphi, col_eta);

  fDPhiNumeratorTrue->Fill(dphi, weight);
  fDCosNumeratorTrue->Fill(cosphi, weight);

  fDPhiPtNumerator->Fill(dphi, pt);
  fDCosPtNumerator->Fill(cosphi, pt);
}

void AliFemtoModelCorrFctnDEtaDPhiAK::AddMixedPair(const AliFemtoPair& pair)
{
  double dphi, deta, col_eta, pt, cosphi;
  std::tie(dphi, deta, col_eta, pt, cosphi) = calculate_quantities(pair);

  Double_t weight = fManager->GetWeight(const_cast<AliFemtoPair*>(&pair));
  fDPhiDEtaNumeratorFake->Fill(dphi, deta, weight);
  fDPhiDEtaDenominator->Fill(dphi, deta);
  fDPhiDEtaColDenominator->Fill(dphi, col_eta);

  fDPhiNumeratorFake->Fill(dphi, weight);
  fDCosNumeratorFake->Fill(cosphi, weight);

  fDPhiDenominator->Fill(dphi);
  fDCosDenominator->Fill(cosphi);

  fDPhiPtDenominator->Fill(dphi, pt);
  fDCosPtDenominator->Fill(cosphi, pt);
}

TList* AliFemtoModelCorrFctnDEtaDPhiAK::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDEtaNumeratorTrue);
  tOutputList->Add(fDPhiDEtaNumeratorFake);
  tOutputList->Add(fDPhiDEtaDenominator);
  tOutputList->Add(fDPhiDEtaColNumerator);
  tOutputList->Add(fDPhiDEtaColDenominator);
  tOutputList->Add(fDPhiNumeratorTrue);
  tOutputList->Add(fDPhiNumeratorFake);
  tOutputList->Add(fDPhiDenominator);
  tOutputList->Add(fDCosNumeratorTrue);
  tOutputList->Add(fDCosNumeratorFake);
  tOutputList->Add(fDCosDenominator);
  tOutputList->Add(fDPhiPtNumerator);
  tOutputList->Add(fDPhiPtDenominator);
  tOutputList->Add(fDCosPtNumerator);
  tOutputList->Add(fDCosPtDenominator);

  return tOutputList;
}
