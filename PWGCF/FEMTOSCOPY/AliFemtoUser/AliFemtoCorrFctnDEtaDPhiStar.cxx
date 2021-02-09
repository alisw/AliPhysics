///
/// \file AliFemtoCorrFctnDEtaDPhiStar.cxx
///

#include "AliFemtoCorrFctnDEtaDPhiStar.h"

#include "AliFemtoPairCutDetaDphi.h"


AliFemtoCorrFctnDEtaDPhiStar
  ::AliFemtoCorrFctnDEtaDPhiStar()
  : AliFemtoCorrFctnDEtaDPhiStar(Params())
{
}

AliFemtoCorrFctnDEtaDPhiStar
  ::AliFemtoCorrFctnDEtaDPhiStar(const char *suffix,
                                 double radius)
  : AliFemtoCorrFctnDEtaDPhiStar(Params().Suffix(suffix).Radius(radius))
{
}

AliFemtoCorrFctnDEtaDPhiStar
  ::AliFemtoCorrFctnDEtaDPhiStar(const Params &p)
  : AliFemtoCorrFctn()
  , fDPhiDEtaNumerator(nullptr)
  , fDPhiDEtaDenominator(nullptr)
  , fRadius(p.radius)
  , fMagneticField(0)
{
  fDPhiDEtaNumerator = new TH2F("Num" + p.suffix,
                                Form("#Delta#eta & #Delta#phi* Numerator - (Radius = %g);"
                                     "#Delta#eta;"
                                     "#Delta#phi*;", fRadius),
                                p.nbins_eta, -p.eta_max, p.eta_max,
                                p.nbins_phi, -p.phi_max, p.phi_max);

  fDPhiDEtaDenominator = new TH2F("Den" + p.suffix,
                                Form("#Delta#eta & #Delta#phi* Denominator - (Radius = %g);"
                                     "#Delta#eta;"
                                     "#Delta#phi*;", fRadius),
                                p.nbins_eta, -p.eta_max, p.eta_max,
                                p.nbins_phi, -p.phi_max, p.phi_max);
}

AliFemtoCorrFctnDEtaDPhiStar
  ::AliFemtoCorrFctnDEtaDPhiStar(const AliFemtoCorrFctnDEtaDPhiStar &orig)
  : AliFemtoCorrFctn(orig)
  , fDPhiDEtaNumerator(static_cast<TH2F*>(orig.fDPhiDEtaNumerator->Clone()))
  , fDPhiDEtaDenominator(static_cast<TH2F*>(orig.fDPhiDEtaDenominator->Clone()))
  , fRadius(orig.fRadius)
  , fMagneticField(0)
{
  fDPhiDEtaNumerator->Reset();
  fDPhiDEtaDenominator->Reset();
}

AliFemtoCorrFctnDEtaDPhiStar&
AliFemtoCorrFctnDEtaDPhiStar::operator=(const AliFemtoCorrFctnDEtaDPhiStar &rhs)
{
  if (this != &rhs) {
    AliFemtoCorrFctn::operator=(rhs);

    *fDPhiDEtaNumerator = *rhs.fDPhiDEtaNumerator;
    *fDPhiDEtaDenominator = *rhs.fDPhiDEtaDenominator;
    fRadius = rhs.fRadius;

    // fDPhiDEtaNumerator->Reset();
    // fDPhiDEtaDenominator->Reset();
  }

  return *this;
}

AliFemtoCorrFctnDEtaDPhiStar::~AliFemtoCorrFctnDEtaDPhiStar()
{
}

void
AliFemtoCorrFctnDEtaDPhiStar::AddPair(TH2 &hist, const AliFemtoPair &pair)
{
  const AliFemtoTrack
    &track1 = *pair.Track1()->Track(),
    &track2 = *pair.Track2()->Track();

  const auto &p1 = track1.P(),
             &p2 = track2.P();

  const double
    phistar = AliFemtoPairCutDetaDphi::CalculateDPhiStar(
                  p1, track1.Charge(),
                  p2, track2.Charge(),
                  fRadius,
                  fMagneticField),
    eta = AliFemtoPairCutDetaDphi::CalculateDEta(p1, p2);

  hist.Fill(eta, phistar);
}

void
AliFemtoCorrFctnDEtaDPhiStar::AddRealPair(AliFemtoPair *pair)
{
  AddPair(*fDPhiDEtaNumerator, *pair);
}

void
AliFemtoCorrFctnDEtaDPhiStar::AddMixedPair(AliFemtoPair *pair)
{
  AddPair(*fDPhiDEtaDenominator, *pair);
}


TList*
AliFemtoCorrFctnDEtaDPhiStar::GetOutputList()
{
  TList *list = new TList();
  list->Add(fDPhiDEtaNumerator);
  list->Add(fDPhiDEtaDenominator);
  return list;
}
