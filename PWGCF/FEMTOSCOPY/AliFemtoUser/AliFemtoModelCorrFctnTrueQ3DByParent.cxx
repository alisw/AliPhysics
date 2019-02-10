///
/// \file AliFemtoModelCorrFctnTrueQ3DByParentByParent.cxx
///

#include "AliFemtoModelCorrFctnTrueQ3DByParent.h"

#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelHiddenInfo.h"


#include <THnSparse.h>
#include <TAxis.h>
#include <TList.h>
#include <TString.h>


#include <tuple>
#include <iostream>
#include <array>

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent():
  AliFemtoModelCorrFctnTrueQ3DByParent("CF_TrueQ3DBP")
{
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const char *title)
  : AliFemtoModelCorrFctnTrueQ3DByParent(title, 57, 0.1425)
{
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const char *title, UInt_t nbins, Double_t qmax)
  : AliFemtoModelCorrFctnTrueQ3DByParent(title, nbins, -qmax, qmax)
{
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(Int_t nbins, Double_t qmax)
  : AliFemtoModelCorrFctnTrueQ3DByParent(nbins, -abs(qmax), abs(qmax))
{
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(Int_t nbins, Double_t qmin, Double_t qmax)
  : AliFemtoModelCorrFctnTrueQ3DByParent("", nbins, -qmax, qmax)
{
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const char *prefix, UInt_t nbins, Double_t qmin, Double_t qmax)
  : AliFemtoCorrFctn()
  , fManager(nullptr)
  , fNumeratorGenerated(nullptr)
  , fNumeratorReconstructed(nullptr)
  , fDenominatorGenerated(nullptr)
  , fDenominatorReconstructed(nullptr)
{
  const TString axis_titles = "; parent_code; q_{out}; q_{side}; q_{long}";

  const UInt_t N_AXES = 4;

  const float
    par_nbins = 4000,
    parent_max = par_nbins + 0.5,
    parent_min = -parent_max;

  Int_t inbins = nbins;
  std::array<Int_t, 4> nbins_v = {static_cast<Int_t>(par_nbins * 2 + 1), inbins, inbins, inbins};
  std::array<Double_t, 4> min_v = {parent_min, qmin, qmin, qmin};
  std::array<Double_t, 4> max_v = {parent_max, qmax, qmax, qmax};

  auto build_sparse_hist = [&] (const TString &name, const char *title)
    {
      return new THnSparseF(prefix + name,
                            title + axis_titles,
                            N_AXES, nbins_v.data(), min_v.data(), max_v.data());
    };

  fNumeratorGenerated = build_sparse_hist("NumGen", "Numerator (MC-Generated Momentum)");
  fNumeratorReconstructed = build_sparse_hist("NumRec", "Numerator (Reconstructed Momentum)");

  fDenominatorGenerated = build_sparse_hist("DenGen", "Denominator (MC-Generated Momentum)");
  fDenominatorReconstructed = build_sparse_hist("DenRec", "Denominator (Reconstructed Momentum)");
}


AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const AliFemtoModelCorrFctnTrueQ3DByParent::Parameters &params):
  AliFemtoModelCorrFctnTrueQ3DByParent(params.prefix.Data(), params.bin_count, params.qmin, params.qmax)
{
  SetManager(params.mc_manager);
}


AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const AliFemtoModelCorrFctnTrueQ3DByParent& orig)
  : AliFemtoCorrFctn(orig)
  , fManager(orig.fManager)
  , fNumeratorGenerated(nullptr)
  , fNumeratorReconstructed(nullptr)
  , fDenominatorGenerated(nullptr)
  , fDenominatorReconstructed(nullptr)
{
  fNumeratorGenerated = static_cast<THnSparseF*>(orig.fNumeratorGenerated->Clone());
  fNumeratorReconstructed = static_cast<THnSparseF*>(orig.fNumeratorReconstructed->Clone());
  fDenominatorGenerated = static_cast<THnSparseF*>(orig.fDenominatorGenerated->Clone());
  fDenominatorReconstructed = static_cast<THnSparseF*>(orig.fDenominatorReconstructed->Clone());
}


AliFemtoModelCorrFctnTrueQ3DByParent&
AliFemtoModelCorrFctnTrueQ3DByParent::operator=(const AliFemtoModelCorrFctnTrueQ3DByParent &rhs)
{
  if (this == &rhs) {
    return *this;
  }

  AliFemtoCorrFctn::operator=(rhs);

  delete fNumeratorGenerated;
  fNumeratorGenerated = static_cast<THnSparseF*>(rhs.fNumeratorGenerated->Clone());

  delete fNumeratorReconstructed;
  fNumeratorReconstructed = static_cast<THnSparseF*>(rhs.fNumeratorReconstructed->Clone());

  delete fDenominatorGenerated;
  fDenominatorGenerated = static_cast<THnSparseF*>(rhs.fDenominatorGenerated->Clone());

  delete fDenominatorReconstructed;
  fDenominatorReconstructed = static_cast<THnSparseF*>(rhs.fDenominatorReconstructed->Clone());

  return *this;
}

AliFemtoModelCorrFctnTrueQ3DByParent::~AliFemtoModelCorrFctnTrueQ3DByParent()
{
  delete fNumeratorGenerated;
  delete fNumeratorReconstructed;
  delete fDenominatorGenerated;
  delete fDenominatorReconstructed;
}


TList*
AliFemtoModelCorrFctnTrueQ3DByParent::GetOutputList()
{
  TList *result = new TList();
  AppendOutputList(*result);
  return result;
}


TList*
AliFemtoModelCorrFctnTrueQ3DByParent::AppendOutputList(TList &list)
{
  AddOutputObjectsTo(list);

  return &list;
}

void
AliFemtoModelCorrFctnTrueQ3DByParent::AddOutputObjectsTo(TCollection &dest)
{
  dest.Add(fNumeratorGenerated);
  dest.Add(fNumeratorReconstructed);
  dest.Add(fDenominatorGenerated);
  dest.Add(fDenominatorReconstructed);
}

static
std::tuple<Double_t, Double_t, Double_t>
Qcms(const AliFemtoLorentzVector &p1, const AliFemtoLorentzVector &p2)
{
  const AliFemtoLorentzVector p = p1 + p2,
                              d = p1 - p2;

  const Double_t
    k1 = p.Perp(),
    k2 = d.x()*p.x() + d.y()*p.y(),
    beta = p.z()/p.t(),
    gamma = 1.0 / TMath::Sqrt((1.0-beta)*(1.0+beta)),

    // relative momentum out component in lab frame
    qout = (k1 == 0) ? 0.0 : k2 / k1,
    qside = (k1 == 0) ? 0.0 : 2.0 * (p2.x()*p1.y() - p1.x()*p2.y()) / k1,
    qlong = gamma * (d.z() - beta*d.t());
    // qlong = (p.t()*d.z() - p.z()*d.t()) / TMath::Sqrt(p.t()*p.t() - p.z()*p.z());

  return std::make_tuple(qout, qside, qlong);
}


static void
AddPair(const AliFemtoParticle &particle1,
        const AliFemtoParticle &particle2,
        const Double_t weight,
        THnSparseF *gen_hist,
        THnSparseF *rec_hist)
{
  // Get generated momentum from hidden info
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(particle1.HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(particle2.HiddenInfo());

  if (info1 == nullptr || info2 == nullptr) {
    return;
  }

  const Float_t mass1 = info1->GetMass(),
                mass2 = info2->GetMass();

  // block all zero-mass particles from the correlation function
  if (mass1 == 0.0 || mass2 == 0.0) {
    return;
  }

  const Int_t parent1 = info1->GetMotherPdgCode(),
              parent2 = info2->GetMotherPdgCode();

  // If the pair is made of two unique parents, store in 1, otherwise
  // use parent PDG code as the first bin.
  const double pbin = (parent1 == 0 || parent2 == 0 || parent1 == parent2)
                    ? parent1 | parent2
                    : 1;

  Double_t q_out, q_side, q_long;

  auto fill_if_within_bounds = [] (THnSparseF &hist,
                                   const std::array<double, 4> &data,
                                   const double weight)
    {
      for (int i=0; i<4; ++i) {
        const auto *axis = hist.GetAxis(i);
        const double x = data[i];

        if (x < axis->GetXmin() || axis->GetXmax() <= x) {
          return;
        }
      }

      hist.Fill(data.data(), weight);
    };

  std::tie(q_out, q_side, q_long) = Qcms(particle1.FourMomentum(), particle2.FourMomentum());
  fill_if_within_bounds(*rec_hist, {pbin, q_out, q_side, q_long}, weight);

  // Fill generated-momentum histogram with "true" particle momentum
  const AliFemtoThreeVector &true_momentum1 = *info1->GetTrueMomentum(),
                            &true_momentum2 = *info2->GetTrueMomentum();

  const Double_t e1 = std::sqrt(mass1 * mass1 + true_momentum1.Mag2()),
                 e2 = std::sqrt(mass2 * mass2 + true_momentum2.Mag2());

  std::tie(q_out, q_side, q_long) = Qcms({e1, true_momentum1}, {e2, true_momentum2});
  fill_if_within_bounds(*gen_hist, {pbin, q_out, q_side, q_long}, weight);
}


void
AliFemtoModelCorrFctnTrueQ3DByParent::AddRealPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  AddPair(*p1, *p2, fManager->GetWeight(pair), fNumeratorGenerated, fNumeratorReconstructed);
}


void
AliFemtoModelCorrFctnTrueQ3DByParent::AddMixedPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  AddPair(*p1, *p2, 1.0, fDenominatorGenerated, fDenominatorReconstructed);
}


AliFemtoString
AliFemtoModelCorrFctnTrueQ3DByParent::Report()
{
  TString report;
  return AliFemtoString(report.Data());
}
