///
/// \file AliFemtoModelCorrFctnTrueQ3DByParentByParent.cxx
///

#include "AliFemtoModelCorrFctnTrueQ3DByParent.h"

#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelHiddenInfo.h"

#include <TRandom.h>

#include <THnSparse.h>
#include <TList.h>
#include <TString.h>
#include <TRandom.h>


#include <tuple>
#include <iostream>

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent():
  AliFemtoModelCorrFctnTrueQ3DByParent("CF_TrueQ3DBP")
{
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const char *title)
  : AliFemtoModelCorrFctnTrueQ3DByParent(title, 56, 0.14)
{
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const char *title, UInt_t nbins, Double_t qmax)
  : AliFemtoModelCorrFctnTrueQ3DByParent(title, nbins, -qmax, qmax)
{
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const char *name, UInt_t nbins, Double_t qmin, Double_t qmax)
  : AliFemtoCorrFctn()
  , fManager(nullptr)
  , fNumeratorGenerated(nullptr)
  , fNumeratorReconstructed(nullptr)
  , fDenominatorGenerated(nullptr)
  , fDenominatorReconstructed(nullptr)
  , fRng(new TRandom())
{
  /*
  fNumeratorGenerated = new THnSparseF(
    TString::Format("%s_NumGen", title), "Numerator (MC-Generated Momentum)",
                                 nbins, qmin, qmax,
                                 nbins, qmin, qmax,
                                 nbins, qmin, qmax);
  fNumeratorReconstructed = new TH3D(TString::Format("%s_NumRec", title), "Numerator (Reconstructed Momentum)",
                                 nbins, qmin, qmax,
                                 nbins, qmin, qmax,
                                 nbins, qmin, qmax);
  fDenominatorGenerated = new TH3D(TString::Format("%s_DenGen", title), "Denominator (MC-Generated Momentum)",
                                   nbins, qmin, qmax,
                                   nbins, qmin, qmax,
                                   nbins, qmin, qmax);
  fDenominatorReconstructed = new TH3D(TString::Format("%s_DenRec", title), "Denominator (Reconstructed Momentum)",
                                   nbins, qmin, qmax,
                                   nbins, qmin, qmax,
                                   nbins, qmin, qmax);
                                   */
}

AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(Int_t nbins, Double_t qmin, Double_t qmax)
  : AliFemtoCorrFctn()
  , fManager(nullptr)
  , fNumeratorGenerated(nullptr)
  , fNumeratorReconstructed(nullptr)
  , fDenominatorGenerated(nullptr)
  , fDenominatorReconstructed(nullptr)
  , fRng(new TRandom())
{
  const TString axis_titles = "; parent_code; q_{out}; q_{side}; q_{long}";

  const UInt_t N_AXES = 4;

  const float
    par_nbins = 4000,
    parent_max = par_nbins + 0.5,
    parent_min = -parent_max;

  std::array<Int_t, 4> nbins_v = {static_cast<int>(par_nbins * 2 + 1), nbins, nbins, nbins};
  std::array<Double_t, 4> min_v = {parent_min, qmin, qmin, qmin};
  std::array<Double_t, 4> max_v = {parent_max, qmax, qmax, qmax};

  fNumeratorGenerated = new THnSparseF(
    "NumGen", "Numerator (MC-Generated Momentum)" + axis_titles,
    N_AXES, nbins_v.data(), min_v.data(), max_v.data());

  fNumeratorReconstructed = new THnSparseF(
    "NumRec", "Numerator (Reconstructed Momentum)" + axis_titles,
    N_AXES, nbins_v.data(), min_v.data(), max_v.data());

  fDenominatorGenerated = new THnSparseF(
    "DenGen", "Denominator (MC-Generated Momentum)" + axis_titles,
    N_AXES, nbins_v.data(), min_v.data(), max_v.data());

  fDenominatorReconstructed = new THnSparseF(
    "%s_DenRec", "Denominator (Reconstructed Momentum)" + axis_titles,
    N_AXES, nbins_v.data(), min_v.data(), max_v.data());
}


AliFemtoModelCorrFctnTrueQ3DByParent::AliFemtoModelCorrFctnTrueQ3DByParent(const AliFemtoModelCorrFctnTrueQ3DByParent::Parameters &params):
  AliFemtoModelCorrFctnTrueQ3DByParent(params.title.Data(), params.bin_count, params.qmin, params.qmax)
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
  , fRng(new TRandom())
{
  fNumeratorGenerated = static_cast<THnSparseF*>(orig.fNumeratorGenerated->Clone());
  // fNumeratorReconstructed = new THnSparseF(*orig.fNumeratorReconstructed);
  // fDenominatorGenerated = new THnSparseF(*orig.fDenominatorGenerated);
  // fDenominatorReconstructed = new THnSparseF(*orig.fDenominatorReconstructed);
}


AliFemtoModelCorrFctnTrueQ3DByParent&
AliFemtoModelCorrFctnTrueQ3DByParent::operator=(const AliFemtoModelCorrFctnTrueQ3DByParent&)
{
  return *this;
}

AliFemtoModelCorrFctnTrueQ3DByParent::~AliFemtoModelCorrFctnTrueQ3DByParent()
{
  delete fNumeratorGenerated;
  delete fNumeratorReconstructed;
  delete fDenominatorGenerated;
  delete fDenominatorReconstructed;
  delete fRng;
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
  list.Add(fNumeratorGenerated);
  list.Add(fNumeratorReconstructed);
  list.Add(fDenominatorGenerated);
  list.Add(fDenominatorReconstructed);

  return &list;
}

//Double_t
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

  {
    std::tie(q_out, q_side, q_long) = Qcms(particle1.FourMomentum(), particle2.FourMomentum());
    const Double_t data[] = {pbin, q_out, q_side, q_long};
    rec_hist->Fill(data, weight);
  }

  // Fill generated-momentum histogram with "true" particle momentum
  {
    const AliFemtoThreeVector &true_momentum1 = *info1->GetTrueMomentum(),
                              &true_momentum2 = *info2->GetTrueMomentum();

    const Double_t e1 = std::sqrt(mass1 * mass1 + true_momentum1.Mag2()),
                   e2 = std::sqrt(mass2 * mass2 + true_momentum2.Mag2());

    std::tie(q_out, q_side, q_long) = Qcms({e1, true_momentum1}, {e2, true_momentum2});
    const Double_t data[] = {pbin, q_out, q_side, q_long};
    gen_hist->Fill(data, weight);
  }
}

void
AliFemtoModelCorrFctnTrueQ3DByParent::AddRealPair(AliFemtoPair *pair)
{
  Double_t weight = fManager->GetWeight(pair);

  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // randomize to avoid ordering biases
  if (fRng->Uniform() >= 0.5) {
    std::swap(p1, p2);
  }
  AddPair(*p1, *p2, weight, fNumeratorGenerated, fNumeratorReconstructed);
}

void
AliFemtoModelCorrFctnTrueQ3DByParent::AddMixedPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // randomize to avoid ordering biases
  if (fRng->Uniform() >= 0.5) {
    std::swap(p1, p2);
  }

  AddPair(*p1, *p2, 1.0, fDenominatorGenerated, fDenominatorReconstructed);
}


AliFemtoString
AliFemtoModelCorrFctnTrueQ3DByParent::Report()
{
  TString report;
  return AliFemtoString(report.Data());
}
