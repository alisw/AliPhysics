///
/// \file AliFemtoModelCorrFctnTrueQ6D.cxx
///

#include "AliFemtoModelCorrFctnTrueQ6D.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelHiddenInfo.h"

#include <TList.h>
#include <TString.h>
#include <TRandom.h>
#include <TAxis.h>

#include <tuple>
#include <iostream>

void
AliFemtoModelCorrFctnTrueQ6D::UpdateQlimits()
{
  const auto aout = fHistogram->GetAxis(3),
             aside = fHistogram->GetAxis(4),
             along = fHistogram->GetAxis(5);

  fQlimits[0] = make_pair(aout->GetXmin(), aout->GetXmax());
  fQlimits[1] = make_pair(aside->GetXmin(), aside->GetXmax());
  fQlimits[2] = make_pair(along->GetXmin(), along->GetXmax());
}

AliFemtoModelCorrFctnTrueQ6D::Builder::Builder(AliFemtoConfigObject cfg)
  : Builder()
{
  cfg.pop_all()
    ("name", title)
    ("qmin", qmin)
    ("qmax", qmax)
    ("nbins", bin_count);
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D()
  : AliFemtoModelCorrFctnTrueQ6D("CF_TrueQ6D")
{
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const HistType &hist, AliFemtoModelManager *mgr)
  : AliFemtoCorrFctn()
  , fManager(mgr)
  , fHistogram(reinterpret_cast<HistType*>(hist.Clone()))
  , fRng(new TRandom())
{
  UpdateQlimits();
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(HistType *&hist, AliFemtoModelManager *mgr)
  : AliFemtoCorrFctn()
  , fManager(mgr)
  , fHistogram(hist)
  , fRng(new TRandom())
{
  // take ownership of the pointer
  hist = nullptr;
  UpdateQlimits();
}


AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const TString &title)
  : AliFemtoModelCorrFctnTrueQ6D(title, 56, 0.14)
{
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const TString &title,
                                                           UInt_t nbins,
                                                           Double_t qmax)
  : AliFemtoModelCorrFctnTrueQ6D(title, nbins, -qmax, qmax)
{
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const TString &title,
                                                           UInt_t nbins,
                                                           Double_t qmin,
                                                           Double_t qmax)
  : AliFemtoCorrFctn()
  , fManager(nullptr)
  , fHistogram(nullptr)
  , fRng(new TRandom())
{
  std::vector<Int_t> nbins_v(6, static_cast<Int_t>(nbins));
  std::vector<double> hist_min(6, qmin), hist_max(6, qmax);
  fHistogram = new HistType(title +"_Histogram", "Histogram; ", 6, nbins_v.data(), hist_min.data(), hist_max.data());
  UpdateQlimits();
}


AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const Builder &params)
  : AliFemtoModelCorrFctnTrueQ6D(params.title.Data(), params.bin_count, params.qmin, params.qmax)
{
  SetManager(params.mc_manager);
  UpdateQlimits();
}


AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const AliFemtoModelCorrFctnTrueQ6D &orig)
  : AliFemtoCorrFctn(orig)
  , fManager(orig.fManager)
  , fHistogram(reinterpret_cast<decltype(fHistogram)>(orig.fHistogram->Clone()))
  , fRng(new TRandom())
{
  UpdateQlimits();
}


AliFemtoModelCorrFctnTrueQ6D::~AliFemtoModelCorrFctnTrueQ6D()
{
  delete fHistogram;
  delete fRng;
}


TList*
AliFemtoModelCorrFctnTrueQ6D::GetOutputList()
{
  TList *result = new TList();
  AppendOutputList(*result);
  return result;
}


TList*
AliFemtoModelCorrFctnTrueQ6D::AppendOutputList(TList &list)
{
  list.Add(fHistogram);
  return &list;
}

//Double_t
static
std::tuple<Double_t, Double_t, Double_t>
Qcms(const AliFemtoLorentzVector &p1, const AliFemtoLorentzVector &p2)
{
  const AliFemtoLorentzVector p = p1 + p2,
                              d = p1 - p2;

  Double_t k1 = p.Perp(),
           k2 = d.x()*p.x() + d.y()*p.y();

  Double_t qout = (k1 == 0) ? 0.0 : k2/k1;

  Double_t qside = (k1 == 0) ? 0.0 : 2.0 * (p2.x()*p1.y() - p1.x()*p2.y())/k1;

  // Double_t beta = p.z()/p.t(),
  //         gamma = 1.0 / TMath::Sqrt((1.0-beta)*(1.0+beta)),
  //         qlong = gamma * (d.z() - beta*d.t());

  Double_t pt = p.t(),
           pz = p.z(),
           dt = d.t(),
           dz = d.z(),
           qlong = (pt*dz - pz*dt) / TMath::Sqrt(pt*pt - pz*pz);

  return std::make_tuple(qout, qside, qlong);
}


void
AliFemtoModelCorrFctnTrueQ6D::AddPair(const AliFemtoParticle &particle1, const AliFemtoParticle &particle2)
{
  // Fill reconstructed histogram with "standard" particle momentum
  Double_t q_out, q_side, q_long;
  std::tie(q_out, q_side, q_long) = Qcms(particle1.FourMomentum(), particle2.FourMomentum());

  auto out_of_bounds = [] (double q, const std::pair<double,double> limit) {
    return q < limit.first || limit.second < q;
  };

  if (out_of_bounds(q_out, fQlimits[0]) ||
      out_of_bounds(q_side, fQlimits[1]) ||
      out_of_bounds(q_long, fQlimits[2])) {
    return;
  }

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

  const AliFemtoThreeVector &true_momentum1 = *info1->GetTrueMomentum(),
                            &true_momentum2 = *info2->GetTrueMomentum();

  const Double_t e1 = sqrt(mass1 * mass1 + true_momentum1.Mag2()),
                 e2 = sqrt(mass2 * mass2 + true_momentum2.Mag2());

  const AliFemtoLorentzVector p1(e1, true_momentum1),
                              p2(e2, true_momentum2);

  // Fill generated-momentum histogram with "true" particle momentum
  Double_t true_q_out, true_q_side, true_q_long;
  std::tie(true_q_out, true_q_side, true_q_long) = Qcms(p1, p2);

  Double_t q[6] = {true_q_out, true_q_side, true_q_long, q_out, q_side, q_long};
  fHistogram->Fill(q);
}

void
AliFemtoModelCorrFctnTrueQ6D::AddRealPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // randomize to avoid ordering biases
  if (fRng->Uniform() >= 0.5) {
    std::swap(p1, p2);
  }
  AddPair(*p1, *p2);
}

void
AliFemtoModelCorrFctnTrueQ6D::AddMixedPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // randomize to avoid ordering biases
  if (fRng->Uniform() >= 0.5) {
    std::swap(p1, p2);
  }

  AddPair(*p1, *p2);
}

AliFemtoString
AliFemtoModelCorrFctnTrueQ6D::Report()
{
  TString report;
  return AliFemtoString(report.Data());
}
